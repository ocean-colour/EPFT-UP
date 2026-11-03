"""ELBO training, augmentation and cross-validation folds (cf_design.md §6.3, §7, §8.1).

Loss (negative ELBO): ``L = (1/N) Σ_i NLL_i + (β/N)·KL(q(w)‖p(w))`` with N the
number of *independent labels* in the training fold (not augmented spectra),
one weight sample per step (local reparameterisation).
"""

import copy
import time

import numpy as np
import tensorflow as tf
import keras

from epft_up.claude_fun import algorithm as alg
from epft_up.claude_fun import pipeline as pl

DEFAULT_TRAIN = {
    'lr': 1e-3,
    'lr_schedule': 'cosine',               # cosine decay over max_epochs
    'batch_size': 16,
    'max_epochs': 3000,
    'patience': 200,                       # epochs without inner-fold NLL improvement
    'eval_every': 5,                       # epochs between inner-fold NLL evaluations
    'inner_val_frac': 0.2,                 # inner split of each training fold (early stopping)
    'inner_val_draws': 16,                 # weight draws for the inner-fold predictive NLL
    'beta': 1.0,
    'kl_warmup': False,                    # optional linear warm-up of β ...
    'kl_warmup_frac': 0.1,                 # ... over the first 10 % of epochs
    'seed': 0,
    'ensemble_seeds': None,                # e.g. (0, 1, 2, 3, 4) for the optional ensemble
    'augment': {
        'replicate': True,                 # draw a member 1-min spectrum from the ±15 min bin
        'noise': True,                     # Rrs += ε·σ_Rrs
        'degrade': True,                   # band dropout + Gaussian SRF resampling
        'p_resample': 0.5,                 # probability of SRF convolution + regular grid
        'p_dropout': 0.5,                  # probability of band dropout
        'keep_frac': (0.3, 1.0),
        'fwhm_nm': (3.5, 15.0),
        'spacing_nm': (2.0, 10.0),
        'min_tokens': 5,
    },
    'cv_lat_blocks': 5,
    'sst_permutations': 5,                 # repeats of the SST permutation test
    'moana': {'n_pc': {'pro': 5, 'syn': 5, 'peuk': 5},   # retrained-MOANA baseline (§8.3)
              'pro_space': 'linear',       # MOANA fits Pro in linear cells mL^-1
              'use_sst': True,
              'log_clip': 1.0},            # floor (cells mL^-1) before log10 of linear Pro
    'published_moana': {'module': 'epft_up.moana.algorithm', 'func': 'run_moana',
                        'lut_module': 'epft_up.moana.io', 'lut_func': 'load_luts',
                        'keys': {'pro': 'pro', 'syn': 'syn', 'peuk': 'peuk'}},
    'verbose_every': 250,
}


def config(**overrides):
    """Deep copy of :data:`DEFAULT_TRAIN` with overrides (``augment`` and
    ``moana`` are merged one level deep). Outputs: dict."""
    cfg = copy.deepcopy(DEFAULT_TRAIN)
    for k, v in overrides.items():
        if k not in cfg:
            raise KeyError(f"unknown train key {k!r}")
        if isinstance(cfg[k], dict) and isinstance(v, dict):
            cfg[k].update(v)
        else:
            cfg[k] = v
    return cfg


# ==========================================================================
# Data
# ==========================================================================

def prepare(ds, pipeline_cfg=None):
    """Training/evaluation view of a pipeline dataset.

    Inputs: ds (training matrix or transfer set); pipeline_cfg (σ floor).
    Outputs: dict ``ds``, ``spec`` (median spectra set), ``labels``,
    ``wave``, ``rrs``, ``sigma``, ``sigma_bin`` (or None), ``members`` (list
    of (m_i, W) arrays, or None), ``lat``, ``floor`` (f_rel, σ_abs).
    """
    pc = pipeline_cfg or pl.DEFAULT_PIPELINE
    members = None
    if 'member_rrs' in ds:
        M = ds['member_rrs'].values
        members = [M[i][np.isfinite(M[i]).any(axis=1)] for i in range(M.shape[0])]
    return {'ds': ds, 'spec': alg.spectra_from_dataset(ds), 'labels': alg.labels_from_dataset(ds),
            'wave': ds['wavelength'].values, 'rrs': ds['rrs'].values, 'sigma': ds['sigma'].values,
            'sigma_bin': ds['sigma_bin'].values if 'sigma_bin' in ds else None,
            'members': members, 'lat': ds['lat'].values, 'n': ds.sizes['sample'],
            'floor': (pc['sigma_f_rel'], pc['sigma_abs'])}


def sub_labels(labels, rows):
    """Rows of a labels dict. Inputs: labels, rows. Outputs: dict."""
    return {k: v[np.asarray(rows)] for k, v in labels.items()}


def _floor(r, floor):
    return np.sqrt((floor[0] * r) ** 2 + floor[1] ** 2)


def gaussian_srf_resample(wave, rrs, sigma, grid, fwhm):
    """Convolve with a Gaussian spectral response and sample on ``grid``;
    σ propagated as sqrt(Σ_j w_j² σ_j²) (band-independent errors).

    Inputs: wave, rrs, sigma (1-D, finite), grid (target centres), fwhm [nm].
    Outputs: (grid_kept, rrs_out, sigma_out) — only centres whose ±FWHM lies
    inside the input coverage are kept.
    """
    sd = fwhm / (2 * np.sqrt(2 * np.log(2)))
    dl = np.gradient(wave) if wave.size > 1 else np.ones(1)
    W = np.exp(-0.5 * ((grid[:, None] - wave[None, :]) / sd) ** 2) * dl[None, :]
    keep = (grid - fwhm >= wave[0]) & (grid + fwhm <= wave[-1]) & (W.sum(1) > 0)
    W = W[keep] / W[keep].sum(1, keepdims=True)
    return grid[keep], W @ rrs, np.sqrt((W ** 2) @ sigma ** 2)


def augment_spectrum(rng, data, i, aug, wave_range):
    """One augmented draw of sample i (§7.2): replicate sampling, noise,
    spectral degradation (Gaussian SRF on a random regular grid; band
    dropout), with the σ floor re-applied after degradation.

    Inputs: rng (numpy Generator), data (prepare), i, aug (augment cfg),
    wave_range (model range, nm).
    Outputs: (wave, rrs, sigma) 1-D arrays.
    """
    w = data['wave']
    fl = data['floor']
    if aug.get('replicate') and data['members'] is not None and len(data['members'][i]):
        mem = data['members'][i]
        r = mem[rng.integers(len(mem))].copy()
        sb = data['sigma_bin'][i] if data['sigma_bin'] is not None else 0.0
        s = np.fmax(np.nan_to_num(sb), _floor(r, fl))
    else:
        r, s = data['rrs'][i].copy(), data['sigma'][i].copy()
    ok = np.isfinite(r) & np.isfinite(s)
    w, r, s = w[ok], r[ok], s[ok]
    if aug.get('noise'):
        r = r + rng.standard_normal(r.size) * s
    if aug.get('degrade'):
        lo, hi = wave_range
        if rng.random() < aug['p_resample']:
            fwhm = rng.uniform(*aug['fwhm_nm'])
            dl = rng.uniform(*aug['spacing_nm'])
            g = np.arange(lo + rng.uniform(0, dl), hi + 1e-9, dl)
            g2, r2, s2 = gaussian_srf_resample(w, r, s, g, fwhm)
            if ((g2 >= lo) & (g2 <= hi)).sum() >= aug['min_tokens']:
                w, r, s = g2, r2, np.fmax(s2, _floor(r2, fl))
        if rng.random() < aug['p_dropout']:
            inr = np.flatnonzero((w >= lo) & (w <= hi))
            k = max(aug['min_tokens'], int(round(rng.uniform(*aug['keep_frac']) * inr.size)))
            if k < inr.size:
                keep = np.sort(rng.choice(inr, k, replace=False))
                w, r, s = w[keep], r[keep], s[keep]
    return w, r, s


def augmented_spec(rng, data, rows, aug, wave_range):
    """Spectra set of augmented draws for ``rows``. Outputs: spectra set."""
    out = {'wave': [], 'rrs': [], 'sigma': [], 'sst': data['spec']['sst'][np.asarray(rows)]}
    for i in rows:
        w, r, s = augment_spectrum(rng, data, i, aug, wave_range)
        out['wave'].append(w)
        out['rrs'].append(r)
        out['sigma'].append(s)
    return out


# ==========================================================================
# Training
# ==========================================================================

def _make_step(model, opt, n_labels):
    """Compiled training step for one batch (one weight sample, LRT)."""
    vars_ = alg.trainable_variables(model)
    counterbin = model['cfg']['encoder'] == 'counterbin'

    def loss_fn(x, mask, sst, y, obs, cens, cvar, beta):
        if counterbin:
            h, _ = model['encoder'](x, training=True)
        else:
            h, _ = model['encoder'](x, mask, training=True)
        mu, sig, L = alg.head_forward(model, h, sst, mode='lrt')
        nll = alg.joint_nll_tf({'y': y, 'obs': obs, 'cens': cens, 'cvar': cvar}, mu, sig, L)
        mean_nll = tf.reduce_mean(nll)
        return mean_nll + beta * tf.cast(alg.kl(model), tf.float64) / n_labels, mean_nll

    @tf.function(reduce_retracing=True)
    def step(x, mask, sst, y, obs, cens, cvar, beta):
        with tf.GradientTape() as tape:
            loss, nll = loss_fn(x, mask, sst, y, obs, cens, cvar, beta)
        g = tape.gradient(loss, vars_)
        opt.apply_gradients([(gi, v) for gi, v in zip(g, vars_) if gi is not None])
        return loss, nll
    return step


def predictive_nll_std(model, batch, lab, n_draws=16, seed=999):
    """Mean NLL of the moment-matched predictive in standardised units
    (early-stopping criterion). Inputs: model, batch, lab (standardised),
    n_draws. Outputs: float."""
    h, _ = alg.encode(model, batch)
    mode = 'sample' if model['cfg']['bayesian_head'] else 'mean'
    mus, covs = [], []
    for t in range(n_draws if mode == 'sample' else 1):
        mu, sig, L = alg.head_forward(model, h, batch['sst'], mode=mode, seed=[seed, t])
        Ls = sig.numpy()[:, :, None] * L.numpy()
        mus.append(mu.numpy())
        covs.append(Ls @ np.transpose(Ls, (0, 2, 1)))
    mus, covs = np.stack(mus), np.stack(covs)
    m = mus.mean(0)
    d = mus - m
    cov = covs.mean(0) + np.einsum('tni,tnj->nij', d, d) / max(len(mus) - 1, 1)
    nll = alg.joint_nll_np(lab['y'], lab['obs'] > 0, lab['cens'] > 0, m, cov, lab['cvar'])
    return float(np.nanmean(nll))


def train_cf(data, rows=None, model_cfg=None, train_cfg=None, seed=None, verbose=False):
    """Train one Claude Fun model on the given rows (§6.3, §7).

    An inner split (``inner_val_frac``) of the training rows is held out for
    early stopping on its predictive NLL; the best weights are restored.
    All standardisation statistics come from the inner training rows.

    Inputs: data (prepare), rows (training rows; default all), model_cfg,
    train_cfg, seed (default train_cfg['seed']), verbose.
    Outputs: model dict (algorithm.build_model) with ``history`` and
    ``meta`` (train/val rows, N, best epoch, seed, seconds).
    """
    mc = copy.deepcopy(model_cfg or alg.DEFAULT_MODEL)
    tc = copy.deepcopy(train_cfg or DEFAULT_TRAIN)
    seed = tc['seed'] if seed is None else seed
    keras.utils.set_random_seed(int(seed))
    rng = np.random.default_rng(int(seed))
    rows = np.arange(data['n']) if rows is None else np.asarray(rows)
    perm = rng.permutation(rows)
    n_val = int(round(tc['inner_val_frac'] * rows.size)) if rows.size >= 5 else 0
    val, tr = np.sort(perm[:n_val]), np.sort(perm[n_val:])
    spec, lab = data['spec'], data['labels']
    stats = alg.fit_stats(alg.subset(spec, tr), sub_labels(lab, tr), mc)
    model = alg.build_model(mc, stats)
    N = int(tr.size)
    bs = int(tc['batch_size'])
    spe = int(np.ceil(N / bs))
    if tc['lr_schedule'] == 'cosine':
        lr = keras.optimizers.schedules.CosineDecay(tc['lr'], decay_steps=tc['max_epochs'] * spe)
    else:
        lr = tc['lr']
    opt = keras.optimizers.Adam(lr)
    step = _make_step(model, opt, float(N))
    vb = vl = None
    if n_val:
        vb = alg.make_batch(alg.subset(spec, val), stats, mc)
        vl = alg.standardise_labels(sub_labels(lab, val), stats, mc['counting_noise'])
    vars_ = alg.trainable_variables(model)
    best, best_w, best_epoch, wait = np.inf, None, -1, 0
    hist = {'epoch': [], 'loss': [], 'val_nll': []}
    t0 = time.time()
    aug = tc['augment']
    for ep in range(tc['max_epochs']):
        beta = tc['beta']
        if tc['kl_warmup']:
            beta *= min(1.0, (ep + 1) / max(1.0, tc['kl_warmup_frac'] * tc['max_epochs']))
        order = rng.permutation(tr)
        losses = []
        for b in range(spe):
            idx = order[b * bs:(b + 1) * bs]
            sb = augmented_spec(rng, data, idx, aug, mc['wave_range_nm'])
            batch = alg.make_batch(sb, stats, mc)
            ls = alg.standardise_labels(sub_labels(lab, idx), stats, mc['counting_noise'])
            mask = batch['mask'] if batch['mask'] is not None else np.ones((len(idx), 1), bool)
            loss, _ = step(tf.constant(batch['x']), tf.constant(mask), tf.constant(batch['sst']),
                           *(tf.constant(ls[k]) for k in ('y', 'obs', 'cens', 'cvar')),
                           tf.constant(beta, tf.float64))
            losses.append(float(loss))
        last = ep == tc['max_epochs'] - 1
        if n_val and (ep % tc['eval_every'] == 0 or last):
            v = predictive_nll_std(model, vb, vl, tc['inner_val_draws'])
            hist['epoch'].append(ep)
            hist['loss'].append(float(np.mean(losses)))
            hist['val_nll'].append(v)
            if v < best - 1e-6:
                best, best_epoch, wait = v, ep, 0
                best_w = [x.numpy() for x in vars_]
            else:
                wait += tc['eval_every']
                if wait >= tc['patience']:
                    break
        elif not n_val:
            hist['epoch'].append(ep)
            hist['loss'].append(float(np.mean(losses)))
        if verbose and ep % tc['verbose_every'] == 0:
            vtxt = f" val {hist['val_nll'][-1]:.3f}" if hist['val_nll'] else ''
            print(f"    epoch {ep:4d} loss {np.mean(losses):.3f}{vtxt}")
    if best_w is not None:
        for x, w in zip(vars_, best_w):
            x.assign(w)
    model['history'] = hist
    model['meta'] = {'train_rows': tr.tolist(), 'val_rows': val.tolist(), 'N': N,
                     'best_epoch': int(best_epoch), 'epochs_run': int(ep + 1), 'seed': int(seed),
                     'seconds': float(time.time() - t0)}
    return model


def train_ensemble(data, rows=None, model_cfg=None, train_cfg=None, seeds=None, verbose=False):
    """Train one model per seed (optional ensemble, §7.1).

    Inputs: as :func:`train_cf`; seeds (default train_cfg['ensemble_seeds']
    or (train_cfg['seed'],)). Outputs: list of models.
    """
    tc = train_cfg or DEFAULT_TRAIN
    seeds = seeds or tc['ensemble_seeds'] or (tc['seed'],)
    return [train_cf(data, rows, model_cfg, tc, s, verbose) for s in seeds]


# ==========================================================================
# Cross-validation folds (§8.1)
# ==========================================================================

def folds_loso(data):
    """Leave-one-station-out (= leave-one-out: one label per station).
    Inputs: data. Outputs: list of (train_rows, test_rows)."""
    n = data['n']
    return [(np.delete(np.arange(n), i), np.array([i])) for i in range(n)]


def folds_latblock(data, k=None):
    """Latitude-blocked k-fold: contiguous blocks along the transect.
    Inputs: data, k (default DEFAULT_TRAIN['cv_lat_blocks']).
    Outputs: list of (train_rows, test_rows)."""
    k = k or DEFAULT_TRAIN['cv_lat_blocks']
    order = np.argsort(data['lat'], kind='stable')
    blocks = np.array_split(order, k)
    n = data['n']
    return [(np.setdiff1d(np.arange(n), b), np.sort(b)) for b in blocks]
