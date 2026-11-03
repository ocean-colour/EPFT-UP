"""Spectral attribution maps (cf_design.md §8.7).

Methods
-------
occlusion_map          band-window deletion (primary): Δμ_k, Δσ_k, ΔNLL
integrated_gradients   per-token IG of μ_k w.r.t. r from the training-fold
                       mean spectrum, per-nm density on a common grid
moana_jacobian         analytic ∂ŷ_k/∂Rrs of retrained MOANA
moana_known_answer     IG and occlusion run on MOANA as a black box
randomization_check    maps from a model with re-initialised head weights

All Claude Fun maps are averaged over T weight draws (each draw is a fixed
function, so IG completeness holds per draw) and carry the draw spread.
Pooling attention weights are a diagnostic only (algorithm.attention_weights).
"""

import copy

import numpy as np
import tensorflow as tf
import xarray as xr

from epft_up.claude_fun import algorithm as alg
from epft_up.claude_fun import pipeline as pl
from epft_up.claude_fun import train as tr
from epft_up.claude_fun import validation as val
from epft_up.claude_fun.io import TAXA

DEFAULT_ATTR = {
    'window_nm': 10.0,
    'step_nm': 5.0,
    'ig_steps': 64,
    'n_draws': 200,                       # T weight samples per map
    'density_grid_nm': (400.0, 700.0, 1.0),
    'completeness_tol': 0.01,
    'lat_bands': (-90.0, -30.0, -10.0, 10.0, 30.0, 90.0),
    'min_per_band': 5,
    'seed': 2024,
}


def window_centres(wave_range, width=10.0, step=5.0):
    """Centres of sliding windows spanning the model range.
    Inputs: wave_range (lo, hi), width, step. Outputs: array."""
    lo, hi = wave_range
    return np.arange(lo + width / 2, hi - width / 2 + 1e-9, step)


def _drop_window(spec, lo, hi):
    out = {'wave': [], 'rrs': [], 'sigma': [], 'sst': spec['sst']}
    for w, r, s in zip(spec['wave'], spec['rrs'], spec['sigma']):
        keep = ~((w >= lo) & (w < hi))
        out['wave'].append(w[keep])
        out['rrs'].append(r[keep])
        out['sigma'].append(s[keep])
    return out


def _marginal_draws(model, spec, T, seed):
    """(mu (T,n,3), sd (T,n,3)) in log10 for each weight draw (sd = √Σ_kk)."""
    mu_t, cov_t = alg._draws(model, spec, T, seed)
    return mu_t, np.sqrt(np.einsum('tnii->tni', cov_t))


def _n_in_range(spec, cfg):
    lo, hi = cfg['wave_range_nm']
    return np.array([np.sum((w >= lo) & (w <= hi)) for w in spec['wave']])


def occlusion_map(models, spec, labels=None, cfg=None):
    """Band-window occlusion: delete all tokens in each window and record the
    change of the per-draw μ_k and σ_k, and of the predictive NLL where labels
    exist. Deletion is in-distribution for a set encoder trained with band
    dropout, so no fill value is needed.

    Inputs: models (model or list), spec, labels (optional), cfg (DEFAULT_ATTR).
    Outputs: xarray.Dataset dims (window, sample, taxon): ``d_mu``,
    ``d_mu_sd`` (spread over draws), ``d_sigma``, ``d_sigma_sd``, ``d_nll``
    (NaN without labels or for empty windows).
    """
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    models = models if isinstance(models, (list, tuple)) else [models]
    mcfg = models[0]['cfg']
    C = window_centres(mcfg['wave_range_nm'], cfg['window_nm'], cfg['step_nm'])
    n = len(spec['rrs'])
    T = cfg['n_draws']
    base = [_marginal_draws(m, spec, T, cfg['seed'] + 17 * k) for k, m in enumerate(models)]
    mu0 = np.concatenate([b[0] for b in base])
    sd0 = np.concatenate([b[1] for b in base])
    dmu = np.full((C.size, n, 3), np.nan)
    dmu_sd, dsg, dsg_sd, dnll = (np.full_like(dmu, np.nan) for _ in range(4))
    nll0 = _pred_nll(mu0, sd0, labels) if labels is not None else None
    for j, c in enumerate(C):
        so = _drop_window(spec, c - cfg['window_nm'] / 2, c + cfg['window_nm'] / 2)
        rows = np.flatnonzero(_n_in_range(so, mcfg) > 0)
        if rows.size == 0:
            continue
        sub = alg.subset(so, rows)
        dr = [_marginal_draws(m, sub, T, cfg['seed'] + 17 * k) for k, m in enumerate(models)]
        mu1 = np.concatenate([d[0] for d in dr])
        sd1 = np.concatenate([d[1] for d in dr])
        d1, d2 = mu1 - mu0[:, rows], sd1 - sd0[:, rows]
        dmu[j, rows], dmu_sd[j, rows] = d1.mean(0), d1.std(0)
        dsg[j, rows], dsg_sd[j, rows] = d2.mean(0), d2.std(0)
        if labels is not None:
            dnll[j, rows] = _pred_nll(mu1, sd1, tr.sub_labels(labels, rows)) - nll0[rows]
    return xr.Dataset({'d_mu': (('window', 'sample', 'taxon'), dmu),
                       'd_mu_sd': (('window', 'sample', 'taxon'), dmu_sd),
                       'd_sigma': (('window', 'sample', 'taxon'), dsg),
                       'd_sigma_sd': (('window', 'sample', 'taxon'), dsg_sd),
                       'd_nll': (('window', 'sample', 'taxon'), dnll)},
                      coords={'window': C, 'taxon': list(TAXA)})


def _pred_nll(mu_t, sd_t, labels):
    """Per-taxon marginal NLL of the moment-matched predictive (observed:
    Gaussian; censored: −log Φ). Outputs: (n,3)."""
    from scipy.special import log_ndtr
    m = mu_t.mean(0)
    s = np.sqrt((sd_t ** 2).mean(0) + mu_t.var(0))
    y, o, c = labels['y'], labels['observed'], labels['censored']
    z = (y - m) / s
    return np.where(o, 0.5 * np.log(2 * np.pi) + np.log(s) + 0.5 * z ** 2,
                    np.where(c, -log_ndtr(z), np.nan))


def integrated_gradients(models, spec, baseline, cfg=None):
    """Integrated gradients of μ_k (log10) w.r.t. the token feature r = Rrs/s_R.

    Path: r(α) = r0 + α(r − r0), r0 the baseline spectrum (training-fold mean)
    interpolated to the sample's own λ; σ and λ features fixed at the input's.
    Midpoint rule with ``ig_steps`` points; one fixed weight draw per pass, so
    completeness Σ attr = μ(x) − μ(baseline) holds per draw (checked).
    The encoder Jacobian along the path is computed once per sample and
    chained with each draw's head Jacobian.

    Inputs: models, spec, baseline (dict ``wave``, ``rrs``), cfg.
    Outputs: xarray.Dataset — ``density`` / ``density_sd`` (sample, grid,
    taxon) per-nm attribution on the common grid (mean / spread over draws),
    ``completeness_err`` (sample, taxon; worst relative error over draws),
    ``delta`` (sample, taxon; mean μ(x) − μ(x0)); attrs ``tokens`` is not
    stored (use :func:`ig_tokens` for raw per-token values).
    """
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    models = models if isinstance(models, (list, tuple)) else [models]
    g = np.arange(cfg['density_grid_nm'][0], cfg['density_grid_nm'][1] + 1e-9, cfg['density_grid_nm'][2])
    n = len(spec['rrs'])
    dens = np.full((n, g.size, 3), np.nan)
    dens_sd = np.full_like(dens, np.nan)
    comp = np.full((n, 3), np.nan)
    delta = np.full((n, 3), np.nan)
    for i in range(n):
        per_draw, errs, dl = [], [], []
        for k, m in enumerate(models):
            w, A, err, d = ig_tokens(m, spec, i, baseline, cfg, seed=cfg['seed'] + 17 * k)
            spacing = np.gradient(w) if w.size > 1 else np.ones(1)
            for a in A:                       # a: (3, n_tok) per draw
                per_draw.append(np.stack([np.interp(g, w, a[q] / spacing, left=np.nan, right=np.nan)
                                          for q in range(3)], axis=1))
            errs.append(err)
            dl.append(d)
        P = np.stack(per_draw)
        dens[i], dens_sd[i] = P.mean(0), P.std(0)
        comp[i] = np.max(np.concatenate(errs), axis=0)
        delta[i] = np.concatenate(dl).mean(0)
    return xr.Dataset({'density': (('sample', 'grid', 'taxon'), dens),
                       'density_sd': (('sample', 'grid', 'taxon'), dens_sd),
                       'completeness_err': (('sample', 'taxon'), comp),
                       'delta': (('sample', 'taxon'), delta)},
                      coords={'grid': g, 'taxon': list(TAXA)})


def ig_tokens(model, spec, i, baseline, cfg=None, seed=0):
    """Per-token IG for one sample and one model.

    Inputs: model, spec, i (sample index), baseline, cfg, seed.
    Outputs: (wave (m,), attr (T, 3, m) log10 per token, completeness
    relative error (T, 3), delta μ(x) − μ(x0) (T, 3)).
    """
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    mc, st = model['cfg'], model['stats']
    if mc['encoder'] != 'set':
        raise ValueError("IG is implemented for the set encoder")
    b = alg.make_batch(alg.subset(spec, [i]), st, mc)
    m = int(b['mask'][0].sum())
    w = b['wave'][0]
    x = b['x'][:, :m]
    r = x[0, :, 0].astype(np.float64)
    r0 = np.interp(w, baseline['wave'], baseline['rrs']) / st['s_R']
    S = cfg['ig_steps']
    a = (np.arange(S) + 0.5) / S
    rest = tf.constant(np.repeat(x[:, :, 1:], S, axis=0))
    Rp = tf.constant((r0[None, :] + a[:, None] * (r - r0)[None, :]).astype(np.float32))
    ones = tf.ones((S, m), tf.bool)
    with tf.GradientTape() as tape:
        tape.watch(Rp)
        hp, _ = model['encoder'](tf.concat([Rp[:, :, None], rest], -1), ones, training=False)
    try:
        Jh = tape.batch_jacobian(hp, Rp)
    except Exception:
        with tf.GradientTape() as tape:
            tape.watch(Rp)
            hp, _ = model['encoder'](tf.concat([Rp[:, :, None], rest], -1), ones, training=False)
        Jh = tape.batch_jacobian(hp, Rp, experimental_use_pfor=False)
    xe = np.concatenate([x, x], 0)
    xe[1, :, 0] = r0
    he, _ = model['encoder'](tf.constant(xe), tf.ones((2, m), tf.bool), training=False)
    sst_p = tf.constant(np.repeat(b['sst'], S, 0))
    sst_e = tf.constant(np.repeat(b['sst'], 2, 0))
    ys = tf.constant(np.asarray(st['y_sd'], np.float32))
    mode = 'sample' if mc['bayesian_head'] else 'mean'
    T = cfg['n_draws'] if mode == 'sample' else 1
    attrs, errs, deltas = [], [], []
    for t in range(T):
        with tf.GradientTape() as t2:
            t2.watch(hp)
            mu, _, _ = alg.head_forward(model, hp, sst_p, mode=mode, seed=[seed, t])
            mu = mu * ys
        G = t2.batch_jacobian(mu, hp)                                    # (S,3,d)
        grad = tf.einsum('skd,sdt->kt', G, Jh).numpy() / S             # path mean
        att = grad * (r - r0)[None, :]
        me, _, _ = alg.head_forward(model, he, sst_e, mode=mode, seed=[seed, t])
        me = (me * ys).numpy()
        d = me[0] - me[1]
        attrs.append(att)
        deltas.append(d)
        errs.append(np.abs(att.sum(1) - d) / np.maximum(np.abs(d), 1e-6))
    return w, np.stack(attrs), np.stack(errs), np.stack(deltas)


# ==========================================================================
# MOANA references
# ==========================================================================

def moana_jacobian(mm, R124, sst):
    """Analytic ∂log10ŷ_k/∂Rrs(λ) of a retrained MOANA (standardise → project
    → evaluate). With z = (R − R̄)/sd and w = Σ_j c_j V_j,
    ``J = [w − mean(w) − z (w·z)/(N−1)] / sd``; linear-space Pro adds the
    factor 1/(ŷ ln 10) (0 where ŷ is floored). Satisfies Σ_λ J·Rrs = 0 and
    Σ_λ J = 0.

    Inputs: mm (validation.moana_fit), R124 (n,124), sst (n,).
    Outputs: (n, 3, 124).
    """
    R = np.atleast_2d(np.asarray(R124, float))
    N = R.shape[1]
    Z = val.moana_standardize(R)
    sd = R.std(1, ddof=1, keepdims=True)
    out = np.zeros((R.shape[0], 3, N))
    for k, t in enumerate(TAXA):
        npc = mm['n_pc'][t]
        w = mm['coef'][t][1:1 + npc] @ mm['V'][:npc]
        J = (w[None, :] - w.mean() - Z * (Z @ w)[:, None] / (N - 1)) / sd
        if mm['space'][t] == 'linear':
            y = val.moana_raw(mm, R, sst, t)
            J = np.where((y > mm['log_clip'])[:, None], J / (np.maximum(y, 1e-30)[:, None] * np.log(10)), 0.0)
        out[:, k] = J
    return out


def ig_blackbox(f, x, x0, steps=64, grad=None, h_rel=1e-5):
    """Integrated gradients of a scalar black-box f (midpoint rule); gradients
    by central differences unless ``grad`` is given.
    Inputs: f (callable (b,m)→(b,)), x, x0 (m,), steps, grad, h_rel.
    Outputs: attributions (m,)."""
    x, x0 = np.asarray(x, float), np.asarray(x0, float)
    a = (np.arange(steps) + 0.5) / steps
    P = x0[None] + a[:, None] * (x - x0)[None]
    if grad is None:
        h = h_rel * max(np.max(np.abs(x)), 1e-12)
        m = x.size
        E = np.eye(m) * h
        G = np.zeros_like(P)
        for s, p in enumerate(P):
            G[s] = (f(p[None] + E) - f(p[None] - E)) / (2 * h)
    else:
        G = np.stack([grad(p) for p in P])
    return G.mean(0) * (x - x0)


def occlusion_blackbox(f, wave, x, centres, width, grid):
    """Occlusion for a fixed-grid black box: delete bands of ``wave`` in each
    window, re-interpolate to ``grid`` (the model interpolates across the gap;
    edges are held constant) and record f(occluded) − f(full).
    Inputs: f ((b,m)→(b,)), wave, x (native), centres, width, grid.
    Outputs: (n_windows,) deltas."""
    full = np.interp(grid, wave, x)
    base = f(full[None])[0]
    out = np.zeros(len(centres))
    for j, c in enumerate(centres):
        keep = ~((wave >= c - width / 2) & (wave < c + width / 2))
        out[j] = f(np.interp(grid, wave[keep], x[keep])[None])[0] - base
    return out


def moana_known_answer(mm, wave, R_native, sst, baseline124, cfg=None, pipeline_cfg=None):
    """Known-answer check (§8.7): run IG and occlusion on retrained MOANA as a
    black box. IG must recover the analytic Jacobian (correlation of IG with
    J(x)⊙(x − x0) > 0.9 over 414–660 nm); occlusion must put its largest
    effects in the same windows as the linearised Jacobian prediction J·ΔR.

    Inputs: mm (moana_fit), wave (native), R_native (n, m), sst (n,),
    baseline124 (124,), cfg, pipeline_cfg.
    Outputs: dict of (n,3) arrays ``ig_corr``, ``ig_completeness``,
    ``occ_corr`` (Δ vs J·ΔR over windows), ``occ_top3_overlap`` (0–3), plus
    ``pass_ig`` (bool, median ig_corr > 0.9 for every taxon) and the sum-rule
    residuals ``sum_JR`` and ``sum_J`` (max |·| relative).
    """
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    g = pl.moana_grid(pipeline_cfg)
    C = window_centres((g[0], g[-1]), cfg['window_nm'], cfg['step_nm'])
    n = len(R_native)
    out = {k: np.full((n, 3), np.nan) for k in ('ig_corr', 'ig_completeness', 'occ_corr', 'occ_top3_overlap')}
    R124 = np.vstack([np.interp(g, wave, r) for r in R_native])
    J = moana_jacobian(mm, R124, sst)
    scale = np.abs(J).max(axis=2) * np.abs(R124).sum(axis=1)[:, None]
    nz = scale > 0                      # floored linear Pro has J = 0
    out['sum_JR'] = float(np.max(np.abs(np.einsum('nkm,nm->nk', J, R124))[nz] / scale[nz])) if nz.any() else 0.0
    out['sum_J'] = float(np.max(np.abs(J.sum(2))[nz] / np.abs(J).sum(2)[nz])) if nz.any() else 0.0
    for i in range(n):
        for k, t in enumerate(TAXA):
            f = lambda X: val.moana_log10(mm, X, np.full(len(X), sst[i]), t)
            ig = ig_blackbox(f, R124[i], baseline124, cfg['ig_steps'])
            target = J[i, k] * (R124[i] - baseline124)
            if np.std(ig) == 0 or np.std(target) == 0:
                continue                # e.g. linear Pro floored at log_clip: no sensitivity
            out['ig_corr'][i, k] = np.corrcoef(ig, target)[0, 1]
            d = f(R124[i][None])[0] - f(baseline124[None])[0]
            out['ig_completeness'][i, k] = abs(ig.sum() - d) / max(abs(d), 1e-9)
            occ = occlusion_blackbox(f, wave, R_native[i], C, cfg['window_nm'], g)
            lin = np.array([J[i, k] @ (np.interp(g, wave[~((wave >= c - cfg['window_nm'] / 2) & (wave < c + cfg['window_nm'] / 2))],
                                                 R_native[i][~((wave >= c - cfg['window_nm'] / 2) & (wave < c + cfg['window_nm'] / 2))])
                                       - R124[i]) for c in C])
            out['occ_corr'][i, k] = np.corrcoef(occ, lin)[0, 1] if np.std(occ) > 0 else np.nan
            out['occ_top3_overlap'][i, k] = len(set(np.argsort(-np.abs(occ))[:3]) & set(np.argsort(-np.abs(lin))[:3]))
    out['pass_ig'] = bool(np.all(np.nanmedian(out['ig_corr'], axis=0) > 0.9))
    out['windows'] = C
    return out


# ==========================================================================
# Sanity checks and orchestration
# ==========================================================================

def randomized_head(model):
    """Copy of a model with the trained encoder and a freshly initialised head
    (model-randomisation check). Inputs: model. Outputs: model dict."""
    new = alg.build_model(model['cfg'], model['stats'])
    for a, b in zip(model['encoder'].weights, new['encoder'].weights):
        b.assign(a)
    new['meta'] = {'randomized_head': True}
    return new


def randomization_check(models, spec, cfg=None):
    """Correlation (over windows) between sample-mean occlusion Δμ maps of the
    trained and head-randomised models; structure should not survive.
    Inputs: models, spec, cfg. Outputs: dict taxon → correlation."""
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    models = models if isinstance(models, (list, tuple)) else [models]
    a = occlusion_map(models, spec, None, cfg)['d_mu'].mean('sample')
    b = occlusion_map([randomized_head(m) for m in models], spec, None, cfg)['d_mu'].mean('sample')
    out = {}
    for t in TAXA:
        x, y = a.sel(taxon=t).values, b.sel(taxon=t).values
        ok = np.isfinite(x) & np.isfinite(y)
        out[t] = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 and np.std(y[ok]) > 0 else np.nan
    return out


def maps_for(models, spec, labels, baseline, mm=None, R124=None, sst=None, lat=None, cfg=None,
             randomize=False):
    """All maps for one set of held-out samples.

    Inputs: models; spec, labels (held-out); baseline (training mean spectrum
    dict); mm (retrained MOANA, optional) with R124, sst for its Jacobian;
    lat (for stratification); cfg; randomize (also run the randomisation
    check).
    Outputs: dict ``occlusion``, ``ig`` (Datasets), ``moana_jacobian``
    (n,3,124) or None, ``attention`` (list), ``lat``, ``randomization``.
    """
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    out = {'occlusion': occlusion_map(models, spec, labels, cfg),
           'ig': integrated_gradients(models, spec, baseline, cfg),
           'moana_jacobian': moana_jacobian(mm, R124, sst) if mm is not None else None,
           'lat': np.asarray(lat) if lat is not None else None,
           'randomization': randomization_check(models, spec, cfg) if randomize else None}
    m0 = models[0] if isinstance(models, (list, tuple)) else models
    out['attention'] = alg.attention_weights(m0, spec) if m0['cfg']['encoder'] == 'set' else None
    return out


def fold_maps(models, data, trn, tst, cfg=None, train_cfg=None, pipeline_cfg=None, randomize=False):
    """``run_cv`` callback: maps for one CV fold (baseline and MOANA from the
    training rows). Inputs: models, data, trn, tst, cfgs. Outputs: maps_for
    dict plus ``test_rows``."""
    tc = train_cfg or tr.DEFAULT_TRAIN
    spec = alg.subset(data['spec'], tst)
    labels = tr.sub_labels(data['labels'], tst)
    baseline = {'wave': data['wave'], 'rrs': np.nanmean(data['rrs'][trn], axis=0)}
    R124 = val.to_moana_grid(data, cfg=pipeline_cfg)
    sst = np.asarray(data['spec']['sst'])
    mm = val.moana_fit(R124[trn], tr.sub_labels(data['labels'], trn), data['ds']['abundance'].values[trn],
                       sst[trn], tc['moana'])
    out = maps_for(models, spec, labels, baseline, mm, R124[tst], sst[tst], data['lat'][tst], cfg, randomize)
    out['test_rows'] = np.asarray(tst)
    return out


def combine_maps(maps, sample_index=None):
    """Stack per-fold (or per-set) map dicts into one Dataset over samples,
    with fold means and spreads.

    Inputs: maps (list of maps_for/fold_maps dicts).
    Outputs: xarray.Dataset — per-sample ``occ_d_mu``, ``occ_d_sigma``,
    ``occ_d_nll`` (sample, window, taxon), ``ig_density`` (sample, grid,
    taxon), ``ig_completeness_err``, ``moana_jacobian`` (sample, taxon,
    moana_wave), ``lat``, ``fold``; fold summaries ``occ_d_mu_fold``
    (fold, window, taxon) and ``ig_density_fold`` (fold, grid, taxon).
    """
    occ = xr.concat([m['occlusion'] for m in maps], dim='sample')
    ig = xr.concat([m['ig'] for m in maps], dim='sample')
    fold = np.concatenate([np.full(m['occlusion'].sizes['sample'], f) for f, m in enumerate(maps)])
    ds = xr.Dataset({'occ_d_mu': occ['d_mu'].transpose('sample', 'window', 'taxon'),
                     'occ_d_mu_sd': occ['d_mu_sd'].transpose('sample', 'window', 'taxon'),
                     'occ_d_sigma': occ['d_sigma'].transpose('sample', 'window', 'taxon'),
                     'occ_d_nll': occ['d_nll'].transpose('sample', 'window', 'taxon'),
                     'ig_density': ig['density'], 'ig_density_sd': ig['density_sd'],
                     'ig_completeness_err': ig['completeness_err'], 'fold_index': ('sample', fold)})
    if all(m['moana_jacobian'] is not None for m in maps):
        J = np.concatenate([m['moana_jacobian'] for m in maps])
        ds['moana_jacobian'] = (('sample', 'taxon', 'moana_wave'), J)
        ds = ds.assign_coords(moana_wave=pl.moana_grid())
    if all(m['lat'] is not None for m in maps):
        ds['lat'] = ('sample', np.concatenate([m['lat'] for m in maps]))
    if 'test_rows' in maps[0]:
        ds['row'] = ('sample', np.concatenate([m['test_rows'] for m in maps]))
    for v in ('occ_d_mu', 'ig_density'):
        g = ds[v].assign_coords(fold=('sample', fold)).groupby('fold').mean('sample')
        ds[v + '_fold'] = g.rename(fold='fold_id').drop_vars('fold_id', errors='ignore')
    rnd = [m['randomization'] for m in maps if m.get('randomization')]
    if rnd:
        for t in TAXA:
            ds.attrs[f'randomization_corr_{t}'] = float(np.nanmean([r[t] for r in rnd]))
    return ds


def stratify_by_latitude(ds, cfg=None):
    """Mean maps per latitude band where n ≥ min_per_band.
    Inputs: combine_maps Dataset, cfg. Outputs: Dataset (band, window/grid, taxon)."""
    cfg = {**DEFAULT_ATTR, **(cfg or {})}
    edges = np.asarray(cfg['lat_bands'])
    band = np.digitize(ds['lat'].values, edges) - 1
    keep = [b for b in np.unique(band) if np.sum(band == b) >= cfg['min_per_band']]
    parts = []
    for b in keep:
        sel = ds.isel(sample=np.flatnonzero(band == b))
        parts.append(xr.Dataset({'occ_d_mu': sel['occ_d_mu'].mean('sample'),
                                 'ig_density': sel['ig_density'].mean('sample'),
                                 'n': sel.sizes['sample']}).assign_coords(
            band=f"{edges[b]:.0f}..{edges[b + 1]:.0f}"))
    return xr.concat(parts, dim='band') if parts else None
