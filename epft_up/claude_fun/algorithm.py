"""Claude Fun model: tokenisation, standardisation, set encoder, Bayesian head
(cf_design.md §5–6).

Input: a set of per-band tokens (λ, Rrs, σ_Rrs) of any length and grid, plus
SST. Output: a joint Gaussian over log10(cells mL⁻¹) of Pro, Syn and Peuk with
a 3×3 covariance Σ = D R D.

State is plain dicts: a *model* is ``{'cfg', 'stats', 'encoder', 'head'}``
where encoder and head are custom Keras layers (the one permitted class
exception, §6.5/§9). A *spectra set* is ``{'wave', 'rrs', 'sigma'}`` (lists of
1-D arrays, one per sample, any grid) plus ``'sst'`` (n,).
"""

import copy
import json
import os

import numpy as np
import tensorflow as tf
import keras
import xarray as xr
from scipy.special import log_ndtr

from epft_up.claude_fun.io import TAXA

DEFAULT_MODEL = {
    # ---- input (§5) -----------------------------------------------------------
    'wave_range_nm': (400.0, 700.0),      # native 400.5–698.0 nm = 86 bands
    'n_fourier': 16,                       # 16 periods → 32 λ features
    'fourier_period_nm': (5.0, 500.0),     # log-spaced
    'lambda0_nm': 300.0,
    'use_sst': True,
    # ---- encoder (§6.1) ---------------------------------------------------------
    'encoder': 'set',                      # 'set' | 'counterbin' (§8.3 ablation)
    'd_model': 64, 'n_heads': 4, 'd_ff': 128, 'n_blocks': 2, 'dropout': 0.1,
    'mc_dropout': False,                   # MC dropout in the encoder at prediction
    'counterbin_step_nm': 1.0,             # fine grid of the counter-binning ablation
    # ---- head (§6.2) ------------------------------------------------------------
    'head_units': (64, 64),
    'bayesian_head': True,                 # False → deterministic head, no KL (ablation)
    'prior_sigma': 1.0,                    # p(w) = N(0, σ_w²)
    'posterior_init_sigma': 0.01,          # initial q(w) standard deviation
    'sigma_min': 0.01,                     # standardised units
    'corr_mode': 'shared',                 # 'shared' | 'input' | 'diag'
    'rho_max': 0.999,                      # |ρ| ≤ rho_max keeps R numerically PD
    # ---- likelihood (§5.4, §6.3) ----------------------------------------------
    'counting_noise': True,                # Σ_eff = Σ + diag(σ_count²)
    # ---- prediction (§6.4) ----------------------------------------------------
    'n_samples': 200,
}

LOG10E = float(np.log10(np.e))
_LOG2PI = float(np.log(2 * np.pi))


def config(**overrides):
    """Deep copy of :data:`DEFAULT_MODEL` with overrides. Outputs: dict."""
    cfg = copy.deepcopy(DEFAULT_MODEL)
    unknown = set(overrides) - set(cfg)
    if unknown:
        raise KeyError(f"unknown model keys {sorted(unknown)}")
    cfg.update(overrides)
    return cfg


# ==========================================================================
# Spectra sets, standardisation, tokens
# ==========================================================================

def spectra_from_dataset(ds, rows=None):
    """Spectra set from a training-matrix / transfer dataset (median spectra).

    Inputs: ds (pipeline dataset with rrs, sigma, sst, wavelength); rows.
    Outputs: dict ``wave``, ``rrs``, ``sigma`` (lists of 1-D arrays, non-finite
    bands dropped per sample) and ``sst`` (n,).
    """
    rows = np.arange(ds.sizes['sample']) if rows is None else np.asarray(rows)
    w = ds['wavelength'].values
    R, S = ds['rrs'].values[rows], ds['sigma'].values[rows]
    out = {'wave': [], 'rrs': [], 'sigma': [], 'sst': ds['sst'].values[rows].astype(float)}
    for r, s in zip(R, S):
        ok = np.isfinite(r) & np.isfinite(s)
        out['wave'].append(w[ok])
        out['rrs'].append(r[ok])
        out['sigma'].append(s[ok])
    return out


def labels_from_dataset(ds, rows=None):
    """Label arrays (log10 space) from a pipeline dataset.

    Inputs: ds, rows. Outputs: dict ``y`` (n,3; log10 DL where censored, NaN
    where missing), ``observed``, ``censored`` (n,3 bool), ``sigma_count``
    (n,3; 0 unless observed).
    """
    rows = np.arange(ds.sizes['sample']) if rows is None else np.asarray(rows)
    return {'y': ds['y_log10'].values[rows].astype(float),
            'observed': ds['observed'].values[rows].astype(bool),
            'censored': ds['censored'].values[rows].astype(bool),
            'sigma_count': np.nan_to_num(ds['sigma_count'].values[rows].astype(float))}


def subset(spec, rows):
    """Rows of a spectra set. Inputs: spec, rows. Outputs: spectra set."""
    rows = np.asarray(rows)
    return {'wave': [spec['wave'][i] for i in rows], 'rrs': [spec['rrs'][i] for i in rows],
            'sigma': [spec['sigma'][i] for i in rows], 'sst': np.asarray(spec['sst'])[rows]}


def fit_stats(spec, labels, cfg):
    """Training-fold standardisation statistics (§5.2–5.4).

    ``s_R`` is one global scalar, the median Rrs over all in-range bands and
    samples (no per-spectrum normalisation: magnitude carries biomass).

    Inputs: spec, labels (training rows only), cfg.
    Outputs: dict ``s_R``, ``sst_mean``, ``sst_sd``, ``y_mean``, ``y_sd`` (3,).
    """
    lo, hi = cfg['wave_range_nm']
    vals = np.concatenate([r[(w >= lo) & (w <= hi)] for w, r in zip(spec['wave'], spec['rrs'])])
    s_R = float(np.median(vals[np.isfinite(vals)]))
    sst = np.asarray(spec['sst'], float)
    y, o = labels['y'], labels['observed']
    ym = np.array([np.mean(y[o[:, k], k]) for k in range(3)])
    ys = np.array([np.std(y[o[:, k], k], ddof=1) if o[:, k].sum() > 1 else 1.0 for k in range(3)])
    return {'s_R': s_R, 'sst_mean': float(np.nanmean(sst)),
            'sst_sd': float(np.nanstd(sst, ddof=1)) if np.isfinite(sst).sum() > 1 else 1.0,
            'y_mean': ym.tolist(), 'y_sd': np.where(ys > 0, ys, 1.0).tolist()}


def fourier_features(wave, cfg):
    """sin/cos(2π(λ−λ0)/P_k) for log-spaced P_k (§5.2). Continuous in λ.

    Inputs: wave (n,), cfg. Outputs: (n, 2·n_fourier) float32.
    """
    p0, p1 = cfg['fourier_period_nm']
    P = p0 * (p1 / p0) ** (np.arange(cfg['n_fourier']) / max(cfg['n_fourier'] - 1, 1))
    a = 2 * np.pi * (np.asarray(wave, float)[:, None] - cfg['lambda0_nm']) / P[None, :]
    return np.concatenate([np.sin(a), np.cos(a)], axis=1).astype(np.float32)


def tokens(wave, rrs, sigma, stats, cfg):
    """Per-band tokens [r, u, λ-features] inside the model's wavelength range.

    r = Rrs/s_R, u = σ_Rrs/s_R (§5.2).
    Inputs: wave, rrs, sigma (1-D), stats, cfg.
    Outputs: (tokens (m, 2+2K) float32, wave_used (m,)).
    """
    lo, hi = cfg['wave_range_nm']
    w = np.asarray(wave, float)
    r = np.asarray(rrs, float)
    s = np.asarray(sigma, float)
    ok = (w >= lo) & (w <= hi) & np.isfinite(r) & np.isfinite(s)
    w, r, s = w[ok], r[ok], s[ok]
    t = np.concatenate([(r / stats['s_R'])[:, None], (s / stats['s_R'])[:, None],
                        fourier_features(w, cfg)], axis=1).astype(np.float32)
    return t, w


def counterbin(wave, rrs, sigma, grid):
    """Counter-binning onto a fixed fine grid (ablation input, §8.3).

    Each input band owns the interval between the midpoints to its neighbours
    (edge bands extend to the grid edges); every fine bin in that interval
    takes the band's value with σ·sqrt(k), k = owned width / fine step, so
    the band's information (Σ 1/σ²) is conserved, not duplicated.

    Inputs: wave, rrs, sigma (1-D, sorted), grid (fine, regular).
    Outputs: (rrs_fine, sigma_fine) on ``grid``.
    """
    w = np.asarray(wave, float)
    edges = np.concatenate([[-np.inf], (w[1:] + w[:-1]) / 2, [np.inf]])
    step = grid[1] - grid[0]
    lo_g, hi_g = grid[0] - step / 2, grid[-1] + step / 2
    width = np.minimum(edges[1:], hi_g) - np.maximum(edges[:-1], lo_g)
    k = np.maximum(width / step, 1.0)
    idx = np.clip(np.searchsorted(edges, grid, side='right') - 1, 0, w.size - 1)
    return np.asarray(rrs)[idx], np.asarray(sigma)[idx] * np.sqrt(k[idx])


def make_batch(spec, stats, cfg, rows=None, pad_to=None):
    """Padded model inputs for a spectra set.

    Inputs: spec, stats, cfg, rows (subset), pad_to (int; default the max
    token count rounded up to a multiple of 32, which limits retracing).
    Outputs: dict ``x`` (B,T,F) float32 tokens or (B, 2G) counter-binned
    features, ``mask`` (B,T) bool, ``sst`` (B,1) float32 (standardised),
    ``wave`` list of token wavelengths per sample.
    """
    rows = range(len(spec['rrs'])) if rows is None else rows
    sst = np.asarray(spec['sst'], float)
    z = ((sst - stats['sst_mean']) / stats['sst_sd']).astype(np.float32)
    if cfg['encoder'] == 'counterbin':
        lo, hi = cfg['wave_range_nm']
        g = np.arange(lo, hi + 1e-9, cfg['counterbin_step_nm'])
        X, W = [], []
        for i in rows:
            t, w = tokens(spec['wave'][i], spec['rrs'][i], spec['sigma'][i], stats, cfg)
            r, s = counterbin(w, t[:, 0], t[:, 1], g)
            X.append(np.concatenate([r, s]).astype(np.float32))
            W.append(w)
        return {'x': np.stack(X), 'mask': None, 'sst': z[list(rows)][:, None], 'wave': W}
    toks, W = [], []
    for i in rows:
        t, w = tokens(spec['wave'][i], spec['rrs'][i], spec['sigma'][i], stats, cfg)
        if t.shape[0] == 0:
            raise ValueError(f"sample {i} has no bands in {cfg['wave_range_nm']}")
        toks.append(t)
        W.append(w)
    T = int(32 * np.ceil(max(t.shape[0] for t in toks) / 32))
    T = max(T, pad_to or 0)
    F = toks[0].shape[1]
    x = np.zeros((len(toks), T, F), np.float32)
    m = np.zeros((len(toks), T), bool)
    for i, t in enumerate(toks):
        x[i, :t.shape[0]] = t
        m[i, :t.shape[0]] = True
    return {'x': x, 'mask': m, 'sst': z[list(rows)][:, None], 'wave': W}


def standardise_labels(labels, stats, counting_noise=True):
    """Labels in standardised log space for the NLL.

    Inputs: labels (labels_from_dataset), stats, counting_noise (bool).
    Outputs: dict ``y`` (n,3; censored entries hold the standardised log10 DL
    threshold; missing → 0), ``obs``, ``cens`` (float masks), ``cvar``
    (σ_count² in standardised units; 0 unless observed and enabled).
    """
    m, s = np.asarray(stats['y_mean']), np.asarray(stats['y_sd'])
    y = (np.nan_to_num(labels['y']) - m) / s
    obs = labels['observed'].astype(np.float32)
    cens = labels['censored'].astype(np.float32)
    cv = (labels['sigma_count'] / s) ** 2 * obs if counting_noise else np.zeros_like(obs)
    return {'y': np.where((obs + cens) > 0, y, 0.0).astype(np.float32), 'obs': obs,
            'cens': cens, 'cvar': cv.astype(np.float32)}


# ==========================================================================
# Custom Keras layers (the permitted class exception, §6.5)
# ==========================================================================

def _softplus_inv(x):
    return float(np.log(np.expm1(x)))


@keras.saving.register_keras_serializable('claude_fun')
class VariationalDense(keras.layers.Layer):
    """Dense layer with a mean-field Gaussian posterior over weights and biases.

    Modes: ``'lrt'`` (local reparameterisation: per-example pre-activation
    noise, used in training), ``'sample'`` (one weight draw shared by the batch,
    from a stateless seed — used for prediction and attribution so a draw is a
    fixed function), ``'mean'`` (posterior means). With ``bayesian=False``
    the layer is an ordinary Dense layer and its KL is zero.
    """

    def __init__(self, units, activation=None, prior_sigma=1.0, init_sigma=0.01,
                 bayesian=True, **kw):
        super().__init__(**kw)
        self.units = int(units)
        self.activation = keras.activations.get(activation)
        self.prior_sigma = float(prior_sigma)
        self.init_sigma = float(init_sigma)
        self.bayesian = bool(bayesian)

    def build(self, input_shape):
        d = int(input_shape[-1])
        self.w_mu = self.add_weight(name='w_mu', shape=(d, self.units), initializer='glorot_uniform')
        self.b_mu = self.add_weight(name='b_mu', shape=(self.units,), initializer='zeros')
        if self.bayesian:
            c = keras.initializers.Constant(_softplus_inv(self.init_sigma))
            self.w_rho = self.add_weight(name='w_rho', shape=(d, self.units), initializer=c)
            self.b_rho = self.add_weight(name='b_rho', shape=(self.units,), initializer=c)

    def call(self, x, mode='lrt', seed=None):
        if not self.bayesian or mode == 'mean':
            return self.activation(tf.matmul(x, self.w_mu) + self.b_mu)
        ws, bs = tf.nn.softplus(self.w_rho), tf.nn.softplus(self.b_rho)
        if mode == 'lrt':
            m = tf.matmul(x, self.w_mu) + self.b_mu
            v = tf.matmul(tf.square(x), tf.square(ws)) + tf.square(bs)
            return self.activation(m + tf.sqrt(v + 1e-12) * tf.random.normal(tf.shape(m)))
        if mode == 'sample':
            s1 = tf.cast(seed, tf.int32)
            ew = tf.random.stateless_normal(tf.shape(self.w_mu), seed=s1)
            eb = tf.random.stateless_normal(tf.shape(self.b_mu), seed=s1 + [0, 7919])
            return self.activation(tf.matmul(x, self.w_mu + ws * ew) + self.b_mu + bs * eb)
        raise ValueError(mode)

    def kl(self):
        """KL(q‖p) with p = N(0, σ_w²), summed over weights and biases."""
        if not self.bayesian:
            return tf.constant(0.0)
        out = 0.0
        p2 = self.prior_sigma ** 2
        for mu, rho in ((self.w_mu, self.w_rho), (self.b_mu, self.b_rho)):
            s = tf.nn.softplus(rho)
            out += tf.reduce_sum(tf.math.log(self.prior_sigma / s)
                                 + (tf.square(s) + tf.square(mu)) / (2 * p2) - 0.5)
        return out

    def get_config(self):
        c = super().get_config()
        c.update(units=self.units, activation=keras.activations.serialize(self.activation),
                 prior_sigma=self.prior_sigma, init_sigma=self.init_sigma, bayesian=self.bayesian)
        return c


@keras.saving.register_keras_serializable('claude_fun')
class SharedCorr(keras.layers.Layer):
    """Three variational parameters z (no input dependence) for the global
    correlation structure of ``corr_mode='shared'``; same prior as the head."""

    def __init__(self, prior_sigma=1.0, init_sigma=0.01, bayesian=True, **kw):
        super().__init__(**kw)
        self.prior_sigma, self.init_sigma, self.bayesian = float(prior_sigma), float(init_sigma), bool(bayesian)

    def build(self, input_shape=None):
        self.z_mu = self.add_weight(name='z_mu', shape=(3,), initializer='zeros')
        if self.bayesian:
            self.z_rho = self.add_weight(name='z_rho', shape=(3,),
                                         initializer=keras.initializers.Constant(_softplus_inv(self.init_sigma)))

    def call(self, mode='lrt', seed=None):
        if not self.bayesian or mode == 'mean':
            return self.z_mu
        s = tf.nn.softplus(self.z_rho)
        if mode == 'sample':
            e = tf.random.stateless_normal((3,), seed=tf.cast(seed, tf.int32) + [0, 104729])
        else:
            e = tf.random.normal((3,))
        return self.z_mu + s * e

    def kl(self):
        if not self.bayesian:
            return tf.constant(0.0)
        s = tf.nn.softplus(self.z_rho)
        p2 = self.prior_sigma ** 2
        return tf.reduce_sum(tf.math.log(self.prior_sigma / s)
                             + (tf.square(s) + tf.square(self.z_mu)) / (2 * p2) - 0.5)


@keras.saving.register_keras_serializable('claude_fun')
class SetEncoder(keras.layers.Layer):
    """Token MLP → pre-LN transformer blocks with masked attention → attention
    pooling with one learned seed query (§6.1). Permutation-invariant: there is
    no positional encoding beyond the λ features."""

    def __init__(self, d_model=64, n_heads=4, d_ff=128, n_blocks=2, dropout=0.1, **kw):
        super().__init__(**kw)
        self.d, self.h, self.ff, self.nb, self.p = d_model, n_heads, d_ff, n_blocks, dropout
        self.tok1 = keras.layers.Dense(d_model, activation='gelu')
        self.tok2 = keras.layers.Dense(d_model)
        self.blocks = []
        for _ in range(n_blocks):
            self.blocks.append({
                'ln1': keras.layers.LayerNormalization(epsilon=1e-5),
                'att': keras.layers.MultiHeadAttention(n_heads, d_model // n_heads, dropout=dropout),
                'dr1': keras.layers.Dropout(dropout),
                'ln2': keras.layers.LayerNormalization(epsilon=1e-5),
                'ff1': keras.layers.Dense(d_ff, activation='gelu'),
                'ff2': keras.layers.Dense(d_model),
                'dr2': keras.layers.Dropout(dropout)})
        self.ln_out = keras.layers.LayerNormalization(epsilon=1e-5)
        self.pool = keras.layers.MultiHeadAttention(n_heads, d_model // n_heads)

    def build(self, input_shape):
        self.seed_query = self.add_weight(name='seed_query', shape=(1, 1, self.d),
                                          initializer='glorot_uniform')

    def call(self, x, mask, training=False):
        mask = tf.cast(mask, tf.bool)
        h = self.tok2(self.tok1(x))
        am = tf.logical_and(mask[:, :, None], mask[:, None, :])
        # padded query rows get a self-edge only, so their softmax is defined
        am = tf.logical_or(am, tf.eye(tf.shape(x)[1], batch_shape=[tf.shape(x)[0]], dtype=tf.bool))
        for b in self.blocks:
            a = b['att'](b['ln1'](h), b['ln1'](h), attention_mask=am, training=training)
            h = h + b['dr1'](a, training=training)
            h = h + b['dr2'](b['ff2'](b['ff1'](b['ln2'](h))), training=training)
        h = self.ln_out(h)
        q = tf.tile(self.seed_query, [tf.shape(x)[0], 1, 1])
        out, sc = self.pool(q, h, attention_mask=mask[:, None, :], return_attention_scores=True,
                            training=training)
        return out[:, 0, :], tf.reduce_mean(sc[:, :, 0, :], axis=1)


@keras.saving.register_keras_serializable('claude_fun')
class GridEncoder(keras.layers.Layer):
    """Counter-binning ablation encoder: fixed fine grid of (r, u) → Dense(d)."""

    def __init__(self, d_model=64, **kw):
        super().__init__(**kw)
        self.d = d_model
        self.dense = keras.layers.Dense(d_model, activation='gelu')

    def call(self, x, mask=None, training=False):
        h = self.dense(x)
        return h, tf.zeros((tf.shape(x)[0], 1))


@keras.saving.register_keras_serializable('claude_fun')
class BayesHead(keras.layers.Layer):
    """Two variational dense layers (GELU) and a variational output layer
    giving μ_k, s_k (σ_k = softplus(s_k) + σ_min) and, for corr_mode='input',
    z per sample; 'shared' uses :class:`SharedCorr` (§6.2)."""

    def __init__(self, units=(64, 64), corr_mode='shared', prior_sigma=1.0, init_sigma=0.01,
                 bayesian=True, **kw):
        super().__init__(**kw)
        self.corr_mode = corr_mode
        kw2 = dict(prior_sigma=prior_sigma, init_sigma=init_sigma, bayesian=bayesian)
        self.hidden = [VariationalDense(u, activation='gelu', **kw2) for u in units]
        n_out = 9 if corr_mode == 'input' else 6
        self.out = VariationalDense(n_out, **kw2)
        self.corr = SharedCorr(**kw2) if corr_mode == 'shared' else None

    def call(self, x, mode='lrt', seed=None):
        for i, layer in enumerate(self.hidden):
            x = layer(x, mode=mode, seed=None if seed is None else seed + [0, 31 * (i + 1)])
        o = self.out(x, mode=mode, seed=None if seed is None else seed + [0, 997])
        if self.corr_mode == 'input':
            z = o[:, 6:9]
        elif self.corr_mode == 'shared':
            z = tf.tile(self.corr(mode=mode, seed=seed)[None, :], [tf.shape(o)[0], 1])
        else:
            z = None
        return o[:, :3], o[:, 3:6], z

    def kl(self):
        k = tf.add_n([l.kl() for l in self.hidden + [self.out]])
        return k + (self.corr.kl() if self.corr is not None else 0.0)

    def build(self, input_shape):
        if self.corr is not None:
            self.corr.build()


# ==========================================================================
# Model assembly and forward pass
# ==========================================================================

def build_model(cfg=None, stats=None, n_features=None):
    """Create a model dict and build its weights with a dummy pass.

    Inputs: cfg (DEFAULT_MODEL-like), stats (fit_stats), n_features (token or
    flattened grid width; inferred from cfg if None).
    Outputs: dict ``cfg``, ``stats``, ``encoder``, ``head``.
    """
    cfg = copy.deepcopy(cfg or DEFAULT_MODEL)
    if cfg['encoder'] == 'set':
        enc = SetEncoder(cfg['d_model'], cfg['n_heads'], cfg['d_ff'], cfg['n_blocks'], cfg['dropout'])
        F = n_features or 2 + 2 * cfg['n_fourier']
        enc(tf.zeros((1, 4, F)), tf.ones((1, 4), tf.bool))
    elif cfg['encoder'] == 'counterbin':
        lo, hi = cfg['wave_range_nm']
        G = np.arange(lo, hi + 1e-9, cfg['counterbin_step_nm']).size
        enc = GridEncoder(cfg['d_model'])
        enc(tf.zeros((1, n_features or 2 * G)))
    else:
        raise ValueError(f"encoder must be 'set' or 'counterbin', got {cfg['encoder']!r}")
    if cfg['corr_mode'] not in ('shared', 'input', 'diag'):
        raise ValueError(f"corr_mode {cfg['corr_mode']!r}")
    head = BayesHead(cfg['head_units'], cfg['corr_mode'], cfg['prior_sigma'],
                     cfg['posterior_init_sigma'], cfg['bayesian_head'])
    head(tf.zeros((1, cfg['d_model'] + int(cfg['use_sst']))), mode='mean')
    return {'cfg': cfg, 'stats': stats, 'encoder': enc, 'head': head}


def trainable_variables(model):
    """All trainable variables of the encoder and head."""
    return model['encoder'].trainable_variables + model['head'].trainable_variables


def kl(model):
    """KL(q(w)‖p(w)) of the head (the encoder is deterministic, §6.2)."""
    return model['head'].kl()


def corr_cholesky(z, rho_max=0.999):
    """Cholesky factor of a 3×3 correlation matrix from canonical partial
    correlations ρ12 = tanh(z1), ρ13 = tanh(z2), ρ23|1 = tanh(z3) (§6.2);
    valid (positive-definite) for any z. Works for tf tensors and numpy.

    Inputs: z (B,3) or None (→ identity). Outputs: L (B,3,3) with R = L Lᵀ.
    """
    if isinstance(z, np.ndarray):
        t = rho_max * np.tanh(z)
        r12, r13, r23 = t[:, 0], t[:, 1], t[:, 2]
        c12, c13, c23 = np.sqrt(1 - r12 ** 2), np.sqrt(1 - r13 ** 2), np.sqrt(1 - r23 ** 2)
        L = np.zeros((z.shape[0], 3, 3))
        L[:, 0, 0] = 1
        L[:, 1, 0], L[:, 1, 1] = r12, c12
        L[:, 2, 0], L[:, 2, 1], L[:, 2, 2] = r13, r23 * c13, c13 * c23
        return L
    t = rho_max * tf.tanh(z)
    r12, r13, r23 = t[:, 0], t[:, 1], t[:, 2]
    c12, c13, c23 = tf.sqrt(1 - r12 ** 2), tf.sqrt(1 - r13 ** 2), tf.sqrt(1 - r23 ** 2)
    o, zz = tf.ones_like(r12), tf.zeros_like(r12)
    return tf.stack([tf.stack([o, zz, zz], -1), tf.stack([r12, c12, zz], -1),
                     tf.stack([r13, r23 * c13, c13 * c23], -1)], axis=1)


def encode(model, batch, training=False):
    """Encoder pass. Inputs: model, batch (make_batch). Outputs: (h, attn)."""
    if model['cfg']['encoder'] == 'counterbin':
        return model['encoder'](tf.convert_to_tensor(batch['x']), training=training)
    return model['encoder'](tf.convert_to_tensor(batch['x']), tf.convert_to_tensor(batch['mask']),
                            training=training)


def head_forward(model, h, sst, mode='lrt', seed=None):
    """Head pass from pooled features.

    Inputs: model, h (B,d), sst (B,1) standardised, mode, seed ([2] int for
    'sample').
    Outputs: (mu (B,3), sigma (B,3), L (B,3,3)) in standardised log space.
    """
    cfg = model['cfg']
    x = tf.concat([h, tf.cast(sst, h.dtype)], axis=1) if cfg['use_sst'] else h
    sd = None if seed is None else tf.constant(seed, tf.int32) if not tf.is_tensor(seed) else seed
    mu, s, z = model['head'](x, mode=mode, seed=sd)
    sigma = tf.nn.softplus(s) + cfg['sigma_min']
    if z is None:
        L = tf.eye(3, batch_shape=[tf.shape(mu)[0]])
    else:
        L = corr_cholesky(z, cfg['rho_max'])
    return mu, sigma, L


# ==========================================================================
# Likelihood (§6.3)
# ==========================================================================

def _log_ndtr_tf(x):
    """log Φ(x), stable for very negative x (asymptotic series below −5)."""
    x = tf.cast(x, tf.float64)
    xs = tf.where(x > -5.0, x, -5.0 * tf.ones_like(x))
    a = tf.math.log(0.5 * tf.math.erfc(-xs / np.sqrt(2.0)))
    xl = tf.where(x <= -5.0, x, -5.0 * tf.ones_like(x))
    x2 = xl * xl
    b = -0.5 * x2 - tf.math.log(-xl) - 0.5 * _LOG2PI + tf.math.log(1 - 1 / x2 + 3 / x2 ** 2 - 15 / x2 ** 3)
    return tf.where(x > -5.0, a, b)


def joint_nll_tf(lab, mu, sigma, L):
    """Per-sample negative log-likelihood of the label vector (§6.3).

    Observed taxa: multivariate Gaussian N(μ, Σ_eff) marginalised to the
    observed set; censored taxon (at most one): + −log P(y_c < t_c | y_O).
    Missing taxa are dropped. Computed in float64.

    Inputs: lab (standardise_labels dict of (B,3) tensors), mu, sigma (B,3),
    L (B,3,3).
    Outputs: (B,) NLL in standardised units.
    """
    f = lambda a: tf.cast(a, tf.float64)
    y, obs, cens, cvar = f(lab['y']), f(lab['obs']), f(lab['cens']), f(lab['cvar'])
    mu, sigma, L = f(mu), f(sigma), f(L)
    Lf = sigma[:, :, None] * L
    S = tf.matmul(Lf, Lf, transpose_b=True) + tf.linalg.diag(cvar)
    M = obs[:, :, None] * obs[:, None, :]
    Sm = S * M + tf.linalg.diag(1.0 - obs)
    r = (y - mu) * obs
    C = tf.linalg.cholesky(Sm)
    alpha = tf.linalg.cholesky_solve(C, r[:, :, None])[:, :, 0]
    nll = 0.5 * (tf.reduce_sum(r * alpha, 1)
                 + 2 * tf.reduce_sum(tf.math.log(tf.linalg.diag_part(C)), 1)
                 + tf.reduce_sum(obs, 1) * _LOG2PI)
    v = tf.linalg.matvec(S, cens) * obs
    a = tf.linalg.cholesky_solve(C, v[:, :, None])[:, :, 0]
    m_c = tf.reduce_sum(mu * cens, 1) + tf.reduce_sum(a * r, 1)
    s2 = tf.reduce_sum(cens * tf.linalg.matvec(S, cens), 1) - tf.reduce_sum(v * a, 1)
    has = tf.reduce_sum(cens, 1) > 0
    zc = (tf.reduce_sum(y * cens, 1) - m_c) / tf.sqrt(tf.maximum(s2, 1e-12))
    nll_c = tf.where(has, -_log_ndtr_tf(tf.where(has, zc, tf.zeros_like(zc))), tf.zeros_like(zc))
    return nll + nll_c


def joint_nll_np(y, observed, censored, mu, cov, cvar=None):
    """NumPy twin of :func:`joint_nll_tf` with a full covariance (any units).

    Inputs: y (n,3) (threshold where censored), observed, censored (n,3
    bool), mu (n,3), cov (n,3,3), cvar (n,3) extra diagonal variance.
    Outputs: (n,) NLL; NaN for rows with no observed or censored taxon.
    """
    n = y.shape[0]
    out = np.full(n, np.nan)
    for i in range(n):
        S = cov[i] + (np.diag(cvar[i]) if cvar is not None else 0)
        O = np.flatnonzero(observed[i])
        Cn = np.flatnonzero(censored[i])
        if O.size + Cn.size == 0:
            continue
        val = 0.0
        if O.size:
            So = S[np.ix_(O, O)]
            r = y[i, O] - mu[i, O]
            sign, ld = np.linalg.slogdet(So)
            val += 0.5 * (r @ np.linalg.solve(So, r) + ld + O.size * _LOG2PI)
        if Cn.size == 1:
            c = Cn[0]
            if O.size:
                a = np.linalg.solve(S[np.ix_(O, O)], S[O, c])
                m = mu[i, c] + a @ (y[i, O] - mu[i, O])
                s2 = S[c, c] - S[c, O] @ a
            else:
                m, s2 = mu[i, c], S[c, c]
            val -= log_ndtr((y[i, c] - m) / np.sqrt(s2))
        elif Cn.size > 1:
            continue
        out[i] = val
    return out


# ==========================================================================
# Prediction (§6.4)
# ==========================================================================

def predict_cf(model, spec, n_samples=None, seed=12345, return_draws=False, batch_size=256):
    """Predictive distribution over log10 abundances from T weight samples.

    Reports predictive mean mean(μ_t), aleatoric covariance mean(Σ_t),
    epistemic covariance cov(μ_t) and their sum; per-taxon σ from the total
    covariance, and its correlation matrix.

    Inputs: model (dict or list of dicts = seed ensemble), spec (spectra
    set), n_samples (T, default cfg), seed, return_draws, batch_size.
    Outputs: xarray.Dataset dims (sample, taxon, taxon_b[, draw]) with
    ``mu``, ``sigma``, ``sigma_aleatoric``, ``sigma_epistemic``,
    ``cov_total``, ``cov_aleatoric``, ``cov_epistemic``, ``corr_total``,
    ``corr_aleatoric``, ``pool_attention`` is not included (see
    :func:`attention_weights`); with return_draws, ``draw_mu`` (draw, sample,
    taxon) and ``draw_cov`` (draw, sample, taxon, taxon_b).
    """
    models = model if isinstance(model, (list, tuple)) else [model]
    mus, covs = [], []
    for k, m in enumerate(models):
        mu, cov = _draws(m, spec, n_samples or m['cfg']['n_samples'], seed + 100003 * k, batch_size)
        mus.append(mu)
        covs.append(cov)
    mu_t, cov_t = np.concatenate(mus), np.concatenate(covs)
    return summarise_draws(mu_t, cov_t, return_draws)


def head_fn(model, mode):
    """Compiled head pass for one weight draw (cached on the model dict).
    Inputs: model, mode ('sample' or 'mean'). Outputs: callable(h, sst, seed)."""
    key = '_head_fn_' + mode
    if key not in model:
        @tf.function(reduce_retracing=True)
        def f(h, sst, seed):
            return head_forward(model, h, sst, mode=mode, seed=seed)
        model[key] = f
    return model[key]


def _draws(model, spec, T, seed, batch_size=256):
    """(mu_t (T,n,3), cov_t (T,n,3,3)) in log10 units."""
    cfg, st = model['cfg'], model['stats']
    n = len(spec['rrs'])
    ym, ys = np.asarray(st['y_mean']), np.asarray(st['y_sd'])
    mu_t = np.zeros((T, n, 3))
    cov_t = np.zeros((T, n, 3, 3))
    for a in range(0, n, batch_size):
        rows = list(range(a, min(n, a + batch_size)))
        b = make_batch(spec, st, cfg, rows)
        h, _ = encode(model, b, training=False)
        mode = 'sample' if cfg['bayesian_head'] else 'mean'
        fn = head_fn(model, mode)
        sst = tf.constant(b['sst'])
        for t in range(T):
            if cfg['mc_dropout']:
                h, _ = encode(model, b, training=True)
            mu, sig, L = fn(h, sst, tf.constant([seed, t], tf.int32))
            mu, sig, L = mu.numpy(), sig.numpy(), L.numpy()
            Ls = (sig * ys)[:, :, None] * L
            mu_t[t, rows] = ym + ys * mu
            cov_t[t, rows] = Ls @ np.transpose(Ls, (0, 2, 1))
    return mu_t, cov_t


def summarise_draws(mu_t, cov_t, return_draws=False):
    """Moment summary of weight draws (§6.4). Inputs: mu_t (T,n,3),
    cov_t (T,n,3,3). Outputs: xarray.Dataset (see :func:`predict_cf`)."""
    T = mu_t.shape[0]
    mu = mu_t.mean(0)
    alea = cov_t.mean(0)
    d = mu_t - mu[None]
    epi = np.einsum('tni,tnj->nij', d, d) / max(T - 1, 1)
    tot = alea + epi

    def corr(c):
        s = np.sqrt(np.einsum('nii->ni', c))
        return c / (s[:, :, None] * s[:, None, :])
    ds = xr.Dataset({
        'mu': (('sample', 'taxon'), mu),
        'sigma': (('sample', 'taxon'), np.sqrt(np.einsum('nii->ni', tot))),
        'sigma_aleatoric': (('sample', 'taxon'), np.sqrt(np.einsum('nii->ni', alea))),
        'sigma_epistemic': (('sample', 'taxon'), np.sqrt(np.einsum('nii->ni', epi))),
        'cov_total': (('sample', 'taxon', 'taxon_b'), tot),
        'cov_aleatoric': (('sample', 'taxon', 'taxon_b'), alea),
        'cov_epistemic': (('sample', 'taxon', 'taxon_b'), epi),
        'corr_total': (('sample', 'taxon', 'taxon_b'), corr(tot)),
        'corr_aleatoric': (('sample', 'taxon', 'taxon_b'), corr(alea)),
    }, coords={'taxon': list(TAXA), 'taxon_b': list(TAXA)})
    if return_draws:
        ds['draw_mu'] = (('draw', 'sample', 'taxon'), mu_t)
        ds['draw_cov'] = (('draw', 'sample', 'taxon', 'taxon_b'), cov_t)
    ds.attrs['n_draws'] = T
    return ds


def attention_weights(model, spec):
    """Pooling attention weights per token (diagnostic only, §8.7: not
    attribution). Inputs: model, spec. Outputs: list of (wave, weights)."""
    b = make_batch(spec, model['stats'], model['cfg'])
    _, a = encode(model, b)
    a = a.numpy()
    return [(w, a[i, :len(w)]) for i, w in enumerate(b['wave'])]


# ==========================================================================
# Persistence
# ==========================================================================

def save_model(model, directory):
    """Write cfg + stats (JSON) and weights (npz). Inputs: model, directory."""
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, 'model.json'), 'w') as fh:
        json.dump({'cfg': model['cfg'], 'stats': model['stats'],
                   'history': model.get('history', {}), 'meta': model.get('meta', {})},
                  fh, indent=1, default=float)
    w = {f'v{i:04d}': v.numpy() for i, v in enumerate(trainable_variables(model))}
    np.savez(os.path.join(directory, 'weights.npz'), **w)
    return directory


def load_model(directory):
    """Inverse of :func:`save_model`. Inputs: directory. Outputs: model dict."""
    with open(os.path.join(directory, 'model.json')) as fh:
        d = json.load(fh)
    cfg = d['cfg']
    for k in ('wave_range_nm', 'fourier_period_nm', 'head_units'):
        cfg[k] = tuple(cfg[k])
    model = build_model(cfg, d['stats'])
    w = np.load(os.path.join(directory, 'weights.npz'))
    vs = trainable_variables(model)
    if len(vs) != len(w.files):
        raise ValueError("weight count mismatch")
    for i, v in enumerate(vs):
        v.assign(w[f'v{i:04d}'])
    model['history'], model['meta'] = d.get('history', {}), d.get('meta', {})
    return model
