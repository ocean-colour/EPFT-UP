"""Metrics, baselines, robustness and transfer tests (cf_design.md §8.1–8.6).

A *prediction* is a dict ``{'mu': (n,3) log10, 'cov': (n,3,3) total}`` with
optional ``cov_aleatoric``, ``cov_epistemic``; point-only baselines carry
``cov=None``. Labels are dicts from :func:`algorithm.labels_from_dataset`.
"""

import copy
import importlib
import os
import time

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats as sps
from scipy.special import log_ndtr

from epft_up.claude_fun import algorithm as alg
from epft_up.claude_fun import io as cfio
from epft_up.claude_fun import pipeline as pl
from epft_up.claude_fun import train as tr
from epft_up.claude_fun.io import TAXA

Z68, Z95 = sps.norm.ppf(0.84), sps.norm.ppf(0.975)
CHI2_68, CHI2_95 = sps.chi2.ppf(0.68, 3), sps.chi2.ppf(0.95, 3)
MULTISPECTRAL_NM = (412.0, 443.0, 490.0, 510.0, 555.0, 670.0)
ROBUST_KINDS = ('moana_grid', 'oci_like', 'multispectral')


def pred_from_dataset(p):
    """Prediction dict from :func:`algorithm.predict_cf` output."""
    return {'mu': p['mu'].values, 'cov': p['cov_total'].values,
            'cov_aleatoric': p['cov_aleatoric'].values, 'cov_epistemic': p['cov_epistemic'].values}


# ==========================================================================
# Metrics (§8.2)
# ==========================================================================

def point_metrics(y, mu):
    """Point skill in log10 (Lange Table 1 style plus RMSE).

    bias and MAE are reported in dex and as multiplicative factors
    10^mean(Δ) and 10^mean|Δ| (Seegers et al. 2018 convention). R² is the
    coefficient of determination; ``r2_pearson`` is the squared correlation.

    Inputs: y, mu (n,) log10 (finite pairs used).
    Outputs: dict n, bias_dex, bias, mae_dex, mae, rmse, r2, r2_pearson.
    """
    ok = np.isfinite(y) & np.isfinite(mu)
    y, mu = y[ok], mu[ok]
    if y.size == 0:
        return {'n': 0}
    d = mu - y
    sst = np.sum((y - y.mean()) ** 2)
    r = np.corrcoef(y, mu)[0, 1] if y.size > 2 and np.std(mu) > 0 else np.nan
    return {'n': int(y.size), 'bias_dex': float(d.mean()), 'bias': float(10 ** d.mean()),
            'mae_dex': float(np.abs(d).mean()), 'mae': float(10 ** np.abs(d).mean()),
            'rmse': float(np.sqrt(np.mean(d ** 2))),
            'r2': float(1 - np.sum(d ** 2) / sst) if sst > 0 else np.nan,
            'r2_pearson': float(r ** 2) if np.isfinite(r) else np.nan}


def crps_gaussian(y, mu, sd):
    """Closed-form CRPS of N(mu, sd²). Inputs/outputs: arrays."""
    z = (y - mu) / sd
    return sd * (z * (2 * sps.norm.cdf(z) - 1) + 2 * sps.norm.pdf(z) - 1 / np.sqrt(np.pi))


def energy_score(y, mu, cov, n_mc=1000, seed=0):
    """Energy score (multivariate CRPS) of N(mu, cov) by Monte Carlo.
    Inputs: y, mu (n,3), cov (n,3,3). Outputs: (n,)."""
    rng = np.random.default_rng(seed)
    out = np.full(len(y), np.nan)
    for i in range(len(y)):
        try:
            Lc = np.linalg.cholesky(cov[i])
        except np.linalg.LinAlgError:
            continue
        X = mu[i] + rng.standard_normal((n_mc, 3)) @ Lc.T
        X2 = mu[i] + rng.standard_normal((n_mc, 3)) @ Lc.T
        out[i] = (np.linalg.norm(X - y[i], axis=1).mean()
                  - 0.5 * np.linalg.norm(X - X2, axis=1).mean())
    return out


def score(pred, labels, counting_noise=True):
    """All §8.2 metrics for one prediction.

    Per taxon (marginals): point metrics on observed labels; mean NLL (censored
    labels contribute −log Φ); CRPS, 68/95 % coverage and PIT on observed
    labels. Joint: joint NLL (missing marginalised, censored conditional),
    68/95 % ellipsoid coverage (χ²₃) and energy score on fully observed rows.
    Counting noise σ_count² is added to the predictive variance of observed
    labels when ``counting_noise``.

    Inputs: pred, labels, counting_noise.
    Outputs: (table DataFrame indexed by taxon + 'joint', pit dict taxon →
    array).
    """
    y, o, c = labels['y'], labels['observed'], labels['censored']
    sc = labels['sigma_count'] if counting_noise else np.zeros_like(y)
    mu, cov = pred['mu'], pred.get('cov')
    rows, pits = [], {}
    for k, t in enumerate(TAXA):
        r = {'taxon': t, **point_metrics(np.where(o[:, k], y[:, k], np.nan), mu[:, k]),
             'n_censored': int(c[:, k].sum())}
        if cov is not None:
            sd = np.sqrt(cov[:, k, k])
            sde = np.sqrt(sd ** 2 + sc[:, k] ** 2)
            z = (y[:, k] - mu[:, k]) / sde
            nll_o = 0.5 * np.log(2 * np.pi) + np.log(sde) + 0.5 * z ** 2
            nll_c = -log_ndtr((y[:, k] - mu[:, k]) / sd)
            nll = np.where(o[:, k], nll_o, np.where(c[:, k], nll_c, np.nan))
            ok = o[:, k] & np.isfinite(z)
            r.update(nll=float(np.nanmean(nll)),
                     crps=float(np.mean(crps_gaussian(y[ok, k], mu[ok, k], sde[ok]))) if ok.any() else np.nan,
                     cov68=float(np.mean(np.abs(z[ok]) <= Z68)) if ok.any() else np.nan,
                     cov95=float(np.mean(np.abs(z[ok]) <= Z95)) if ok.any() else np.nan,
                     sigma_median=float(np.median(sd[ok])) if ok.any() else np.nan)
            if 'cov_epistemic' in pred and pred['cov_epistemic'] is not None:
                r['sigma_epistemic_median'] = float(np.median(np.sqrt(pred['cov_epistemic'][ok, k, k]))) if ok.any() else np.nan
                r['sigma_aleatoric_median'] = float(np.median(np.sqrt(pred['cov_aleatoric'][ok, k, k]))) if ok.any() else np.nan
            pits[t] = sps.norm.cdf(z[ok])
        rows.append(r)
    if cov is not None:
        j = {'taxon': 'joint'}
        cv = sc ** 2 * o
        nll = alg.joint_nll_np(np.where(o | c, y, 0.0), o, c, mu, cov, cv)
        j['nll'] = float(np.nanmean(nll))
        full = o.all(1) & np.isfinite(mu).all(1)
        if full.any():
            S = cov[full] + np.stack([np.diag(v) for v in cv[full]])
            d = y[full] - mu[full]
            m2 = np.einsum('ni,ni->n', d, np.linalg.solve(S, d[:, :, None])[:, :, 0])
            j.update(n=int(full.sum()), cov68=float(np.mean(m2 <= CHI2_68)),
                     cov95=float(np.mean(m2 <= CHI2_95)),
                     energy=float(np.nanmean(energy_score(y[full], mu[full], S))))
        rows.append(j)
    return pd.DataFrame(rows).set_index('taxon'), pits


def reliability(pred, labels, levels=np.linspace(0.1, 0.9, 9), counting_noise=True):
    """Empirical coverage of central intervals at nominal levels (reliability
    plot data). Inputs: pred, labels, levels. Outputs: DataFrame levels × taxa."""
    y, o = labels['y'], labels['observed']
    sc = labels['sigma_count'] if counting_noise else 0 * y
    out = {}
    for k, t in enumerate(TAXA):
        sd = np.sqrt(pred['cov'][:, k, k] + sc[:, k] ** 2)
        z = np.abs((y[:, k] - pred['mu'][:, k]) / sd)[o[:, k]]
        out[t] = [np.mean(z <= sps.norm.ppf(0.5 + p / 2)) for p in levels]
    return pd.DataFrame(out, index=pd.Index(levels, name='nominal'))


# ==========================================================================
# Baselines (§8.3)
# ==========================================================================

def climatology(labels_train, n):
    """Training-fold mean and SD per taxon (observed labels), independent taxa.
    Inputs: labels_train, n (rows to predict). Outputs: prediction dict."""
    y, o = labels_train['y'], labels_train['observed']
    m = np.array([y[o[:, k], k].mean() for k in range(3)])
    s = np.array([y[o[:, k], k].std(ddof=1) for k in range(3)])
    return {'mu': np.tile(m, (n, 1)), 'cov': np.tile(np.diag(s ** 2), (n, 1, 1))}


def moana_standardize(R):
    """MOANA per-spectrum standardisation (subtract mean, divide by sample SD,
    N−1), which makes the model invariant to scaling and offsets.
    Inputs: R (n, m). Outputs: Z (n, m)."""
    R = np.atleast_2d(np.asarray(R, float))
    return (R - R.mean(1, keepdims=True)) / R.std(1, ddof=1, keepdims=True)


def to_moana_grid(data_or_spec, rows=None, cfg=None):
    """Spectra on the 414–660 @ 2 nm MOANA grid (Stage 4 resampler).
    Inputs: prepare() data or spectra set; rows; pipeline cfg.
    Outputs: (n, 124) Rrs."""
    g = pl.moana_grid(cfg)
    if 'spec' in data_or_spec:
        spec = data_or_spec['spec']
    else:
        spec = data_or_spec
    rows = range(len(spec['rrs'])) if rows is None else rows
    return np.vstack([pl.resample_linear(spec['wave'][i], spec['rrs'][i], g) for i in rows])


def moana_fit(R124, labels, abundance, sst, cfg=None):
    """Retrain the MOANA recipe (standardise → PCA → linear model on PC scores
    and log10 SST) on given rows (§8.3 baseline 2).

    Pro is fitted in linear cells mL⁻¹ (zeros included, as MOANA) when
    ``pro_space='linear'``; Syn and Peuk in log10 on observed labels.

    Inputs: R124 (n,124), labels, abundance (n,3 cells mL⁻¹), sst (n,),
    cfg (DEFAULT_TRAIN['moana']).
    Outputs: model dict ``zbar``, ``V`` (k,124), ``coef`` (taxon → array),
    ``n_pc``, ``space``, ``use_sst``, ``resid_sd`` (3,) log10, ``log_clip``.
    """
    cfg = cfg or tr.DEFAULT_TRAIN['moana']
    Z = moana_standardize(R124)
    ok_rows = np.isfinite(Z).all(1)
    zbar = Z[ok_rows].mean(0)
    _, _, Vt = np.linalg.svd(Z[ok_rows] - zbar, full_matrices=False)
    model = {'zbar': zbar, 'V': Vt, 'coef': {}, 'n_pc': dict(cfg['n_pc']), 'space': {},
             'use_sst': cfg['use_sst'], 'log_clip': cfg['log_clip'], 'resid_sd': np.full(3, np.nan)}
    for k, t in enumerate(TAXA):
        X = _moana_design(model, Z, sst, t)
        if t == 'pro' and cfg['pro_space'] == 'linear':
            target = np.where(labels['censored'][:, k], 0.0, abundance[:, k])
            use = (labels['observed'][:, k] | labels['censored'][:, k]) & ok_rows
            model['space'][t] = 'linear'
        else:
            target = labels['y'][:, k]
            use = labels['observed'][:, k] & ok_rows
            model['space'][t] = 'log'
        use &= np.isfinite(X).all(1) & np.isfinite(target)
        model['coef'][t] = np.linalg.lstsq(X[use], target[use], rcond=None)[0]
        p = moana_log10(model, Z[use], sst[use], t, standardized=True)
        o = labels['observed'][use, k]
        d = p[o] - labels['y'][use, k][o]
        dof = max(o.sum() - X.shape[1], 1)
        model['resid_sd'][k] = np.sqrt(np.sum(d ** 2) / dof)
    return model


def _moana_design(model, Z, sst, taxon):
    s = (Z - model['zbar']) @ model['V'][:model['n_pc'][taxon]].T
    cols = [np.ones((Z.shape[0], 1)), s]
    if model['use_sst']:
        cols.append(np.log10(np.asarray(sst, float))[:, None])
    return np.concatenate(cols, axis=1)


def moana_raw(model, R124, sst, taxon, standardized=False):
    """MOANA output of one taxon in its fitted space (linear cells or log10).
    Inputs: model, R124 (n,124) (or Z if standardized), sst, taxon."""
    Z = R124 if standardized else moana_standardize(R124)
    return _moana_design(model, Z, sst, taxon) @ model['coef'][taxon]


def moana_log10(model, R124, sst, taxon, standardized=False):
    """MOANA output of one taxon in log10 (linear Pro floored at log_clip)."""
    v = moana_raw(model, R124, sst, taxon, standardized)
    if model['space'][taxon] == 'linear':
        return np.log10(np.maximum(v, model['log_clip']))
    return v


def moana_predict(model, R124, sst):
    """Prediction dict for retrained MOANA with a homoscedastic σ equal to the
    training residual SD (log10), so probabilistic scores can be computed.
    Inputs: model, R124, sst. Outputs: prediction dict (+ ``n_clipped``)."""
    mu = np.stack([moana_log10(model, R124, sst, t) for t in TAXA], axis=1)
    n_clip = int(np.sum(moana_raw(model, R124, sst, 'pro') < model['log_clip'])) \
        if model['space']['pro'] == 'linear' else 0
    return {'mu': mu, 'cov': np.tile(np.diag(model['resid_sd'] ** 2), (len(mu), 1, 1)),
            'n_clipped': n_clip}


def published_moana(wave, R, sst, cfg=None):
    """Published MOANA (vendored LUTs) through ``epft_up.moana`` (§8.3 baseline 3).

    The call is configurable (``DEFAULT_TRAIN['published_moana']``): module,
    function, LUT loader and output keys. Outputs are cells mL⁻¹; values below
    1 cell mL⁻¹ (MOANA clamps negative Pro to 0) are floored there before log10.

    Inputs: wave (m,), R (n,m) on that grid, sst (n,), cfg.
    Outputs: prediction dict with ``cov=None`` (point-only) and ``n_clipped``.
    """
    cfg = cfg or tr.DEFAULT_TRAIN['published_moana']
    mod = importlib.import_module(cfg['module'])
    lut = getattr(importlib.import_module(cfg['lut_module']), cfg['lut_func'])()
    out = getattr(mod, cfg['func'])(np.asarray(wave), np.asarray(R), sst=np.asarray(sst), lut=lut)
    cols, clip = [], 0
    for t in TAXA:
        key = cfg['keys'][t]
        v = np.asarray(out[key] if isinstance(out, dict) else getattr(out, key), float)
        clip += int(np.sum(v < 1.0))
        cols.append(np.log10(np.maximum(v, 1.0)))
    return {'mu': np.stack(cols, 1), 'cov': None, 'n_clipped': clip}


# ==========================================================================
# Degraded inputs (§8.4)
# ==========================================================================

def degrade_spec(spec, kind, pipeline_cfg=None):
    """Re-grid a spectra set for the resolution-robustness test.

    kinds: ``moana_grid`` (linear to 414–660 @ 2 nm), ``oci_like`` (Gaussian
    5 nm FWHM, 2.5 nm sampling), ``multispectral`` (linear at 412, 443, 490,
    510, 555, 670 nm). σ is propagated (independent bands) and the σ floor
    re-applied.

    Inputs: spec, kind, pipeline_cfg. Outputs: spectra set.
    """
    pc = pipeline_cfg or pl.DEFAULT_PIPELINE
    fl = (pc['sigma_f_rel'], pc['sigma_abs'])
    out = {'wave': [], 'rrs': [], 'sigma': [], 'sst': np.asarray(spec['sst'])}
    for w, r, s in zip(spec['wave'], spec['rrs'], spec['sigma']):
        if kind == 'moana_grid':
            g = pl.moana_grid(pc)
            r2, s2 = pl.resample_linear(w, r, g, s)
            r2, s2, g2 = r2[0], s2[0], g
        elif kind == 'multispectral':
            g = np.array(MULTISPECTRAL_NM)
            r2, s2 = pl.resample_linear(w, r, g, s)
            r2, s2, g2 = r2[0], s2[0], g
        elif kind == 'oci_like':
            g = np.arange(np.ceil(w[0] / 2.5) * 2.5, w[-1] + 1e-9, 2.5)
            g2, r2, s2 = tr.gaussian_srf_resample(w, r, s, g, 5.0)
        else:
            raise ValueError(kind)
        ok = np.isfinite(r2) & np.isfinite(s2)
        out['wave'].append(g2[ok])
        out['rrs'].append(r2[ok])
        out['sigma'].append(np.fmax(s2[ok], tr._floor(r2[ok], fl)))
    return out


def with_sst(spec, sst):
    """Copy of a spectra set with SST replaced. Inputs: spec, sst (n,)."""
    out = dict(spec)
    out['sst'] = np.asarray(sst, float)
    return out


# ==========================================================================
# Cross-validation runner (§8.1–8.4)
# ==========================================================================

def _empty(n):
    return {'mu': np.full((n, 3), np.nan), 'cov': np.full((n, 3, 3), np.nan),
            'cov_aleatoric': np.full((n, 3, 3), np.nan), 'cov_epistemic': np.full((n, 3, 3), np.nan)}


def _put(dst, src, rows):
    for k in dst:
        if src.get(k) is not None:
            dst[k][rows] = src[k]


def run_cv(data, folds, model_cfg=None, train_cfg=None, pipeline_cfg=None, robustness=True,
           sst_perm=True, baselines=True, attribution_fn=None, n_samples=None, verbose=True,
           label='cv'):
    """Cross-validate Claude Fun and the baselines on identical folds.

    Every fitted statistic (s_R, label and SST standardisation, MOANA PCA) is
    computed inside each training fold.

    Inputs
    ------
    data : train.prepare output.
    folds : list of (train_rows, test_rows).
    robustness, sst_perm, baselines : bool — extra predictions per fold.
    attribution_fn : callable(models, data, train_rows, test_rows) → object,
        collected per fold (see attribution.fold_maps).
    n_samples : weight draws per prediction (default model cfg).

    Outputs
    -------
    dict — ``cf``, ``climatology``, ``moana`` (prediction dicts over all
    samples), ``robust`` (kind → prediction), ``sst_perm`` (list of
    predictions), ``fold`` (n,), ``attribution`` (list), ``models_meta``
    (list), ``labels``.
    """
    mc = model_cfg or alg.DEFAULT_MODEL
    tc = train_cfg or tr.DEFAULT_TRAIN
    pc = pipeline_cfg or pl.DEFAULT_PIPELINE
    n = data['n']
    res = {'cf': _empty(n), 'climatology': _empty(n), 'moana': _empty(n),
           'robust': {k: _empty(n) for k in ROBUST_KINDS} if robustness else {},
           'sst_perm': [_empty(n) for _ in range(tc['sst_permutations'])] if sst_perm else [],
           'fold': np.full(n, -1), 'attribution': [], 'models_meta': [], 'labels': data['labels'],
           'moana_clipped': 0}
    rng = np.random.default_rng(tc['seed'] + 7)
    perms = [rng.permutation(n) for _ in res['sst_perm']]
    sst_all = np.asarray(data['spec']['sst'])
    R124 = to_moana_grid(data, cfg=pc) if baselines else None
    ab = data['ds']['abundance'].values
    for f, (trn, tst) in enumerate(folds):
        t0 = time.time()
        models = tr.train_ensemble(data, trn, mc, tc)
        st = alg.subset(data['spec'], tst)
        _put(res['cf'], pred_from_dataset(alg.predict_cf(models, st, n_samples)), tst)
        res['fold'][tst] = f
        for r, perm in enumerate(perms):
            p = alg.predict_cf(models, with_sst(st, sst_all[perm[tst]]), n_samples)
            _put(res['sst_perm'][r], pred_from_dataset(p), tst)
        for kind in res['robust']:
            p = alg.predict_cf(models, degrade_spec(st, kind, pc), n_samples)
            _put(res['robust'][kind], pred_from_dataset(p), tst)
        if baselines:
            _put(res['climatology'], climatology(tr.sub_labels(data['labels'], trn), len(tst)), tst)
            mm = moana_fit(R124[trn], tr.sub_labels(data['labels'], trn), ab[trn], sst_all[trn], tc['moana'])
            pm = moana_predict(mm, R124[tst], sst_all[tst])
            res['moana_clipped'] += pm['n_clipped']
            _put(res['moana'], pm, tst)
        if attribution_fn is not None:
            res['attribution'].append(attribution_fn(models, data, trn, tst))
        res['models_meta'].append([m['meta'] for m in models])
        if verbose:
            mm_ = models[0]['meta']
            print(f"  [{label}] fold {f + 1}/{len(folds)}: n_train {len(trn)}, test {len(tst)}, "
                  f"best epoch {mm_['best_epoch']}/{mm_['epochs_run']}, {time.time() - t0:.0f} s", flush=True)
    return res


def robustness_table(native, robust, labels, counting_noise=True):
    """|Δμ|/σ_total relative to the native-grid prediction and the change in
    skill, per degradation kind and taxon (§8.4; pass: median < 0.25 for
    moana_grid and oci_like).

    Inputs: native (prediction), robust (kind → prediction), labels.
    Outputs: DataFrame (kind, taxon) → median_shift, p84_shift, d_mae_dex,
    d_rmse, d_nll, pass.
    """
    base, _ = score(native, labels, counting_noise)
    rows = []
    
    # Helper to safely compute metric differences without crashing on missing columns
    def safe_diff(s_df, base_df, taxon, col):
        try:
            if col in s_df.columns and col in base_df.columns:
                val_s = s_df.loc[taxon, col]
                val_base = base_df.loc[taxon, col]
                if pd.notna(val_s) and pd.notna(val_base):
                    return float(val_s - val_base)
        except (KeyError, TypeError):
            pass
        return np.nan

    for kind, p in robust.items():
        s, _ = score(p, labels, counting_noise)
        for k, t in enumerate(TAXA):
            sh = np.abs(p['mu'][:, k] - native['mu'][:, k]) / np.sqrt(native['cov'][:, k, k])
            sh = sh[np.isfinite(sh)]
            med = float(np.median(sh)) if sh.size else np.nan
            
            rows.append({
                'kind': kind, 
                'taxon': t, 
                'median_shift': med,
                'p84_shift': float(np.percentile(sh, 84)) if sh.size else np.nan,
                'd_mae_dex': safe_diff(s, base, t, 'mae_dex'),
                'd_rmse': safe_diff(s, base, t, 'rmse'),
                'd_nll': safe_diff(s, base, t, 'nll'),
                'pass': (med < 0.25) if kind in ('moana_grid', 'oci_like') else np.nan
            })
            
    return pd.DataFrame(rows).set_index(['kind', 'taxon'])


def sst_reliance_table(native, permuted, labels, counting_noise=True):
    """Permutation importance of SST: ΔNLL and ΔMAE (mean over repeats)
    when SST is shuffled across held-out samples (§8.2).
    Inputs: native, permuted (list of predictions), labels.
    Outputs: DataFrame per taxon (+ joint NLL): d_nll, d_mae_dex, sd over repeats."""
    base, _ = score(native, labels, counting_noise)
    d = []
    for p in permuted:
        s, _ = score(p, labels, counting_noise)
        d.append(pd.DataFrame({'d_nll': s['nll'] - base['nll'],
                               'd_mae_dex': s.get('mae_dex', np.nan) - base.get('mae_dex', np.nan)}))
    D = pd.concat(d, keys=range(len(d)))
    out = D.groupby(level=1).agg(['mean', 'std'])
    out.columns = ['_'.join(c) for c in out.columns]
    return out


def cv_report(res, counting_noise=True):
    """Metric tables for a CV result.
    Inputs: run_cv output. Outputs: dict of DataFrames ``metrics`` (model,
    taxon), ``robustness``, ``sst_reliance``, ``reliability`` (CF), ``pit``."""
    lab = res['labels']
    tabs, pit = [], None
    for name in ('cf', 'moana', 'climatology'):
        if np.isfinite(res[name]['mu']).any():
            s, pits = score(res[name], lab, counting_noise)
            tabs.append(pd.concat({name: s}, names=['model']))
            if name == 'cf':
                pit = pits
    out = {'metrics': pd.concat(tabs), 'pit': pit,
           'reliability': reliability(res['cf'], lab, counting_noise=counting_noise)}
    if res['robust']:
        out['robustness'] = robustness_table(res['cf'], res['robust'], lab, counting_noise)
    if res['sst_perm']:
        out['sst_reliance'] = sst_reliance_table(res['cf'], res['sst_perm'], lab, counting_noise)
    return out


def predictions_to_dataset(res, data):
    """Held-out predictions of a CV/transfer result as an xarray.Dataset
    (for the report script). Inputs: result dict, data. Outputs: Dataset."""
    ds = xr.Dataset(coords={'sample': np.arange(data['n']), 'taxon': list(TAXA),
                            'taxon_b': list(TAXA)})
    for name in ('cf', 'moana', 'climatology', 'published_moana'):
        if name in res and res[name] is not None:
            ds[f'{name}_mu'] = (('sample', 'taxon'), res[name]['mu'])
            if res[name].get('cov') is not None:
                ds[f'{name}_cov'] = (('sample', 'taxon', 'taxon_b'), res[name]['cov'])
    for k in ('cov_aleatoric', 'cov_epistemic'):
        ds[f'cf_{k}'] = (('sample', 'taxon', 'taxon_b'), res['cf'][k])
    for kind, p in res.get('robust', {}).items():
        ds[f'robust_{kind}_mu'] = (('sample', 'taxon'), p['mu'])
        ds[f'robust_{kind}_cov'] = (('sample', 'taxon', 'taxon_b'), p['cov'])
    lab = data['labels']
    ds['y_log10'] = (('sample', 'taxon'), lab['y'])
    ds['observed'] = (('sample', 'taxon'), lab['observed'])
    ds['censored'] = (('sample', 'taxon'), lab['censored'])
    ds['sigma_count'] = (('sample', 'taxon'), lab['sigma_count'])
    if 'fold' in res:
        ds['fold'] = ('sample', res['fold'])
    for v in ('lat', 'lon', 'sst', 'cruise', 'station', 'time'):
        if v in data['ds']:
            ds[v] = ('sample', data['ds'][v].values)
    return ds


# ==========================================================================
# Transfer (§8.5)
# ==========================================================================

def transfer_test(train_data, test_data, model_cfg=None, train_cfg=None, pipeline_cfg=None,
                  models=None, n_samples=None, robustness=True, sst_perm=True, published=True,
                  verbose=True):
    """Train on all AMT24 labels and evaluate on held-out cruises.

    Compares Claude Fun with AMT24-retrained MOANA, published MOANA and the
    AMT24 climatology; everything is out-of-sample for every model.

    Inputs: train_data, test_data (train.prepare), configs, models (reuse a
    trained list), n_samples, robustness, sst_perm, published.
    Outputs: dict like :func:`run_cv` (no folds) plus ``published_moana``,
    ``models`` and ``per_cruise`` metric tables.
    """
    mc = model_cfg or alg.DEFAULT_MODEL
    tc = train_cfg or tr.DEFAULT_TRAIN
    pc = pipeline_cfg or pl.DEFAULT_PIPELINE
    if models is None:
        models = tr.train_ensemble(train_data, None, mc, tc, verbose=verbose)
    st = test_data['spec']
    n = test_data['n']
    res = {'cf': pred_from_dataset(alg.predict_cf(models, st, n_samples)), 'labels': test_data['labels'],
           'models': models, 'robust': {}, 'sst_perm': []}
    sst_tr = np.asarray(train_data['spec']['sst'])
    R_tr, R_te = to_moana_grid(train_data, cfg=pc), to_moana_grid(test_data, cfg=pc)
    mm = moana_fit(R_tr, train_data['labels'], train_data['ds']['abundance'].values, sst_tr, tc['moana'])
    res['moana'] = moana_predict(mm, R_te, np.asarray(st['sst']))
    res['climatology'] = climatology(train_data['labels'], n)
    res['published_moana'] = None
    if published:
        try:
            res['published_moana'] = published_moana(pl.moana_grid(pc), R_te, np.asarray(st['sst']),
                                                      tc['published_moana'])
        except Exception as e:      # the in-repo package may not be importable here
            res['published_moana_error'] = f"{type(e).__name__}: {e}"
            if verbose:
                print(f"  published MOANA unavailable: {res['published_moana_error']}")
    if sst_perm:
        rng = np.random.default_rng(tc['seed'] + 11)
        for _ in range(tc['sst_permutations']):
            res['sst_perm'].append(pred_from_dataset(
                alg.predict_cf(models, with_sst(st, np.asarray(st['sst'])[rng.permutation(n)]), n_samples)))
    if robustness:
        for kind in ROBUST_KINDS:
            res['robust'][kind] = pred_from_dataset(alg.predict_cf(models, degrade_spec(st, kind, pc), n_samples))
    cruises = test_data['ds']['cruise'].values.astype(str)
    per = []
    for cr in list(np.unique(cruises)) + ['all']:
        rows = np.arange(n) if cr == 'all' else np.flatnonzero(cruises == cr)
        lab = tr.sub_labels(test_data['labels'], rows)
        for name in ('cf', 'moana', 'published_moana', 'climatology'):
            p = res.get(name)
            if p is None:
                continue
            sub = {k: (v[rows] if isinstance(v, np.ndarray) else v) for k, v in p.items()}
            s, _ = score(sub, lab)
            per.append(pd.concat({(cr, name): s}, names=['cruise', 'model']))
    res['per_cruise'] = pd.concat(per)
    return res


# ==========================================================================
# SeaFlow + PACE (§8.6, test only; minimal)
# ==========================================================================

def _retry(fn, *a, tries=3, wait=5.0, **k):
    for i in range(tries):
        try:
            return fn(*a, **k)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(wait * 2 ** i)


def seaflow_pace_matchups(seaflow, pipeline_cfg=None, download_dir='pace_l2', short_name='PACE_OCI_L2_AOP',
                          max_days=None, verbose=True):
    """Match SeaFlow samples to PACE OCI L2 Rrs (§8.6).

    For each PACE granule overlapping a SeaFlow day, the SeaFlow sample
    nearest in time within ±3 h is matched (one matchup per overpass, to
    avoid replicating one pixel box across many 3-min samples). The 5×5 box
    around the ship position must be ≥ half valid after the L2 flags; the
    spectrum is the box median and σ_Rrs = max(box spread (p84−p16)/2, median
    Rrs_unc, σ_floor). SST is the SeaFlow underway value. Non-positive
    abundances are treated as missing (no SeaFlow detection limit is set).

    Inputs: seaflow (io.read_seaflow DataFrame), pipeline_cfg, download_dir,
    short_name, max_days, verbose.
    Outputs: xarray.Dataset in the transfer-set layout (cruise = 'SeaFlow'),
    or None if nothing matched.
    """
    import earthaccess
    pc = pipeline_cfg or pl.DEFAULT_PIPELINE
    earthaccess.login()
    ds_ = _retry(earthaccess.search_datasets, short_name=short_name, cloud_hosted=True) \
        or _retry(earthaccess.search_datasets, short_name=short_name)
    concept = ds_[0]['meta']['concept-id']
    sf = seaflow.dropna(subset=['time', 'lat', 'lon']).sort_values('time').reset_index(drop=True)
    days = sf['time'].dt.floor('D').unique()[:max_days]
    out = []
    for d in days:
        day = sf[sf['time'].dt.floor('D') == d]
        bbox = (day['lon'].min() - 1, day['lat'].min() - 1, day['lon'].max() + 1, day['lat'].max() + 1)
        t0 = (d - pd.Timedelta(hours=pc['pace_max_dt_h'])).isoformat()
        t1 = (d + pd.Timedelta(days=1, hours=pc['pace_max_dt_h'])).isoformat()
        res = _retry(earthaccess.search_data, concept_id=concept, temporal=(t0, t1), bounding_box=bbox)
        if not res:
            continue
        files = _retry(earthaccess.download, res, local_path=download_dir)
        for f in files:
            try:
                with xr.open_dataset(f) as root:
                    tg = pd.Timestamp(root.attrs['time_coverage_start']).tz_localize(None)
            except Exception:
                continue
            dt = (day['time'] - tg).dt.total_seconds().abs() / 3600
            if not (dt <= pc['pace_max_dt_h']).any():
                continue
            j = dt.idxmin()
            s = sf.loc[j]
            box = cfio.read_pace_box(f, s['lat'], s['lon'], pc['pace_box_half'])
            if box is None or box['valid_fraction'] < pc['pace_min_valid_frac']:
                continue
            R = box['rrs'][box['valid']]
            med = np.nanmedian(R, 0)
            p16, p84 = np.nanpercentile(R, [16, 84], axis=0)
            unc = np.nanmedian(box['rrs_unc'][box['valid']], 0)
            sig = np.fmax(np.fmax((p84 - p16) / 2, np.nan_to_num(unc)), pl.sigma_floor(med, pc))
            out.append({'wave': box['wave'], 'rrs': med, 'sigma': sig, 'sst': s['sst'], 'lat': s['lat'],
                        'lon': s['lon'], 'time': s['time'], 'dt_h': float(dt[j]), 'granule': os.path.basename(f),
                        'valid_fraction': box['valid_fraction'],
                        **{t: s[t] for t in TAXA}})
        if verbose:
            print(f"  SeaFlow/PACE {pd.Timestamp(d).date()}: {len(out)} matchups so far", flush=True)
    if not out:
        return None
    wave = out[0]['wave']
    n = len(out)
    ab = np.array([[o[t] for t in TAXA] for o in out], float)
    obs = np.isfinite(ab) & (ab > 0)
    st = ('sample', 'taxon')
    with np.errstate(divide='ignore', invalid='ignore'):
        y = np.where(obs, np.log10(ab), np.nan)
    ds = xr.Dataset({
        'rrs': (('sample', 'wavelength'), np.stack([o['rrs'] for o in out])),
        'sigma': (('sample', 'wavelength'), np.stack([o['sigma'] for o in out])),
        'sst': ('sample', np.array([o['sst'] for o in out], float)),
        'lat': ('sample', np.array([o['lat'] for o in out], float)),
        'lon': ('sample', np.array([o['lon'] for o in out], float)),
        'time': ('sample', pd.to_datetime([o['time'] for o in out]).values),
        'dt_h': ('sample', np.array([o['dt_h'] for o in out])),
        'granule': ('sample', np.array([o['granule'] for o in out])),
        'cruise': ('sample', np.full(n, 'SeaFlow')),
        'station': ('sample', np.array([str(i) for i in range(n)])),
        'abundance': (st, ab), 'y_log10': (st, y), 'observed': (st, obs),
        'censored': (st, np.zeros((n, 3), bool)), 'missing': (st, ~obs),
        'dl': (st, np.full((n, 3), np.nan)), 'n_events': (st, np.zeros((n, 3))),
        'sigma_count': (st, np.zeros((n, 3))),
    }, coords={'wavelength': wave, 'taxon': list(TAXA)})
    ds.attrs.update(sst_source='seaflow_underway', grid='pace_oci')
    return ds


def load_lange_reference(path=None):
    """Lange et al. 2020 Table 1/2 numbers transcribed to
    ``reference/lange2020_tables.csv`` (reference only, in-sample).
    Inputs: path. Outputs: DataFrame (empty values allowed)."""
    path = path or os.path.join(os.path.dirname(__file__), 'reference', 'lange2020_tables.csv')
    return pd.read_csv(path, comment='#')
