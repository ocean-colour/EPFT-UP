"""Validation schemes, baselines and scores for the SDP pigment models.

Three questions (Execution #4):

1. **Skill beyond Tchla.** Accessory pigments co-vary with Tchla (Kramer &
   Siegel 2019), so a model can score well on absolute concentrations just by
   retrieving Tchla. The *Tchla null model* predicts each pigment from a
   retrieved Tchla alone, via an OLS fit in log10 space trained on the same
   split: log10 P = a + b log10 Tchla_ret. ``Tchla_ret`` comes from the GSM
   fit (Execution #2) or from the OC4v6 band ratio (multispectral-like). An
   *oracle* version uses HPLC Tchla: the most that Tchla covariance alone
   can give, which no retrieval can reach.
2. **Composition.** Pigment:Tchla ratios remove the shared biomass signal.
   They are scored in linear space and in log10 space. Zeros (below LOD)
   are replaced through one swappable function, :func:`replace_zeros`.
3. **Leakage.** Kramer's random 75/25 splits mix samples of the same cruise
   between training and validation. *Leave-one-campaign-out* (LOCO) holds
   out each of the 8 campaigns in turn.

Every model is a *fitter*: ``fitter(train_idx, test_idx) -> predictions``.
:func:`cross_validate` runs a fitter over a list of splits, so PCR and every
baseline see identical splits and can be compared pairwise.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from epft_up.sdp import models

#: Reporting resolution of the deposited pigments [mg m^-3].
PIGMENT_RESOLUTION = 0.001
#: OC4v6 (O'Reilly et al. 1998; SeaWiFS v6) polynomial, max(443,490,510)/555.
OC4V6 = (0.3272, -2.9940, 2.7218, -1.2259, -0.5683)


# --------------------------------------------------------------------------
# Detection limits and zero handling
# --------------------------------------------------------------------------
def lod_proxy(pigments):
    """Per-pigment detection-limit proxy [mg m^-3] (report §4.2).

    For a pigment with below-LOD zeros in the deposit, use the smallest
    non-zero value (an upper bound on the effective LOD). For one that is
    never below LOD (Tchla, Zea, Chlc12) the smallest value is only the sample
    minimum, so use the 0.001 reporting resolution.

    Parameters
    ----------
    pigments : DataFrame

    Returns
    -------
    dict
    """
    out = {}
    for p in pigments.columns:
        v = pigments[p]
        out[p] = float(v[v > 0].min()) if (v == 0).any() else PIGMENT_RESOLUTION
    return out


def replace_zeros(x, lod, frac=0.5):
    """Replace values ≤ 0 by ``frac * lod`` (multiplicative replacement; Q&A #20).

    This is the single place where below-detection values become positive
    numbers. A censored-likelihood treatment (Tom Jordan's log-normal model)
    can replace this function later.

    Parameters
    ----------
    x : ndarray
    lod : float
    frac : float, optional

    Returns
    -------
    ndarray
    """
    x = np.asarray(x, dtype=float)
    return np.where(x > 0, x, frac * lod)


# --------------------------------------------------------------------------
# A band-ratio Tchla for the null model
# --------------------------------------------------------------------------
def oc4v6(wave, Rrs):
    """OC4v6 Tchla [mg m^-3] from hyperspectral Rrs at 443/490/510/555 nm."""
    wave = np.asarray(wave)
    Rrs = np.atleast_2d(Rrs)

    def at(w):
        return Rrs[:, int(np.searchsorted(wave, w))]

    x = np.log10(np.maximum.reduce([at(443), at(490), at(510)]) / at(555))
    return 10.0**np.polyval(OC4V6[::-1], x)


# --------------------------------------------------------------------------
# Splits
# --------------------------------------------------------------------------
def random_splits(n, n_perm=100, train_frac=0.75, seed=1):
    """Kramer-style random splits: ``n_perm`` × (floor(train_frac n) train, rest test).

    Returns
    -------
    list of (ndarray, ndarray)
    """
    rng = np.random.default_rng(seed)
    n_tr = int(np.floor(n * train_frac))
    out = []
    for _ in range(n_perm):
        tr = rng.permutation(n)[:n_tr]
        out.append((np.sort(tr), np.setdiff1d(np.arange(n), tr)))
    return out


def loco_splits(groups):
    """Leave-one-group-out splits, one per distinct value of ``groups``.

    Returns
    -------
    names : list
    splits : list of (ndarray, ndarray)
    """
    groups = np.asarray(groups)
    names = list(dict.fromkeys(groups))
    idx = np.arange(len(groups))
    return names, [(idx[groups != g], idx[groups == g]) for g in names]


# --------------------------------------------------------------------------
# Fitters
# --------------------------------------------------------------------------
class PCRFitter:
    """Kramer's PCR (``models.fit_pcr_one``) as a fitter.

    Parameters
    ----------
    X : ndarray, shape (n, p)
    y : ndarray, shape (n,)
    constraint : {'pigment', None}
        ``'pigment'`` clips at 0 (concentrations); ``None`` for log targets.
    seed : int
        Seeds the inner k-fold partitions (one generator over all splits).
    """

    def __init__(self, X, y, constraint='pigment', seed=1, **kw):
        self.X, self.y = np.asarray(X, float), np.asarray(y, float)
        self.constraint = constraint
        self.rng = np.random.default_rng(seed)
        self.kw = kw

    def __call__(self, tr, te):
        coef, icpt, _ = models.fit_pcr_one(self.X[tr], self.y[tr], self.rng,
                                           constraint=self.constraint, **self.kw)
        pred = self.X[te] @ coef + icpt
        return np.maximum(pred, 0.0) if self.constraint == 'pigment' else pred


class LogLinearFitter:
    """OLS of ``log10 y`` on ``log10 x`` (or a supplied transform), back-transformed.

    The Tchla null model when ``x`` is a retrieved Tchla.

    Parameters
    ----------
    x : ndarray, shape (n,)
        Predictor (> 0), e.g. GSM or OC4 Tchla.
    y_log : ndarray, shape (n,)
        log10 of the target (zeros already replaced).
    output : {'linear', 'log'}
        Return ``10**prediction`` or the log10 prediction.
    """

    def __init__(self, x, y_log, output='linear'):
        x = np.asarray(x, float)
        if np.any(~np.isfinite(x)) or np.any(x <= 0):
            raise ValueError('LogLinearFitter needs a finite, positive predictor')
        self.lx = np.log10(x)
        self.ly = np.asarray(y_log, float)
        self.output = output

    def __call__(self, tr, te):
        b, a = np.polyfit(self.lx[tr], self.ly[tr], 1)
        pred = a + b * self.lx[te]
        return 10.0**pred if self.output == 'linear' else pred


class ConstantFitter:
    """Predict the training mean (no-skill reference)."""

    def __init__(self, y):
        self.y = np.asarray(y, float)

    def __call__(self, tr, te):
        return np.full(len(te), self.y[tr].mean())


def cross_validate(fitter, splits):
    """Run ``fitter`` on every split.

    Returns
    -------
    list of (test_idx, prediction)
    """
    return [(te, np.asarray(fitter(tr, te), float)) for tr, te in splits]


# --------------------------------------------------------------------------
# Scores
# --------------------------------------------------------------------------
def score_linear(obs, pred):
    """Table 2-style scores in concentration units.

    ``R2`` and ``RMSE`` come from an OLS of observed on modelled (as
    ``rrsModelTrain.m``); ``MAE``; ``MADn`` = MAE / mean modelled (Table 2's
    "normalized MAD"); ``SS_MAE`` is filled in by :func:`paired`.
    """
    obs = np.asarray(obs, float)
    pred = np.asarray(pred, float)
    r2, rmse = models._ols_r2_rmse(pred, obs)
    mae = float(np.mean(np.abs(pred - obs)))
    mp = float(np.mean(pred))
    return {'R2': r2, 'RMSE': rmse, 'MAE': mae,
            'MADn': mae / mp if mp > 0 else np.nan}


def score_log(obs_log, pred_log):
    """Scores in log10 space: R² (squared correlation), OLS slope, RMS and bias."""
    o = np.asarray(obs_log, float)
    p = np.asarray(pred_log, float)
    r = np.corrcoef(o, p)[0, 1] if np.std(p) > 0 and np.std(o) > 0 else 0.0
    slope = np.polyfit(o, p, 1)[0] if np.std(o) > 0 else np.nan
    return {'R2': float(r**2), 'slope': float(slope),
            'RMS': float(np.sqrt(np.mean((p - o)**2))),
            'bias': float(np.mean(p - o)),
            'MAE': float(np.mean(np.abs(p - o)))}


def per_split_scores(results, obs, scorer):
    """Score each split separately (random-split scheme).

    Returns
    -------
    dict of ndarray
    """
    rows = [scorer(obs[te], pred) for te, pred in results]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def pooled_predictions(results, n):
    """Assemble one prediction per sample from a partition (LOCO).

    Raises if a sample is predicted zero times or more than once.
    """
    out = np.full(n, np.nan)
    seen = np.zeros(n, dtype=int)
    for te, pred in results:
        out[te] = pred
        seen[te] += 1
    if not np.all(seen == 1):
        raise ValueError('splits are not a partition of the samples')
    return out


def paired(model_scores, null_scores, higher_is_better=('R2',)):
    """Paired comparison of two per-split score dicts on identical splits.

    Returns, per score: mean difference (model − null), the fraction of splits
    in which the model is better, and the MAE skill score
    ``1 − MAE_model / MAE_null`` (mean over splits).
    """
    out = {}
    for k in model_scores:
        if k not in null_scores:
            continue
        d = model_scores[k] - null_scores[k]
        better = d > 0 if k in higher_is_better else d < 0
        out[k] = {'mean_diff': float(np.nanmean(d)),
                  'frac_model_better': float(np.mean(better))}
    if 'MAE' in model_scores:
        out['SS_MAE'] = float(np.nanmean(1.0 - model_scores['MAE'] / null_scores['MAE']))
    return out


def bootstrap_paired(obs, pred_model, pred_null, scorer, n_boot=2000, seed=0):
    """Bootstrap (over samples) the paired R² and MAE differences of pooled predictions.

    Returns
    -------
    dict
        ``dR2`` and ``SS_MAE`` (point estimates) with 2.5/97.5 percentiles and
        the fraction of resamples in which the model beats the null on R².
    """
    rng = np.random.default_rng(seed)
    n = len(obs)
    d_r2, ss = [], []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        sm, sn = scorer(obs[i], pred_model[i]), scorer(obs[i], pred_null[i])
        d_r2.append(sm['R2'] - sn['R2'])
        ss.append(1.0 - sm['MAE'] / sn['MAE'])
    sm, sn = scorer(obs, pred_model), scorer(obs, pred_null)
    d_r2, ss = np.array(d_r2), np.array(ss)
    return {'dR2': float(sm['R2'] - sn['R2']),
            'dR2_ci': np.percentile(d_r2, [2.5, 97.5]).tolist(),
            'SS_MAE': float(1.0 - sm['MAE'] / sn['MAE']),
            'SS_MAE_ci': np.percentile(ss, [2.5, 97.5]).tolist(),
            'frac_boot_R2_better': float(np.mean(d_r2 > 0))}


# --------------------------------------------------------------------------
# Benchmark: the standard scorecard every model is judged on (from #4 on)
# --------------------------------------------------------------------------
class Benchmark:
    """Fixed targets, splits and baselines against which any model is scored.

    Every candidate is judged on two criteria at once (user instruction,
    2026-10-02): skill **above the Tchla-only null models**, and skill under
    **leave-one-campaign-out** (LOCO) as well as Kramer's random splits.

    Targets:

    * ``abs:<P>``: concentration of each of the 13 pigments, scored with the
      Table 2 linear R² (OLS of observed on modelled) and log10 R²;
    * ``ratio:<P>``: log10(P/Tchla) for the 12 accessory pigments (zeros →
      ``zero_frac`` × LOD), scored with log10 R² and RMS.

    Baselines (computed once and cached): ``null_GSM`` and ``null_OC4``
    (two-parameter log-log fits on retrieved Tchla, trained per split) and,
    for ratios, ``const`` (the training mean).

    Parameters
    ----------
    pigments : DataFrame
        The 13 pigments (columns as ``data.PIGMENTS_13``), one row per sample.
    tchla_gsm, tchla_oc4 : ndarray
        Retrieved Tchla for the null models.
    campaigns : ndarray
        Campaign label per sample (LOCO groups).
    n_perm, seed, zero_frac : optional
    """

    NULLS = ('null_GSM', 'null_OC4')

    def __init__(self, pigments, tchla_gsm, tchla_oc4, campaigns, n_perm=100, seed=1,
                 zero_frac=0.5):
        self.pig = pigments.reset_index(drop=True)
        self.n = len(self.pig)
        self.lods = lod_proxy(self.pig)
        self.tchla_gsm = np.asarray(tchla_gsm, float)
        self.tchla_oc4 = np.asarray(tchla_oc4, float)
        self.zero_frac = zero_frac
        self.seed = seed
        names, loco = loco_splits(np.asarray(campaigns))
        self.loco_names = names
        self.splits = {'random': random_splits(self.n, n_perm, 0.75, seed), 'loco': loco}
        tchla = self.pig['Tchla'].to_numpy()
        self.targets = {}
        for p in self.pig.columns:
            y = self.pig[p].to_numpy()
            self.targets[f'abs:{p}'] = {
                'kind': 'abs', 'pigment': p, 'y': y,
                'ylog': np.log10(replace_zeros(y, self.lods[p], zero_frac))}
            if p != 'Tchla':
                rl = np.log10(replace_zeros(y, self.lods[p], zero_frac) / tchla)
                self.targets[f'ratio:{p}'] = {'kind': 'ratio', 'pigment': p,
                                              'y': rl, 'ylog': rl}
        self._baseline = {t: self._run_baselines(t) for t in self.targets}

    # -- internals ---------------------------------------------------------
    def _baseline_fitters(self, t):
        tg = self.targets[t]
        out = 'linear' if tg['kind'] == 'abs' else 'log'
        f = {'null_GSM': LogLinearFitter(self.tchla_gsm, tg['ylog'], output=out),
             'null_OC4': LogLinearFitter(self.tchla_oc4, tg['ylog'], output=out)}
        if tg['kind'] == 'ratio':
            f['const'] = ConstantFitter(tg['y'])
        return f

    def _run_baselines(self, t):
        return {scheme: {m: cross_validate(f, sp) for m, f in self._baseline_fitters(t).items()}
                for scheme, sp in self.splits.items()}

    def _to_log(self, t, results):
        tg = self.targets[t]
        if tg['kind'] == 'ratio':
            return results
        lod = self.lods[tg['pigment']]
        return [(te, np.log10(replace_zeros(pr, lod, self.zero_frac))) for te, pr in results]

    def _score_scheme(self, t, scheme, results_by_model):
        """Scores of every model for one target and scheme (random: mean ± SD)."""
        tg = self.targets[t]
        primary = []
        if tg['kind'] == 'abs':
            primary.append(('lin', tg['y'], score_linear, lambda r: r))
        primary.append(('log', tg['ylog'], score_log, lambda r, t=t: self._to_log(t, r)))
        out = {}
        for space, obs, scorer, tf in primary:
            res = {m: tf(r) for m, r in results_by_model.items()}
            if scheme == 'random':
                per = {m: per_split_scores(r, obs, scorer) for m, r in res.items()}
                out[space] = {'per': per}
            else:
                out[space] = {'pooled': {m: pooled_predictions(r, self.n)
                                         for m, r in res.items()}, 'obs': obs,
                              'scorer': scorer}
        return out

    # -- public ------------------------------------------------------------
    def evaluate(self, fitter_factory, label, targets=None, n_boot=2000):
        """Score a model on every target, both schemes, against the baselines.

        Parameters
        ----------
        fitter_factory : callable
            ``fitter_factory(y, constraint) -> fitter``; ``constraint`` is
            ``'pigment'`` for concentrations (clip at 0) and ``None`` for log
            ratios. A fitter maps ``(train_idx, test_idx)`` to predictions in
            target units.
        label : str
            Model name used in the output.
        targets : list of str, optional
            Subset of target keys; default all.
        n_boot : int, optional
            Bootstrap resamples for the LOCO paired comparison.

        Returns
        -------
        DataFrame
            One row per target with, per scheme: the model's primary R²
            (linear for ``abs``, log for ``ratio``) and log R²; the best null
            and its R²; ``dR2`` (model − best null); the random-split fraction
            of paired wins or the LOCO bootstrap 95% CI; ``beats_null_random``
            (wins in ≥ 90% of paired splits); ``beats_null_loco`` (CI lower
            bound > 0); and, for ratios, the LOCO RMS of model and ``const``.
        """
        rows = []
        for t in (targets or list(self.targets)):
            tg = self.targets[t]
            constraint = 'pigment' if tg['kind'] == 'abs' else None
            fitter = fitter_factory(tg['y'], constraint)
            row = {'model': label, 'target': t, 'kind': tg['kind'], 'pigment': tg['pigment']}
            for scheme, sp in self.splits.items():
                res = dict(self._baseline[t][scheme])
                res[label] = cross_validate(fitter, sp)
                sc = self._score_scheme(t, scheme, res)
                space = 'lin' if tg['kind'] == 'abs' else 'log'
                if scheme == 'random':
                    per = sc[space]['per']
                    r2 = {m: float(np.nanmean(v['R2'])) for m, v in per.items()}
                    best = max(self.NULLS, key=r2.get)
                    pr = paired(per[label], per[best])
                    row.update({
                        'R2_random': r2[label],
                        'R2_random_sd': float(np.nanstd(per[label]['R2'], ddof=1)),
                        'logR2_random': float(np.nanmean(sc['log']['per'][label]['R2'])),
                        'null_random': best, 'R2_null_random': r2[best],
                        'dR2_random': pr['R2']['mean_diff'],
                        'frac_wins_random': pr['R2']['frac_model_better'],
                        'beats_null_random': pr['R2']['frac_model_better'] >= 0.9})
                else:
                    pooled, obs, scorer = (sc[space]['pooled'], sc[space]['obs'],
                                           sc[space]['scorer'])
                    r2 = {m: scorer(obs, p)['R2'] for m, p in pooled.items()}
                    best = max(self.NULLS, key=r2.get)
                    bs = bootstrap_paired(obs, pooled[label], pooled[best], scorer,
                                          n_boot=n_boot, seed=self.seed)
                    row.update({
                        'R2_loco': r2[label],
                        'logR2_loco': score_log(sc['log']['obs'],
                                                sc['log']['pooled'][label])['R2'],
                        'null_loco': best, 'R2_null_loco': r2[best],
                        'dR2_loco': bs['dR2'], 'dR2_loco_lo': bs['dR2_ci'][0],
                        'dR2_loco_hi': bs['dR2_ci'][1],
                        'beats_null_loco': bs['dR2_ci'][0] > 0})
                    if tg['kind'] == 'ratio':
                        row['RMS_loco'] = score_log(obs, pooled[label])['RMS']
                        row['RMS_const_loco'] = score_log(obs, pooled['const'])['RMS']
            rows.append(row)
        return pd.DataFrame(rows).set_index('target')

    def pcr_factory(self, X, seed=None, **kw):
        """``fitter_factory`` for Kramer's PCR on predictor matrix ``X``."""
        s = self.seed if seed is None else seed
        return lambda y, constraint: PCRFitter(X, y, constraint, seed=s, **kw)
