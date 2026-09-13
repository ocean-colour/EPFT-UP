"""Retraining MOANA from a matchup table (design doc §6).

Reproduces Lange et al. (2020) §2.4 as literally as the documents allow:

- :func:`fit_pca` — ``prcomp``-equivalent SVD PCA with the **centring
  convention configurable** (design §4.4): R's default column-centres, but
  the operational retrieval projects raw standardised spectra onto the
  loadings with no mean subtracted, so which convention Lange used is
  unknown — :func:`compare_loadings` settles it empirically.
- :func:`truncate_pcs` — Lange's first cut (PC sd ≥ 0.1 % of PC1's sd).
- :func:`backward_stepwise_aic` — Lange's selection: repeatedly drop the
  highest-p predictor, keeping the removal only if AIC decreases.
- :func:`train_moana` — the whole §6 recipe, with a per-taxon ``use_sst``
  switch because the paper's Eq. 7 gives picoeukaryotes an SST term that the
  ATBD and the operational LUT do not have (design §4.3).
- :func:`compare_loadings` — |cosine| similarity of our basis against the
  vendored NASA LUT, sign/order-invariant.

Plain numpy/scipy OLS; no statsmodels dependency. Stepwise selection at
n ≈ 70 is statistically shaky (report §9.5) — implemented anyway because the
goal is *reproduction*, not improvement.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from epft_up.moana.algorithm import standardize
from epft_up.moana.io import load_luts


def fit_pca(rrs_std, center=False):
    """``prcomp``-equivalent PCA of standardised spectra (design §6.1).

    Parameters
    ----------
    rrs_std : (n, w) array — spectra from
        :func:`epft_up.moana.algorithm.standardize`; rows with NaN dropped.
    center : bool, optional
        Column-centre before the SVD (R ``prcomp`` default). False by
        default because the operational retrieval implies no training mean
        (design §4.4); fit both and let :func:`compare_loadings` decide.

    Returns
    -------
    dict with
        ``V`` : (w, p) float64 — loadings (columns = components),
        ``scores`` : (n, p) float64 — projections of the training spectra,
        ``sdev`` : (p,) float64 — component standard deviations
            (singular values / sqrt(n−1), matching prcomp's ``sdev``),
        ``mean`` : (w,) float64 — column means subtracted (zeros if
            ``center=False``),
        ``center`` : bool — echo of the convention used.
    """
    X = np.asarray(rrs_std, dtype=np.float64)
    X = X[np.isfinite(X).all(axis=1)]
    n = X.shape[0]
    mean = X.mean(axis=0) if center else np.zeros(X.shape[1])
    Xc = X - mean
    # Economy SVD: Xc = U S Vt; loadings are Vt.T, sdev = S/sqrt(n-1).
    _, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return {'V': Vt.T, 'scores': Xc @ Vt.T, 'sdev': S / np.sqrt(n - 1),
            'mean': mean, 'center': center}


def truncate_pcs(sdev, frac=1e-3):
    """Lange's first component cut: keep PCs with sd ≥ ``frac``·sd(PC1).

    (Their text: "PCs with a standard deviation lower than 0.1% of the
    standard deviation of the first PC were discarded", leaving 20 of the
    hyperspectral components.)

    Parameters
    ----------
    sdev : (p,) array — from :func:`fit_pca`.
    frac : float, optional — 1e-3 per the paper.

    Returns
    -------
    int — number of components retained.
    """
    sdev = np.asarray(sdev)
    return int((sdev >= frac * sdev[0]).sum())


def _ols(X, y):
    """Ordinary least squares with the statistics stepwise selection needs.

    Parameters
    ----------
    X : (n, k) array — design matrix (no intercept column; one is added).
    y : (n,) array — response.

    Returns
    -------
    dict — ``beta`` ((k+1,), intercept first), ``pvalues`` ((k,), slope
    p-values only), ``aic`` (Gaussian log-likelihood AIC, the quantity R's
    ``AIC`` reports up to a constant), ``resid`` ((n,)), ``r2_adj``.
    """
    n, k = X.shape
    A = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    rss = float(resid @ resid)
    dof = n - (k + 1)
    # Slope p-values via the usual t statistics.
    sigma2 = rss / dof if dof > 0 else np.nan
    cov = sigma2 * np.linalg.pinv(A.T @ A)
    se = np.sqrt(np.diag(cov))
    with np.errstate(invalid='ignore', divide='ignore'):
        tstat = beta / se
    pvals = 2 * stats.t.sf(np.abs(tstat[1:]), dof)
    # Gaussian AIC as R computes it: n*log(RSS/n) + 2*(k+2) + const; the
    # constant cancels in comparisons, k+2 counts intercept + sigma.
    aic = n * np.log(rss / n) + 2 * (k + 2)
    tss = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - rss / tss if tss > 0 else np.nan
    r2_adj = 1 - (1 - r2) * (n - 1) / dof if dof > 0 else np.nan
    return {'beta': beta, 'pvalues': pvals, 'aic': aic, 'resid': resid,
            'r2_adj': r2_adj}


def backward_stepwise_aic(X, y, names):
    """Lange's predictor selection (design §6.2).

    Start from all candidates; repeatedly refit, find the highest-p
    predictor, and remove it **only if the model without it has lower AIC**;
    stop the first time a removal fails to lower AIC.

    Parameters
    ----------
    X : (n, k) array — candidate predictors, one column each.
    y : (n,) array — response (already log10-transformed where applicable).
    names : list of str — column labels (e.g. ``['logSST', 'U1', ...]``).

    Returns
    -------
    dict — ``names`` (retained labels), ``fit`` (the final :func:`_ols`
    result), ``dropped`` (labels in removal order).
    """
    keep = list(range(X.shape[1]))
    dropped = []
    fit = _ols(X[:, keep], y)
    while len(keep) > 1:
        if np.isfinite(fit['pvalues']).any():
            worst = int(np.nanargmax(fit['pvalues']))
            trial = keep[:worst] + keep[worst + 1:]
            trial_fit = _ols(X[:, trial], y)
            if trial_fit['aic'] >= fit['aic']:
                break
        else:
            # Saturated model (dof <= 0): p-values and AIC are undefined, so
            # shed the last (lowest-variance PC) predictor unconditionally
            # until the fit is estimable, then resume the AIC rule.
            worst = len(keep) - 1
            trial = keep[:worst]
            trial_fit = _ols(X[:, trial], y)
        dropped.append(names[keep[worst]])
        keep, fit = trial, trial_fit
    return {'names': [names[j] for j in keep], 'fit': fit,
            'dropped': dropped}


def train_moana(rrs, counts, sst=None, center=False, sd_frac=1e-3,
                use_sst=None):
    """The full §6 retraining recipe on a matchup table.

    Parameters
    ----------
    rrs : (n, 124) array — matched Rrs spectra on the MOANA grid (raw, not
        yet standardised).
    counts : dict — ``{'pro': (n,), 'syn': (n,), 'peuk': (n,)}`` cell
        abundances [cells mL⁻¹]; NaNs allowed (per-taxon rows drop).
    sst : (n,) array, optional — SST [°C] for the taxa whose model uses it.
    center : bool, optional — PCA centring convention (design §4.4).
    sd_frac : float, optional — component cut (Lange: 1e-3 → 20 PCs).
    use_sst : dict, optional — per-taxon SST switch; default
        ``{'pro': True, 'syn': False, 'peuk': False}`` (the operational
        model). ``{'peuk': True}`` tests the paper's Eq. 7 variant
        (design §4.3).

    Returns
    -------
    dict with
        ``pca`` : dict — from :func:`fit_pca`,
        ``n_pc`` : int — components surviving the sd cut,
        ``models`` : dict per taxon — ``names``, ``coef`` (aligned with
            names, intercept first), ``r2_adj``, ``dropped``, ``response``
            ('linear' for pro, 'log10' otherwise), ``n``,
        ``rrs_std`` : (n, 124) — the standardised training spectra.
    """
    use = {'pro': True, 'syn': False, 'peuk': False, **(use_sst or {})}
    rrs_std = standardize(rrs)
    pca = fit_pca(rrs_std, center=center)
    n_pc = truncate_pcs(pca['sdev'], sd_frac)
    # Ours, not Lange's: cap the candidate pool so the initial regression
    # keeps >= 10 residual dof (predictors = n_pc + SST + intercept). In the
    # Lange regime (n ~ 73, 20 PCs) this never binds; it only prevents a
    # saturated fit when n is small or the sd cut is loose.
    n = np.isfinite(rrs_std).all(axis=1).sum()
    n_pc = min(n_pc, max(1, int(n) - 12))
    # Scores of *all* rows (fit_pca dropped NaNs internally; recompute here
    # so row alignment with counts/sst is exact).
    scores = (rrs_std - pca['mean']) @ pca['V'][:, :n_pc]

    models = {}
    for taxon in ('pro', 'syn', 'peuk'):
        y_raw = np.asarray(counts[taxon], dtype=np.float64)
        # Pro is fitted on linear cells/mL (Lange: log-transform "reduced
        # the performance"); the other two on log10 (design §3 step 5).
        response = 'linear' if taxon == 'pro' else 'log10'
        y = y_raw if response == 'linear' else np.log10(y_raw)
        cols, names = [], []
        if use[taxon]:
            if sst is None:
                raise ValueError(f"use_sst[{taxon!r}] is True but sst=None")
            cols.append(np.log10(np.asarray(sst, dtype=np.float64)))
            names.append('logSST')
        for j in range(n_pc):
            cols.append(scores[:, j])
            names.append(f'U{j + 1}')
        X = np.column_stack(cols)
        ok = np.isfinite(X).all(axis=1) & np.isfinite(y)
        n_ok = int(ok.sum())
        if n_ok < 15:
            raise ValueError(
                f"train_moana: only {n_ok} usable rows for {taxon!r} "
                "(need >= 15) — check counts/SST NaNs and the matchup window.")
        # Per-taxon dof guard, mirroring the global n_pc cap above: NaNs in
        # this taxon's counts/SST may leave fewer rows than the matrix had.
        max_cols = max(1, n_ok - 12)
        if X.shape[1] > max_cols:
            X, names = X[:, :max_cols], names[:max_cols]
        sel = backward_stepwise_aic(X[ok], y[ok], names)
        models[taxon] = {'names': sel['names'],
                         'coef': sel['fit']['beta'],
                         'r2_adj': sel['fit']['r2_adj'],
                         'dropped': sel['dropped'],
                         'response': response,
                         'n': int(ok.sum())}
    return {'pca': pca, 'n_pc': n_pc, 'models': models, 'rrs_std': rrs_std}


def compare_loadings(V_ours, V_nasa=None, n_compare=20):
    """Sign/order-invariant comparison of a retrained basis to NASA's (§6.3).

    Parameters
    ----------
    V_ours : (124, p) array — retrained loadings.
    V_nasa : (124, 45) array, optional — defaults to the vendored LUT's.
    n_compare : int, optional — leading NASA components to match.

    Returns
    -------
    dict with
        ``similarity`` : (n_compare, p) float64 — |cosine| matrix,
        ``best_match`` : (n_compare,) int — our column best matching each
            NASA component,
        ``best_cos`` : (n_compare,) float64 — the |cosine| of that match
            (success criterion: > 0.99 on the leading components under one
            of the two centring conventions).
    """
    if V_nasa is None:
        V_nasa = load_luts()['V']
    A = np.asarray(V_nasa)[:, :n_compare]
    B = np.asarray(V_ours)
    # Columns are unit vectors up to float error; normalise anyway.
    A = A / np.linalg.norm(A, axis=0)
    B = B / np.linalg.norm(B, axis=0)
    sim = np.abs(A.T @ B)
    best = sim.argmax(axis=1)
    return {'similarity': sim, 'best_match': best,
            'best_cos': sim[np.arange(sim.shape[0]), best]}
