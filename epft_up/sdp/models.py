"""Pigment models on spectral predictors: the Kramer 2022 principal components regression.

:func:`train_pcr_kramer` re-implements ``rrsModelTrain.m`` (Catlett; edited by
Kramer, 2021), the training procedure of Kramer et al. (2022, §2.5), from its
source. For each of ``n_perm`` permutations:

1. Outer split: a random 75% training / 25% validation partition.
2. Inner k-fold CV on the 75% (k = 5). In each fold, the training spectra
   are z-scored with their own mean and SD (ddof = 1), and the held-out fold
   is z-scored with **its own** mean and SD. That is a quirk of the
   original, kept here (``valid_scaling='own'``). Uncentred PCA (the data
   are already centred) gives up to ``max_pcs`` components. For each
   l = 1..max_pcs an OLS model of the pigment on the first l scores is
   formed, back-projected to spectral coefficients, applied to the held-out
   fold, and clipped at 0. The l with the smallest validation MAE is kept
   (ties go to the smaller l).
3. The k fold coefficient vectors and intercepts are averaged, still in
   standardized units, then converted to raw-predictor units with the
   outer-training mean and SD.
4. The permutation's model is validated on the 25%: predictions are clipped
   at 0, R² and RMSE come from an OLS of observed on modelled, and MAE is
   computed after zeros in the observations are replaced by 1e-4 (as in the
   original).

Because the PC scores within a fold are orthogonal and centred, the OLS
coefficients of the first l components do not depend on l. All max_pcs
candidate models are therefore computed from one fit (exact, not an
approximation).

Kramer's settings (``Kramer_Rrs_pigments.m``): n_perm = 100, max_pcs = 30,
k = 5, metric = 'MAE', outputs constrained ≥ 0, and predictors
``diff(δRrs, 2)`` on the 1 nm grid.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class PCRResult:
    """Output of :func:`train_pcr_kramer` for one target.

    Attributes
    ----------
    coefs : ndarray, shape (n_perm, n_features)
        Spectral coefficients A_m(λ) of each permutation (raw predictor units).
    intercepts : ndarray, shape (n_perm,)
        C_m of each permutation.
    n_pcs : ndarray of int, shape (n_perm, k)
        Number of PCs selected in each inner fold.
    valid_idx : list of ndarray
        Indices of the 25% validation samples of each permutation.
    stats : dict of ndarray, each shape (n_perm,)
        ``R2``, ``RMSE`` (as ``rrsModelTrain.m``), ``MAE``, ``MAE_norm_pred``
        (MAE / mean modelled value), ``MAE_norm_obs`` (MAE / mean observed),
        ``mean_pct_error``, ``median_pct_error``, ``pct_bias``.
    settings : dict
    """
    coefs: np.ndarray
    intercepts: np.ndarray
    n_pcs: np.ndarray
    valid_idx: list
    stats: dict
    settings: dict = field(default_factory=dict)

    def summary(self):
        """Mean and SD (ddof = 1, as MATLAB ``std``) of each statistic."""
        return {k: (float(np.mean(v)), float(np.std(v, ddof=1)))
                for k, v in self.stats.items()}


def _zscore(X):
    """Column z-score with ddof = 1; returns (Z, mean, sd)."""
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    return (X - mu) / sd, mu, sd


def _clip(y, constraint):
    if constraint == 'pigment':
        return np.maximum(y, 0.0)
    if constraint == 'compositions':
        return np.clip(y, 0.0, 1.0)
    if constraint in (None, 'EOFs'):
        return y
    raise ValueError(f'unknown constraint {constraint!r}')


def _ols_r2_rmse(x, y):
    """R² and RMSE of an OLS fit y = a + b x (MATLAB ``fitlm(x, y)``)."""
    n = len(x)
    xm = x - x.mean()
    sxx = np.dot(xm, xm)
    if sxx == 0.0:
        return 0.0, float(np.sqrt(np.sum((y - y.mean())**2) / max(n - 2, 1)))
    b = np.dot(xm, y - y.mean()) / sxx
    resid = y - (y.mean() + b * xm)
    sst = np.sum((y - y.mean())**2)
    r2 = 1.0 - np.dot(resid, resid) / sst if sst > 0 else 0.0
    return float(r2), float(np.sqrt(np.dot(resid, resid) / max(n - 2, 1)))


def _kfold_indices(n, k, rng):
    """Random partition of range(n) into k folds of near-equal size."""
    return np.array_split(rng.permutation(n), k)


def _fit_fold(Xtr, ytr, Xva, yva, max_pcs, constraint, valid_scaling):
    """One inner fold: choose the number of PCs by validation MAE.

    Returns the standardized-space coefficients, intercept and chosen l.
    """
    Ztr, mu, sd = _zscore(Xtr)
    if valid_scaling == 'own':
        Zva = _zscore(Xva)[0]
    elif valid_scaling == 'train':
        Zva = (Xva - mu) / sd
    else:
        raise ValueError(f'unknown valid_scaling {valid_scaling!r}')
    # Uncentred PCA of the (already centred) standardized spectra.
    _, s, Vt = np.linalg.svd(Ztr, full_matrices=False)
    L = min(max_pcs, int(np.sum(s > s[0] * 1e-12)))
    V = Vt[:L].T                                  # loadings (n_feat, L)
    scores = Ztr @ V                              # = U S, orthogonal, centred
    b = (scores.T @ ytr) / np.sum(scores**2, axis=0)
    alpha = float(ytr.mean())
    # Validation predictions of all L nested models at once.
    contrib = (Zva @ V) * b                       # (n_va, L)
    preds = _clip(np.cumsum(contrib, axis=1) + alpha, constraint)
    mae = np.mean(np.abs(preds - yva[:, None]), axis=0)
    l_best = int(np.argmin(mae)) + 1
    beta = V[:, :l_best] @ b[:l_best]
    return beta, alpha, l_best


def train_pcr_kramer(X, y, n_perm=100, max_pcs=30, k=5, constraint='pigment',
                     train_frac=0.75, seed=1, valid_scaling='own'):
    """Kramer 2022 / ``rrsModelTrain.m`` PCR training (see module docstring).

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features)
        Spectral predictors (e.g. ``diff(δRrs, 2)``).
    y : ndarray, shape (n_samples,)
        Target (pigment concentration).
    n_perm, max_pcs, k : int, optional
    constraint : {'pigment', 'compositions', None}, optional
        Output constraint (``pft_index`` of the original).
    train_frac : float, optional
    seed : int, optional
        Seed for :class:`numpy.random.Generator` (MATLAB used ``rng(1)``; the
        splits themselves cannot be reproduced across languages).
    valid_scaling : {'own', 'train'}, optional
        How inner-fold validation spectra are standardized: by their own
        statistics (the original) or by the fold's training statistics.

    Returns
    -------
    PCRResult
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, p = X.shape
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError('train_pcr_kramer: X and y must be finite')
    rng = np.random.default_rng(seed)
    n_tr = int(np.floor(n * train_frac))

    coefs = np.empty((n_perm, p))
    icpts = np.empty(n_perm)
    n_pcs = np.empty((n_perm, k), dtype=int)
    valid_idx = []
    names = ('R2', 'RMSE', 'MAE', 'MAE_norm_pred', 'MAE_norm_obs',
             'mean_pct_error', 'median_pct_error', 'pct_bias')
    stats = {s: np.empty(n_perm) for s in names}

    for i in range(n_perm):
        tr = rng.permutation(n)[:n_tr]
        va = np.setdiff1d(np.arange(n), tr)
        valid_idx.append(va)
        Xtr, ytr = X[tr], y[tr]

        betas = np.empty((k, p))
        alphas = np.empty(k)
        for j, fold in enumerate(_kfold_indices(n_tr, k, rng)):
            mask = np.ones(n_tr, dtype=bool)
            mask[fold] = False
            betas[j], alphas[j], n_pcs[i, j] = _fit_fold(
                Xtr[mask], ytr[mask], Xtr[fold], ytr[fold], max_pcs, constraint,
                valid_scaling)

        mb, ma = betas.mean(axis=0), alphas.mean()
        mu = Xtr.mean(axis=0)
        sd = Xtr.std(axis=0, ddof=1)
        coefs[i] = mb / sd
        icpts[i] = ma - np.sum(mb * mu / sd)

        pred = _clip(X[va] @ coefs[i] + icpts[i], constraint)
        obs = y[va].copy()
        stats['R2'][i], stats['RMSE'][i] = _ols_r2_rmse(pred, obs)
        obs[obs == 0] = 1e-4
        err = pred - obs
        stats['MAE'][i] = np.mean(np.abs(err))
        stats['MAE_norm_pred'][i] = stats['MAE'][i] / np.mean(pred) \
            if np.mean(pred) > 0 else np.nan
        stats['MAE_norm_obs'][i] = stats['MAE'][i] / np.mean(obs)
        pct = np.abs(err / obs) * 100
        stats['mean_pct_error'][i] = pct.mean()
        stats['median_pct_error'][i] = np.median(pct)
        stats['pct_bias'][i] = np.mean(err / obs * 100)

    return PCRResult(coefs=coefs, intercepts=icpts, n_pcs=n_pcs, valid_idx=valid_idx,
                     stats=stats,
                     settings={'n_perm': n_perm, 'max_pcs': max_pcs, 'k': k,
                               'constraint': constraint, 'train_frac': train_frac,
                               'seed': seed, 'valid_scaling': valid_scaling,
                               'selection_metric': 'MAE'})


def predict_ensemble(X, coefs, intercepts, constraint='pigment', lod=None):
    """Median over an ensemble of linear models, as in Kramer 2022 §2.5.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features)
    coefs : ndarray, shape (n_models, n_features)
    intercepts : ndarray, shape (n_models,)
    constraint : {'pigment', 'compositions', None}, optional
        Applied to the median (as the port's ``run_sdp``).
    lod : float, optional
        Values below this detection limit are set to 0 (Kramer: "modeled
        pigment values below the standard HPLC pigment detection limits ...
        were again set to zero").

    Returns
    -------
    median : ndarray, shape (n_samples,)
    all_runs : ndarray, shape (n_samples, n_models)
    """
    runs = np.asarray(X, dtype=float) @ np.asarray(coefs).T + np.asarray(intercepts)
    med = _clip(np.median(runs, axis=1), constraint)
    if lod is not None:
        med = np.where(med < lod, 0.0, med)
    return med, runs
