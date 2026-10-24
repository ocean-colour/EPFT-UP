"""Bayesian shared low-rank weighting with predictive uncertainties (Execution #8).

:class:`BayesianSharedW`:

1. Learns the shared latent weight spectra exactly as Execution #7's best
   model: :class:`epft_up.sdp.weighting.PenalizedRRR` (noise + smoothness
   penalty, reduced rank), with hyperparameters chosen by leave-one-campaign-out
   inside the training set. Optionally it adds a known test-noise covariance
   to the penalty (``Sigma_abs``; noise-aware training).
2. Projects every spectrum onto the r latent directions,
   T = X B_pen V_r (n × r).
3. For each target, fits an evidence-maximized Bayesian linear regression on
   T (sklearn ``BayesianRidge``: Gaussian prior on the r coefficients,
   Gaussian noise; α, λ by type-II maximum likelihood). Its predictive
   distribution is Gaussian, with variance 1/α + tᵀ S t: residual scatter
   plus coefficient uncertainty.

With ``calibration='insample'`` the predictive variance is conditional on
the learned latent directions, which were fitted to the same targets, so the
residual variance is optimistic (≈ 64%/90% coverage of nominal 68/95% under
LOCO). The default ``calibration='oof'`` fixes this by cross-fitting. Out-of-fold
predictions ŷ_oof of the shared W are made by leave-one-campaign-out inside
the training set at the selected hyperparameters. Each target's Bayesian
regression is then y = a + b ŷ_oof + ε, whose residual variance is an honest,
between-campaign error. At prediction time the full-training-set shared W
supplies ŷ.

Input-noise propagation is done by the caller, by Monte Carlo through the
GSM residual (see ``scripts/sdp/uncertainty.py``).
"""

from __future__ import annotations

import numpy as np
from sklearn.linear_model import BayesianRidge

from epft_up.sdp import weighting as W


class BayesianSharedW:
    """Shared low-rank W + per-target Bayesian regression on its latent scores.

    Parameters
    ----------
    lam_noise, lam_smooth, ranks : sequences
        Hyperparameter grids for the inner leave-campaign-out selection.
    Sigma_n : ndarray, optional
        Noise-shape matrix of the tuned noise penalty (as in #7).
    Sigma_abs : ndarray, optional
        Known test-noise covariance added untuned (noise-aware training).
    """

    def __init__(self, lam_noise, lam_smooth, ranks, Sigma_n=None, Sigma_abs=None,
                 lam_aux=1.0, calibration='oof'):
        if calibration not in ('oof', 'insample'):
            raise ValueError("calibration must be 'oof' or 'insample'")
        self.grid = (tuple(lam_noise), tuple(lam_smooth), tuple(ranks))
        self.Sigma_n, self.Sigma_abs, self.lam_aux = Sigma_n, Sigma_abs, lam_aux
        self.calibration = calibration

    def fit(self, Xs, Xa, Y, groups):
        """Fit on training spectra.

        Parameters
        ----------
        Xs : ndarray, shape (n, p_s)
            Spectral predictors (e.g. δRrs).
        Xa : ndarray, shape (n, p_a) or None
            Auxiliary predictors (GSM parameters).
        Y : ndarray, shape (n, q)
            Targets in their own units (log ratios or concentrations).
        groups : ndarray, shape (n,)
            Campaign labels for the inner leave-one-campaign-out selection.

        Returns
        -------
        self
        """
        self.model = W.PenalizedRRR(Xs, Xa, self.Sigma_n, *self.grid, lam_aux=self.lam_aux,
                                    Sigma_abs=self.Sigma_abs)
        idx = np.arange(Y.shape[0])
        g = np.asarray(groups)
        folds = [idx[g == c] for c in dict.fromkeys(g)]
        self.hp = self.model.select(idx, Y, folds)
        B, self.st, self.Ym, self.Ysd = self.model.fit_full(idx, Y, self.hp[0], self.hp[1])
        X, _ = self.model._design(idx, self.st)
        _, _, Vt = np.linalg.svd(X @ B, full_matrices=False)
        self.L = B @ Vt[:self.hp[2]].T                 # (p, r) latent directions
        self._latent_coef = Vt[:self.hp[2]]            # (r, q): ŷ_z = (X L) Vrᵀ
        Yz = (Y - self.Ym) / self.Ysd
        if self.calibration == 'insample':
            T = X @ self.L
            self.regs = [BayesianRidge(fit_intercept=True).fit(T, Yz[:, j])
                         for j in range(Y.shape[1])]
            return self
        # cross-fitted (out-of-fold) predictions of the shared W, standardized units
        oof = np.full_like(Yz, np.nan)
        for f in folds:
            trf = np.setdiff1d(idx, f)
            Bf, stf, Ymf, Ysdf = self.model.fit_full(trf, Y, self.hp[0], self.hp[1])
            Xf, _ = self.model._design(trf, stf)
            Brf, _ = W.PenalizedRRR.reduce_rank(Xf, Bf, self.hp[2])
            pred = self.model.predict(f, Brf, stf, Ymf, Ysdf)
            oof[f] = (pred - self.Ym) / self.Ysd
        self.regs = [BayesianRidge(fit_intercept=True).fit(oof[:, [j]], Yz[:, j])
                     for j in range(Y.shape[1])]
        self.oof_, self.Yz_ = oof, Yz
        return self

    def oof_residuals(self):
        """Cross-fitted residuals y - E[y | ŷ_oof] on the training set, in target units.

        These are leave-one-campaign-out errors of the calibrated model (the
        hyperparameters were selected on all training campaigns). Grouping them
        by campaign separates a between-campaign component from the
        within-campaign scatter (e.g. :func:`campaign_variance`).
        Only available for ``calibration='oof'``.

        Returns
        -------
        ndarray, shape (n, q)
        """
        if self.calibration != 'oof':
            raise ValueError("oof_residuals needs calibration='oof'")
        R = np.empty_like(self.Yz_)
        for j, r in enumerate(self.regs):
            R[:, j] = (self.Yz_[:, j] - r.predict(self.oof_[:, [j]])) * self.Ysd[j]
        return R

    def _latent(self, Xs, Xa):
        X = np.asarray(Xs, float) - self.st['mu_s']
        if Xa is not None:
            X = np.hstack([X, (np.atleast_2d(Xa) - self.st['mu_a']) / self.st['sd_a']])
        T = X @ self.L
        if self.calibration == 'oof':
            # full-fit shared-W prediction (standardized), the regressor of each target
            return T @ self._latent_coef
        return T

    def predict(self, Xs, Xa):
        """Predictive mean and SD in target units, each of shape (n, q)."""
        T = self._latent(Xs, Xa)
        mu = np.empty((T.shape[0], len(self.regs)))
        sd = np.empty_like(mu)
        for j, r in enumerate(self.regs):
            Tj = T[:, [j]] if self.calibration == 'oof' else T
            m, s = r.predict(Tj, return_std=True)
            mu[:, j], sd[:, j] = m * self.Ysd[j] + self.Ym[j], s * self.Ysd[j]
        return mu, sd

    def effective_weights(self):
        """Spectral weights on raw Xs per target, (p_s, q)."""
        p_s = len(self.st['mu_s'])
        if self.calibration == 'oof':
            b = np.array([r.coef_[0] for r in self.regs])
            return (self.L[:p_s] @ self._latent_coef) * b * self.Ysd
        C = np.column_stack([r.coef_ for r in self.regs])          # (r, q)
        return (self.L[:p_s] @ C) * self.Ysd


def campaign_variance(resid, groups):
    """Method-of-moments split of residuals into between- and within-campaign variance.

    With e_ij = b_j + ε_ij, b_j ~ N(0, τ²) and ε_ij ~ N(0, σ_w²):
    τ² = max(0, Var_j(ē_j) - mean_j(s_j² / n_j)) and σ_w² is the pooled
    within-campaign variance.

    Parameters
    ----------
    resid : ndarray, shape (n,) or (n, q)
    groups : ndarray, shape (n,)

    Returns
    -------
    tau2, sigw2 : ndarray, shape (q,)
    """
    R = np.atleast_2d(np.asarray(resid, float).T).T
    g = np.asarray(groups)
    labs = list(dict.fromkeys(g))
    means = np.array([R[g == c].mean(axis=0) for c in labs])
    n = np.array([(g == c).sum() for c in labs])
    s2 = np.array([R[g == c].var(axis=0, ddof=1) if k > 1 else np.zeros(R.shape[1])
                   for c, k in zip(labs, n)])
    sigw2 = np.sum((n[:, None] - 1) * s2, axis=0) / max(np.sum(n - 1), 1)
    tau2 = np.maximum(0.0, means.var(axis=0, ddof=1) - np.mean(s2 / n[:, None], axis=0))
    return tau2, sigw2


def gaussian_scores(y, mu, sd):
    """Coverage, PIT and CRPS of Gaussian predictive distributions.

    Returns
    -------
    dict
        ``cov68``, ``cov95`` (fraction inside ±1σ / ±1.96σ), ``pit`` (array),
        ``crps`` (mean; Gneiting & Raftery 2007 closed form), ``z_sd`` (SD of
        standardized errors; 1 if calibrated).
    """
    from scipy.stats import norm
    y, mu, sd = (np.asarray(a, float).ravel() for a in (y, mu, sd))
    z = (y - mu) / sd
    crps = sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
    return {'cov68': float(np.mean(np.abs(z) <= 1.0)),
            'cov95': float(np.mean(np.abs(z) <= 1.96)),
            'pit': norm.cdf(z), 'crps': float(np.mean(crps)),
            'z_sd': float(np.std(z))}
