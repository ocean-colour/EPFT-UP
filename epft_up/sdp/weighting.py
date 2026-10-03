"""Learned wavelength weightings (Execution #7).

Estimators, all linear in the predictors and usable as
:class:`epft_up.sdp.validate.Benchmark` fitters:

* :class:`RidgeGCVFitter`: ridge on z-scored predictors, λ by generalized
  cross-validation (closed form via the SVD). A per-pigment control.
* :class:`PLSFitter`: PLS1 on z-scored predictors, number of components by
  inner k-fold CV. A per-pigment control.
* :class:`PenalizedRRR`: the **shared low-rank weighting** W. A
  multi-output regression of all targets of one kind on raw-unit spectral
  predictors (⊕ standardized auxiliaries), with the penalty

      λ_n wᵀ Σ_n w   (noise: the prediction variance Rrs noise would cause)
    + λ_s ‖D₂ w‖²    (smoothness of each weight spectrum)
    + ε ‖w‖²         (small ridge for conditioning)

  on the spectral block and a ridge on the auxiliaries. It is followed by
  reduced-rank projection (Izenman 1975; Mukherjee & Zhu 2011): with
  B_pen = (XᵀX + P)⁻¹ XᵀY and the right singular vectors V_r of the fitted
  values X B_pen, B = B_pen V_r V_rᵀ. The r columns of B_pen V_r are the
  shared latent weight spectra. Hyperparameters (λ_n, λ_s, r) are chosen by
  inner CV on the mean squared error of the standardized targets: random
  k-fold, or **leave-one-campaign-out within the training set**
  (``inner='group'``), which selects for settings that transfer between
  campaigns.

The noise term is the natural form of "noise whitening" for a linear
predictor: if the predictors carry noise of covariance Σ_n, the induced
prediction variance is wᵀ Σ_n w (report §5.2), so penalizing it
down-weights bands in proportion to the noise they inject.
"""

from __future__ import annotations

import numpy as np

from epft_up.sdp import theory


def _clip(pred, constraint):
    return np.maximum(pred, 0.0) if constraint == 'pigment' else pred


# --------------------------------------------------------------------------
# Per-pigment controls
# --------------------------------------------------------------------------
class RidgeGCVFitter:
    """Ridge on z-scored X, λ chosen by GCV on the training set.

    Parameters
    ----------
    X : ndarray, shape (n, p)
    y : ndarray, shape (n,)
    constraint : {'pigment', None}
    lams : ndarray, optional
        Grid of λ relative to the mean squared singular value of the training Z.
    """

    def __init__(self, X, y, constraint='pigment', lams=None):
        self.X, self.y, self.constraint = np.asarray(X, float), np.asarray(y, float), constraint
        self.lams = np.logspace(-4, 2, 31) if lams is None else np.asarray(lams)
        self.chosen = []

    def __call__(self, tr, te):
        Xt, yt = self.X[tr], self.y[tr]
        mu, sd = Xt.mean(0), Xt.std(0, ddof=1)
        sd[sd == 0] = 1.0
        Z = (Xt - mu) / sd
        yc = yt - yt.mean()
        U, s, Vt = np.linalg.svd(Z, full_matrices=False)
        uy = U.T @ yc
        n = len(yt)
        best, best_g = None, np.inf
        scale = np.mean(s**2)
        for lr in self.lams:
            lam = lr * scale
            f = s**2 / (s**2 + lam)
            resid = yc - U @ (f * uy)
            g = n * np.dot(resid, resid) / (n - f.sum())**2
            if g < best_g:
                best, best_g = lam, g
        self.chosen.append(best / scale)
        beta = Vt.T @ ((s / (s**2 + best)) * uy)
        pred = ((self.X[te] - mu) / sd) @ beta + yt.mean()
        return _clip(pred, self.constraint)


class PLSFitter:
    """PLS1 on z-scored X; number of components by inner k-fold CV."""

    def __init__(self, X, y, constraint='pigment', max_comp=15, k=5, seed=1):
        self.X, self.y, self.constraint = np.asarray(X, float), np.asarray(y, float), constraint
        self.max_comp, self.k = max_comp, k
        self.rng = np.random.default_rng(seed)
        self.chosen = []

    @staticmethod
    def _fit(Z, y, nc):
        from sklearn.cross_decomposition import PLSRegression
        m = PLSRegression(n_components=nc, scale=False).fit(Z, y)
        return m

    def __call__(self, tr, te):
        Xt, yt = self.X[tr], self.y[tr]
        mu, sd = Xt.mean(0), Xt.std(0, ddof=1)
        sd[sd == 0] = 1.0
        Z = (Xt - mu) / sd
        folds = np.array_split(self.rng.permutation(len(tr)), self.k)
        err = np.zeros(self.max_comp)
        for f in folds:
            m = np.ones(len(tr), bool)
            m[f] = False
            for nc in range(1, self.max_comp + 1):
                pr = self._fit(Z[m], yt[m], nc).predict(Z[f]).ravel()
                err[nc - 1] += np.sum((_clip(pr, self.constraint) - yt[f])**2)
        nc = int(np.argmin(err)) + 1
        self.chosen.append(nc)
        pred = self._fit(Z, yt, nc).predict((self.X[te] - mu) / sd).ravel()
        return _clip(pred, self.constraint)


# --------------------------------------------------------------------------
# Shared low-rank weighting
# --------------------------------------------------------------------------
class PenalizedRRR:
    """Penalized reduced-rank regression (see module docstring).

    Parameters
    ----------
    Xs : ndarray, shape (n, p_s)
        Spectral predictors in physical units (e.g. δRrs [sr⁻¹]).
    Xa : ndarray, shape (n, p_a) or None
        Auxiliary predictors (e.g. GSM parameters); standardized internally.
    Sigma_n : ndarray, shape (p_s, p_s)
        Noise covariance of the spectral predictors.
    lam_noise, lam_smooth : sequence of float
        Grids, relative to tr(XsᵀXs) after each penalty matrix is normalized
        to the same trace.
    ranks : sequence of int
    lam_aux : float
        Ridge on the standardized auxiliaries (relative to n).
    eps : float
        Relative conditioning ridge on the spectral block.
    Sigma_abs : ndarray, shape (p_s, p_s), optional
        Known test-time noise covariance of the spectral predictors, in their
        physical units. If given, ``n_train · wᵀ Σ_abs w`` is added to the
        objective with no tuning: for standardized targets this is exactly
        the expected extra squared error that noise of covariance Σ_abs adds
        to a linear prediction (the errors-in-variables / noise-augmented
        ridge; Execution #8).
    """

    def __init__(self, Xs, Xa=None, Sigma_n=None, lam_noise=(0.0, 0.1, 1.0),
                 lam_smooth=(1e-3, 1e-2, 1e-1, 1.0), ranks=range(1, 11), lam_aux=1.0,
                 eps=1e-6, Sigma_abs=None):
        self.Xs = np.asarray(Xs, float)
        self.Xa = None if Xa is None else np.atleast_2d(np.asarray(Xa, float))
        p = self.Xs.shape[1]
        D = theory.difference_matrix(p, 2)
        self.R = D.T @ D
        self.R /= np.trace(self.R)
        self.N = np.eye(p) if Sigma_n is None else np.asarray(Sigma_n, float)
        self.N = self.N / np.trace(self.N)
        self.I = np.eye(p) / p
        self.lam_noise, self.lam_smooth = tuple(lam_noise), tuple(lam_smooth)
        self.ranks, self.lam_aux, self.eps = tuple(ranks), lam_aux, eps
        self.Sigma_abs = None if Sigma_abs is None else np.asarray(Sigma_abs, float)

    # ---- core fit -------------------------------------------------------
    def _design(self, idx, stats=None):
        Xs = self.Xs[idx]
        if stats is None:
            stats = {'mu_s': Xs.mean(0)}
            if self.Xa is not None:
                Xa = self.Xa[idx]
                stats['mu_a'] = Xa.mean(0)
                sd = Xa.std(0, ddof=1)
                stats['sd_a'] = np.where(sd > 0, sd, 1.0)
        X = Xs - stats['mu_s']
        if self.Xa is not None:
            X = np.hstack([X, (self.Xa[idx] - stats['mu_a']) / stats['sd_a']])
        return X, stats

    def _penalty(self, trXs, n, lam_n, lam_s):
        P = trXs * (lam_n * self.N + lam_s * self.R + self.eps * self.I)
        if self.Sigma_abs is not None:
            P = P + n * self.Sigma_abs
        if self.Xa is None:
            return P
        pa = self.Xa.shape[1]
        out = np.zeros((P.shape[0] + pa, P.shape[0] + pa))
        out[:P.shape[0], :P.shape[0]] = P
        out[P.shape[0]:, P.shape[0]:] = self.lam_aux * n * np.eye(pa)
        return out

    def fit_full(self, tr, Y, lam_n, lam_s):
        """Penalized multi-output coefficients B_pen (centred), with stats."""
        X, st = self._design(tr)
        Ym = Y[tr].mean(0)
        Ysd = Y[tr].std(0, ddof=1)
        Ysd = np.where(Ysd > 0, Ysd, 1.0)
        Yz = (Y[tr] - Ym) / Ysd
        p_s = self.Xs.shape[1]
        trXs = np.sum(X[:, :p_s]**2)
        P = self._penalty(trXs, len(tr), lam_n, lam_s)
        B = np.linalg.solve(X.T @ X + P, X.T @ Yz)
        return B, st, Ym, Ysd

    @staticmethod
    def reduce_rank(X, B, r):
        """B V_r V_rᵀ with V_r from the SVD of the fitted values X B."""
        _, _, Vt = np.linalg.svd(X @ B, full_matrices=False)
        Vr = Vt[:r].T
        return B @ Vr @ Vr.T, Vr

    def predict(self, idx, B, st, Ym, Ysd):
        X, _ = self._design(idx, st)
        return X @ B * Ysd + Ym

    # ---- hyperparameter selection --------------------------------------
    def select(self, tr, Y, inner_folds):
        """Grid search of (λ_n, λ_s, r) on standardized-target MSE; returns the best."""
        best, best_err = None, np.inf
        for lam_n in self.lam_noise:
            for lam_s in self.lam_smooth:
                err = np.zeros(len(self.ranks))
                for f in inner_folds:
                    trf = np.setdiff1d(tr, f)
                    B, st, Ym, Ysd = self.fit_full(trf, Y, lam_n, lam_s)
                    Xtr, _ = self._design(trf, st)
                    Xf, _ = self._design(f, st)
                    _, _, Vt = np.linalg.svd(Xtr @ B, full_matrices=False)
                    Yf = (Y[f] - Ym) / Ysd
                    for j, r in enumerate(self.ranks):
                        rr = min(r, Vt.shape[0])
                        Vr = Vt[:rr].T
                        err[j] += np.sum((Xf @ B @ Vr @ Vr.T - Yf)**2)
                j = int(np.argmin(err))
                if err[j] < best_err:
                    best_err, best = err[j], (lam_n, lam_s, self.ranks[j])
        return best


class SharedRRRFitter:
    """Benchmark adaptor: one joint PenalizedRRR fit per training set, all targets.

    ``factory = SharedRRRFitter(model, Y, campaigns, ...)``; then
    ``factory(y, constraint)`` returns a fitter for the column of ``Y`` equal
    to ``y``. Fits are cached per training set, so each split is solved once
    for all targets.

    Parameters
    ----------
    model : PenalizedRRR
    Y : ndarray, shape (n, q)
        All targets of one kind (in Benchmark units).
    campaigns : ndarray, shape (n,)
        Group labels for ``inner='group'``.
    inner : {'random', 'group'}
    k : int
        Inner folds for ``inner='random'``.
    fixed : tuple, optional
        (λ_n, λ_s, r) to skip selection.
    """

    def __init__(self, model, Y, campaigns, inner='random', k=5, seed=1, fixed=None):
        self.model, self.Y = model, np.asarray(Y, float)
        self.campaigns = np.asarray(campaigns)
        self.inner, self.k, self.fixed = inner, k, fixed
        self.rng = np.random.default_rng(seed)
        self.cache, self.chosen = {}, []

    def _folds(self, tr):
        if self.inner == 'group':
            g = self.campaigns[tr]
            return [tr[g == c] for c in dict.fromkeys(g)]
        return [tr[f] for f in np.array_split(self.rng.permutation(len(tr)), self.k)]

    def _fit(self, tr):
        key = tr.tobytes()
        if key not in self.cache:
            hp = self.fixed or self.model.select(tr, self.Y, self._folds(tr))
            B, st, Ym, Ysd = self.model.fit_full(tr, self.Y, hp[0], hp[1])
            X, _ = self.model._design(tr, st)
            Br, _ = PenalizedRRR.reduce_rank(X, B, hp[2])
            allidx = np.arange(self.Y.shape[0])
            self.cache[key] = self.model.predict(allidx, Br, st, Ym, Ysd)
            self.chosen.append(hp)
        return self.cache[key]

    def __call__(self, y, constraint):
        col = [j for j in range(self.Y.shape[1]) if np.array_equal(self.Y[:, j], y)]
        if len(col) != 1:
            raise ValueError('target not found (exactly once) in Y')
        j = col[0]

        def fitter(tr, te):
            return _clip(self._fit(np.asarray(tr))[te, j], constraint)
        return fitter
