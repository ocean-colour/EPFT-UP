"""Linear-algebra tools for the maths/statistics of the SDP approach (Execution #6).

* :func:`difference_matrix`, :func:`boxcar_matrix`: the finite-difference and
  moving-mean operators as explicit matrices, so a whole preprocessing chain
  is one matrix ``L`` and ``predictors = L @ spectrum``.
* :func:`effective_weights`: a linear model on transformed predictors,
  ``p = Aᵀ (L x) + C``, is the linear model ``(Lᵀ A)ᵀ x + C`` on the original
  spectrum. That gives the *effective weight spectrum* on δRrs.
* :func:`implied_prior_covariance`: the Gaussian prior on the effective
  weights implied by "transform, z-score, then shrink isotropically".
* :func:`propagate_covariance`, :func:`dof_signal`: noise propagation and
  Rodgers' (2000) degrees of freedom for signal,
  d_s = tr[C_s (C_s + C_n)⁻¹].
* :func:`filter_factors`: the spectral-filter view of PCR, ridge and PLS
  (Hansen 1998; Frank & Friedman 1993).
"""

from __future__ import annotations

import numpy as np


def difference_matrix(n, order=2):
    """Matrix of MATLAB-style ``diff(·, order)``: shape (n - order, n).

    Parameters
    ----------
    n : int
    order : {1, 2}

    Returns
    -------
    ndarray
    """
    if order == 1:
        D = np.zeros((n - 1, n))
        i = np.arange(n - 1)
        D[i, i], D[i, i + 1] = -1.0, 1.0
        return D
    if order == 2:
        D = np.zeros((n - 2, n))
        i = np.arange(n - 2)
        D[i, i], D[i, i + 1], D[i, i + 2] = 1.0, -2.0, 1.0
        return D
    raise ValueError('order must be 1 or 2')


def boxcar_matrix(n, width=5):
    """Centred moving mean as an (n, n) matrix; edge rows are the identity.

    Matches :func:`epft_up.sdp.spectral.moving_mean` in the interior. The
    ``width // 2`` edge points are left unsmoothed rather than NaN.
    """
    if width % 2 == 0:
        raise ValueError('width must be odd')
    M = np.eye(n)
    h = width // 2
    for i in range(h, n - h):
        M[i] = 0.0
        M[i, i - h:i + h + 1] = 1.0 / width
    return M


def subsample_matrix(n, step):
    """Selection matrix taking every ``step``-th element: shape (ceil(n/step), n)."""
    idx = np.arange(0, n, step)
    S = np.zeros((len(idx), n))
    S[np.arange(len(idx)), idx] = 1.0
    return S


def effective_weights(coefs, L):
    """Weights on the original spectrum of a linear model on ``L @ x``.

    Parameters
    ----------
    coefs : ndarray, shape (..., m)
        Coefficients A on the transformed predictors.
    L : ndarray, shape (m, n)
        The linear transform (e.g. :func:`difference_matrix`).

    Returns
    -------
    ndarray, shape (..., n)
        ``Lᵀ A``, so that ``A · (L x) = (Lᵀ A) · x`` exactly.
    """
    return np.asarray(coefs) @ np.asarray(L)


def implied_prior_covariance(L, scale):
    """Prior covariance on effective weights implied by transform + z-score + isotropic shrinkage.

    A model ``p = βᵀ Z``, with ``Z = S⁻¹ L x`` (z-scored transformed
    predictors, ``S = diag(scale)``) and an isotropic prior ``β ~ N(0, I)``
    (the ridge / Bayesian reading of a shrinkage estimator), has effective
    weights ``w = Lᵀ S⁻¹ β`` with covariance ``Lᵀ S⁻² L``.

    Parameters
    ----------
    L : ndarray, shape (m, n)
    scale : ndarray, shape (m,)
        Per-feature standard deviation used for z-scoring.

    Returns
    -------
    ndarray, shape (n, n)
    """
    Linv = np.asarray(L) / np.asarray(scale)[:, None]
    return Linv.T @ Linv


def propagate_covariance(L, C):
    """Covariance of ``L x`` given ``Cov(x) = C``: ``L C Lᵀ``."""
    return L @ C @ L.T


def white_covariance(sigma):
    """Diagonal covariance from per-band 1-sigma ``sigma``."""
    return np.diag(np.asarray(sigma, float)**2)


def correlated_covariance(wave, sigma, ell):
    """Gaussian-correlated noise covariance: ``σ_i σ_j exp(-(λ_i - λ_j)² / 2ℓ²)``.

    A model for spectrally smooth errors (e.g. imperfect atmospheric
    correction) with correlation length ``ell`` [nm].
    """
    w = np.asarray(wave, float)
    s = np.asarray(sigma, float) * np.ones_like(w)
    return np.outer(s, s) * np.exp(-0.5 * ((w[:, None] - w[None, :]) / ell)**2)


def dof_signal(C_signal, C_noise, rcond=1e-12):
    """Rodgers (2000) degrees of freedom for signal: tr[C_s (C_s + C_n)⁻¹].

    Computed from the generalized eigenvalues λ_i of (C_s, C_n) as
    Σ λ_i / (1 + λ_i), on the subspace where C_n is positive definite (rank
    deficient noise covariances are handled with a pseudo-inverse square
    root). It is invariant under any invertible linear transform of the data,
    and can only decrease under a non-invertible one.

    Parameters
    ----------
    C_signal, C_noise : ndarray, shape (m, m)

    Returns
    -------
    d_s : float
    lam : ndarray
        Generalized eigenvalues (signal-to-noise variance ratio per independent
        direction), descending.
    """
    evals, evecs = np.linalg.eigh(C_noise)
    keep = evals > rcond * evals.max()
    W = evecs[:, keep] / np.sqrt(evals[keep])          # whitening: Wᵀ C_n W = I
    lam = np.linalg.eigvalsh(W.T @ C_signal @ W)[::-1]
    lam = np.clip(lam, 0.0, None)
    return float(np.sum(lam / (1.0 + lam))), lam


def filter_factors(Z, y, beta):
    """Spectral filter factors f_i of a linear estimator β on centred data Z.

    With the SVD Z = U S Vᵀ, the OLS solution is Σ (u_iᵀ y / s_i) v_i. Any
    estimator in the row space of Z can be written Σ f_i (u_iᵀ y / s_i) v_i.
    PCR has f_i ∈ {0, 1}; ridge has f_i = s_i² / (s_i² + λ); PLS has
    polynomial factors that can exceed 1.

    Returns
    -------
    f : ndarray, shape (rank,)
    s : ndarray
        Singular values.
    """
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    keep = s > s[0] * 1e-10
    U, s, Vt = U[:, keep], s[keep], Vt[keep]
    uy = U.T @ (y - y.mean())
    with np.errstate(divide='ignore', invalid='ignore'):
        f = (Vt @ beta) * s / uy
    return f, s


def ridge_coefficients(Z, y, lam):
    """Ridge solution on centred Z: (ZᵀZ + λI)⁻¹ Zᵀ (y - ȳ), via the SVD."""
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    return Vt.T @ ((s / (s**2 + lam)) * (U.T @ (y - y.mean())))


def pcr_coefficients(Z, y, k):
    """PCR solution with the first k components on centred Z."""
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    return Vt[:k].T @ ((U[:, :k].T @ (y - y.mean())) / s[:k])


def frequency_power(w, d=1.0):
    """Fraction of a weight spectrum's power above each frequency (cycles per unit of d).

    Returns
    -------
    f : ndarray
        Frequencies [cycles nm⁻¹].
    frac_above : ndarray
        Fraction of Σ|W(f)|² at frequencies ≥ f.
    """
    w = np.asarray(w, float) - np.mean(w)
    P = np.abs(np.fft.rfft(w * np.hanning(len(w))))**2
    f = np.fft.rfftfreq(len(w), d=d)
    c = np.cumsum(P[::-1])[::-1]
    return f, c / c[0]
