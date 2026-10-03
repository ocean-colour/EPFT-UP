"""Tests for :mod:`epft_up.sdp.theory` and :mod:`epft_up.sdp.noise`."""

import numpy as np
import pytest

from epft_up.sdp import noise, spectral, theory


def test_difference_matrix_matches_spectral():
    rng = np.random.default_rng(0)
    w = np.arange(400.0, 431.0)
    y = rng.normal(size=(3, w.size))
    for order in (1, 2):
        D = theory.difference_matrix(w.size, order)
        np.testing.assert_allclose(y @ D.T, spectral.difference(w, y, order)[1])
    D2 = theory.difference_matrix(w.size, 2)
    # null space: constants and straight lines
    np.testing.assert_allclose(D2 @ np.ones(w.size), 0)
    np.testing.assert_allclose(D2 @ w, 0, atol=1e-9)


def test_boxcar_matrix_matches_moving_mean_interior():
    rng = np.random.default_rng(1)
    y = rng.normal(size=40)
    M = theory.boxcar_matrix(40, 5)
    np.testing.assert_allclose((M @ y)[2:-2], spectral.moving_mean(y, 5)[2:-2])
    np.testing.assert_allclose((M @ y)[:2], y[:2])


def test_effective_weights_identity():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(10, 50))
    D = theory.difference_matrix(50, 2)
    A = rng.normal(size=(4, 48))
    W = theory.effective_weights(A, D)
    np.testing.assert_allclose((X @ D.T) @ A.T, X @ W.T, atol=1e-12)
    np.testing.assert_allclose(W.sum(axis=1), 0, atol=1e-10)   # blind to offsets


def test_implied_prior_is_rough():
    n = 128
    D = theory.difference_matrix(n, 2)
    C = theory.implied_prior_covariance(D, np.ones(n - 2))
    np.testing.assert_allclose(C, D.T @ D)
    evals = np.linalg.eigvalsh(C)
    assert evals.max() == pytest.approx(16.0, rel=0.01)   # (2 sin(pi f))^4 at Nyquist
    assert np.sum(evals < 1e-9) == 2                      # constant and linear


def test_white_noise_amplification():
    D = theory.difference_matrix(30, 2)
    C = theory.propagate_covariance(D, theory.white_covariance(np.ones(30)))
    np.testing.assert_allclose(np.diag(C), 6.0)


def test_dof_signal_limits_and_invariance():
    rng = np.random.default_rng(3)
    B = rng.normal(size=(20, 5))
    Cs = B @ B.T * 100.0                  # 5 strong signal directions
    Cn = np.eye(20)
    d, lam = theory.dof_signal(Cs, Cn)
    assert 4.5 < d < 5.0 and np.sum(lam > 1) == 5
    T = rng.normal(size=(20, 20))         # invertible transform: invariant
    dT, _ = theory.dof_signal(T @ Cs @ T.T, T @ Cn @ T.T)
    assert dT == pytest.approx(d, rel=1e-6)
    D = theory.difference_matrix(20, 2)   # non-invertible: can only decrease
    dD, _ = theory.dof_signal(D @ Cs @ D.T, D @ Cn @ D.T)
    assert dD <= d + 1e-9


def test_filter_factors_pcr_ridge():
    rng = np.random.default_rng(4)
    Z = rng.normal(size=(40, 15))
    Z -= Z.mean(0)
    y = Z @ rng.normal(size=15) + rng.normal(size=40)
    f, s = theory.filter_factors(Z, y, theory.pcr_coefficients(Z, y, 6))
    np.testing.assert_allclose(f, np.r_[np.ones(6), np.zeros(len(s) - 6)], atol=1e-10)
    lam = 3.0
    fr, s = theory.filter_factors(Z, y, theory.ridge_coefficients(Z, y, lam))
    np.testing.assert_allclose(fr, s**2 / (s**2 + lam), rtol=1e-10)


def test_noise_models():
    w = np.arange(400.0, 701.0)
    s = noise.pace_sigma(w)
    assert s.shape == w.shape and np.all(s > 0) and s[0] > s[-1]
    s1 = noise.pace_sigma(w, include_sampling=True)
    # ocpy rescales by sqrt(file spacing / grid spacing); PACE_error.csv is at ~2 nm
    np.testing.assert_allclose(s1 / s, np.sqrt(2.0), rtol=0.05)
    r = np.linspace(0.01, 1e-4, 301)
    p = noise.pct_floor_sigma(r, 0.02)
    np.testing.assert_allclose(p, 0.02 * np.maximum(r, np.median(r)))
    assert noise.rounding_sigma(3)[0] == pytest.approx(1e-6 / np.sqrt(12))
