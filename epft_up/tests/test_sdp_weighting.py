"""Tests for :mod:`epft_up.sdp.weighting` (learned wavelength weightings)."""

import numpy as np
import pytest

from epft_up.sdp import weighting as W


def _lowrank_data(n=90, p=60, q=4, r=2, seed=0, noise=0.01):
    rng = np.random.default_rng(seed)
    wave = np.arange(p)
    # two smooth "latent" weight spectra
    L = np.column_stack([np.exp(-0.5 * ((wave - 20) / 6)**2),
                         np.exp(-0.5 * ((wave - 42) / 5)**2)])[:, :r]
    X = rng.normal(size=(n, p)) @ np.diag(np.linspace(1, 0.5, p))
    mix = rng.normal(size=(r, q))
    Y = X @ L @ mix + noise * rng.normal(size=(n, q))
    return X, Y, L, np.repeat(np.arange(5), n // 5 + 1)[:n]


def test_ridge_and_pls_fitters_learn_linear_signal():
    X, Y, _, _ = _lowrank_data()
    y = Y[:, 0]
    tr, te = np.arange(70), np.arange(70, 90)
    pr = W.RidgeGCVFitter(X, y, None)(tr, te)
    assert np.corrcoef(pr, y[te])[0, 1] > 0.95
    pl = W.PLSFitter(X, y, None, max_comp=8)(tr, te)
    assert np.corrcoef(pl, y[te])[0, 1] > 0.95
    pos = W.RidgeGCVFitter(X, y - y.min() + 0.01, 'pigment')(tr, te)
    assert np.all(pos >= 0)


def test_reduce_rank_exact_rank():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(30, 10))
    B = rng.normal(size=(10, 5))
    Br, Vr = W.PenalizedRRR.reduce_rank(X, B, 2)
    assert np.linalg.matrix_rank(X @ Br) == 2
    assert Vr.shape == (5, 2)


def test_penalized_rrr_selects_true_rank_and_predicts():
    X, Y, L, groups = _lowrank_data(r=2)
    m = W.PenalizedRRR(X, None, None, lam_noise=(0.0,), lam_smooth=(1e-6, 1e-4),
                       ranks=range(1, 5))
    idx = np.arange(len(Y))
    folds = np.array_split(np.random.default_rng(2).permutation(idx), 5)
    lam_n, lam_s, r = m.select(idx, Y, folds)
    assert r == 2
    fac = W.SharedRRRFitter(m, Y, groups, inner='group')
    tr, te = np.arange(70), np.arange(70, 90)
    preds = [fac(Y[:, j], None)(tr, te) for j in range(Y.shape[1])]
    for j, pr in enumerate(preds):
        assert np.corrcoef(pr, Y[te, j])[0, 1] > 0.95
    assert len(fac.chosen) == 1          # one joint fit, cached across targets


def test_smoothness_penalty_smooths_weights():
    X, Y, _, _ = _lowrank_data(r=1, q=1, noise=0.3)
    idx = np.arange(len(Y))
    rough = W.PenalizedRRR(X, None, None, lam_noise=(0.0,), lam_smooth=(1e-8,), ranks=(1,))
    smooth = W.PenalizedRRR(X, None, None, lam_noise=(0.0,), lam_smooth=(1.0,), ranks=(1,))
    Br, *_ = rough.fit_full(idx, Y, 0.0, 1e-8)
    Bs, *_ = smooth.fit_full(idx, Y, 0.0, 1.0)
    rough_d2 = np.sum(np.diff(Br[:, 0], 2)**2) / np.sum(Br[:, 0]**2)
    smooth_d2 = np.sum(np.diff(Bs[:, 0], 2)**2) / np.sum(Bs[:, 0]**2)
    assert smooth_d2 < 0.1 * rough_d2


def test_noise_penalty_downweights_noisy_bands():
    X, Y, _, _ = _lowrank_data(r=1, q=1)
    p = X.shape[1]
    Sn = np.diag(np.r_[np.ones(p // 2) * 1e-6, np.ones(p - p // 2)])  # 2nd half noisy
    idx = np.arange(len(Y))
    m = W.PenalizedRRR(X, None, Sn, ranks=(1,))
    B0, *_ = m.fit_full(idx, Y, 0.0, 1e-6)
    B1, *_ = m.fit_full(idx, Y, 10.0, 1e-6)
    share = lambda B: np.sum(B[p // 2:, 0]**2) / np.sum(B[:, 0]**2)  # noqa: E731
    assert share(B1) < share(B0)


def test_shared_fitter_rejects_unknown_target():
    X, Y, _, groups = _lowrank_data()
    fac = W.SharedRRRFitter(W.PenalizedRRR(X, None, None), Y, groups)
    with pytest.raises(ValueError):
        fac(np.zeros(len(Y)), None)


def test_sigma_abs_shrinks_weights_on_noisy_bands():
    X, Y, _, _ = _lowrank_data(r=1, q=1)
    p = X.shape[1]
    S = np.diag(np.r_[np.zeros(p // 2), np.ones(p - p // 2)])
    idx = np.arange(len(Y))
    m0 = W.PenalizedRRR(X, None, None, ranks=(1,))
    m1 = W.PenalizedRRR(X, None, None, ranks=(1,), Sigma_abs=S)
    B0, *_ = m0.fit_full(idx, Y, 0.0, 1e-6)
    B1, *_ = m1.fit_full(idx, Y, 0.0, 1e-6)
    assert np.sum(B1[p // 2:]**2) < 0.5 * np.sum(B0[p // 2:]**2)
