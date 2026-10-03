"""Tests for :mod:`epft_up.sdp.bayes` and the Execution #8 noise helpers."""

import numpy as np
import pytest

from epft_up.sdp import bayes, noise


def _data(n=120, p=40, seed=0, sigma=0.2):
    rng = np.random.default_rng(seed)
    wave = np.arange(p)
    L = np.exp(-0.5 * ((wave - 15) / 5)**2)
    X = rng.normal(size=(n, p))
    lat = X @ L
    Y = np.column_stack([lat, -0.5 * lat]) + sigma * rng.normal(size=(n, 2))
    groups = np.repeat(np.arange(6), n // 6)
    return X, Y, groups


def test_bayesian_shared_w_calibrated_on_synthetic():
    X, Y, g = _data()
    tr, te = np.arange(90), np.arange(90, 120)
    m = bayes.BayesianSharedW((0.0,), (1e-6, 1e-3), (1, 2, 3)).fit(X[tr], None, Y[tr], g[tr])
    mu, sd = m.predict(X[te], None)
    assert mu.shape == sd.shape == (30, 2)
    assert m.hp[2] == 1                              # one shared latent direction
    s = bayes.gaussian_scores(Y[te], mu, sd)
    assert 0.55 < s['cov68'] < 0.85 and s['cov95'] > 0.88
    rmse = np.sqrt(np.mean((Y[te] - mu)**2, axis=0))
    np.testing.assert_allclose(np.median(sd, axis=0), rmse, rtol=0.45)
    assert m.effective_weights().shape == (40, 2)


def test_insample_calibration_is_overconfident():
    """The latent directions are fitted to the targets, so in-sample residuals understate error."""
    X, Y, g = _data()
    tr, te = np.arange(90), np.arange(90, 120)
    kw = dict(lam_noise=(0.0,), lam_smooth=(1e-6, 1e-3), ranks=(1, 2, 3))
    ins = bayes.BayesianSharedW(**kw, calibration='insample').fit(X[tr], None, Y[tr], g[tr])
    oof = bayes.BayesianSharedW(**kw, calibration='oof').fit(X[tr], None, Y[tr], g[tr])
    sd_in = ins.predict(X[te], None)[1]
    sd_oof = oof.predict(X[te], None)[1]
    assert np.median(sd_in) < np.median(sd_oof)
    with pytest.raises(ValueError):
        bayes.BayesianSharedW(**kw, calibration='bogus')


def test_gaussian_scores_exact():
    y = np.zeros(10000)
    rng = np.random.default_rng(1)
    mu = rng.normal(size=10000)
    s = bayes.gaussian_scores(y, mu, np.ones(10000))
    assert s['cov68'] == pytest.approx(0.683, abs=0.02)
    assert s['cov95'] == pytest.approx(0.95, abs=0.01)
    assert s['z_sd'] == pytest.approx(1.0, abs=0.03)
    # CRPS of a perfect point at the mean of N(0,1): 1/sqrt(pi)*(sqrt(2)-1)
    s0 = bayes.gaussian_scores([0.0], [0.0], [1.0])
    assert s0['crps'] == pytest.approx((np.sqrt(2) - 1) / np.sqrt(np.pi), rel=1e-6)


@pytest.mark.parametrize('model', noise.MODELS)
def test_draws_match_covariance(model):
    w = np.arange(400.0, 461.0)
    R = np.linspace(0.006, 0.002, w.size)[None, :]
    rng = np.random.default_rng(2)
    d = noise.draw(w, R, model, 3000, rng)[:, 0, :] - R[0]
    C = noise.covariance(w, R[0], model)
    np.testing.assert_allclose(d.std(axis=0), np.sqrt(np.diag(C)), rtol=0.08)


def test_five_nm_grid_has_no_extra_smoothing():
    w5 = np.arange(400.0, 701.0, 5.0)
    C = noise.covariance(w5, np.full(w5.size, 0.003), 'pace_white')
    np.testing.assert_allclose(C, np.diag(np.diag(C)))
    with pytest.raises(ValueError):
        noise.band_sigma(w5, np.ones(w5.size), 'bogus')
