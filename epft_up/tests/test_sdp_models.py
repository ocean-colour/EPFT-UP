"""Tests for :mod:`epft_up.sdp.models` and the finite differences in
:mod:`epft_up.sdp.spectral`.

Tier 1 is synthetic; Tier 2 (``needs_kramer2022``) checks the Tchla row of
Kramer 2022 Table 2 against the Execution #2 δRrs product.
"""

import os

import numpy as np
import pytest

from epft_up.sdp import models, spectral
from epft_up.sdp import data as sdpdata
from epft_up.tests.conftest import needs_kramer2022


def test_difference_orders_and_step():
    w = np.arange(400.0, 421.0)
    y = 3e-6 * (w - 405.0)**2 + 1e-4 * w
    wd, d2 = spectral.difference(w, y, order=2)
    np.testing.assert_allclose(d2, 6e-6, rtol=1e-6)
    assert wd[0] == 401 and wd[-1] == 419 and d2.size == 19
    wd1, d1 = spectral.difference(w, y, order=1)
    assert wd1[0] == 400.5 and d1.size == 20
    wd5, d5 = spectral.difference(w, y, order=2, step=5)
    np.testing.assert_allclose(wd5, [405.0, 410.0, 415.0])
    np.testing.assert_allclose(d5, 6e-6 * 25, rtol=1e-6)
    _, d5s = spectral.difference(w, y, order=2, step=5, scale=True)
    np.testing.assert_allclose(d5s, 6e-6, rtol=1e-6)
    with pytest.raises(ValueError):
        spectral.difference(w, y, order=3)


def test_derivative_features_concat():
    w = np.arange(400.0, 411.0)
    y = np.vstack([np.sin(w / 3.0), np.cos(w / 3.0)])
    labels, X = spectral.derivative_features(w, y, orders=(1, 2))
    assert X.shape == (2, 10 + 9)
    assert labels[0] == 'd1_400.5' and labels[10] == 'd2_401'


def _synthetic(n=120, p=60, rank=6, seed=0):
    rng = np.random.default_rng(seed)
    basis = rng.normal(size=(rank, p))
    amp = rng.normal(size=(n, rank)) * np.array([5, 4, 3, 2, 1, 0.5])[:rank]
    X = amp @ basis + 0.01 * rng.normal(size=(n, p))
    w = np.linalg.pinv(basis) @ np.array([1.0, -0.5, 0.3, 0.0, 0.2, 0.1])[:rank]
    y = X @ w + 2.0
    return X, y


def test_nested_pc_shortcut_matches_explicit_ols():
    """b_l from one fit equals an explicit OLS refit on the first l scores."""
    rng = np.random.default_rng(3)
    X = rng.normal(size=(40, 25))
    y = rng.normal(size=40)
    Z, _, _ = models._zscore(X)
    _, _, Vt = np.linalg.svd(Z, full_matrices=False)
    scores = Z @ Vt.T
    b_all = (scores.T @ y) / np.sum(scores**2, axis=0)
    for l in (1, 5, 12, 24):
        A = np.column_stack([np.ones(40), scores[:, :l]])
        coef = np.linalg.lstsq(A, y, rcond=None)[0]
        assert coef[0] == pytest.approx(y.mean())
        np.testing.assert_allclose(coef[1:], b_all[:l], rtol=1e-8, atol=1e-10)


def test_train_pcr_recovers_linear_truth():
    X, y = _synthetic()
    r = models.train_pcr_kramer(X, y, n_perm=10, max_pcs=15, seed=7, constraint=None)
    s = r.summary()
    assert r.coefs.shape == (10, X.shape[1]) and r.intercepts.shape == (10,)
    assert s['R2'][0] > 0.999
    assert np.all(r.n_pcs <= 15) and np.all(r.n_pcs >= 1)
    med, runs = models.predict_ensemble(X, r.coefs, r.intercepts, constraint=None)
    assert runs.shape == (X.shape[0], 10)
    np.testing.assert_allclose(med, y, atol=0.05)
    # each validation set is the complement of a 75% training draw
    assert all(len(v) == X.shape[0] - int(0.75 * X.shape[0]) for v in r.valid_idx)


def test_train_pcr_deterministic_and_options():
    X, y = _synthetic(seed=1)
    a = models.train_pcr_kramer(X, y, n_perm=3, max_pcs=10, seed=5)
    b = models.train_pcr_kramer(X, y, n_perm=3, max_pcs=10, seed=5)
    np.testing.assert_array_equal(a.coefs, b.coefs)
    c = models.train_pcr_kramer(X, y, n_perm=3, max_pcs=10, seed=5, valid_scaling='train')
    assert c.settings['valid_scaling'] == 'train'
    with pytest.raises(ValueError):
        models.train_pcr_kramer(X, y, n_perm=1, valid_scaling='bogus')
    with pytest.raises(ValueError):
        models.train_pcr_kramer(np.where(X > 3, np.nan, X), y, n_perm=1)


def test_predict_ensemble_clip_and_lod():
    X = np.array([[1.0], [2.0], [3.0]])
    coefs = np.array([[1.0], [2.0], [3.0]])
    icpt = np.array([-2.5, -2.5, -2.5])
    med, runs = models.predict_ensemble(X, coefs, icpt)
    np.testing.assert_allclose(runs[:, 1], [-0.5, 1.5, 3.5])
    np.testing.assert_allclose(med, [0.0, 1.5, 3.5])        # median, clipped
    med_l, _ = models.predict_ensemble(X, coefs, icpt, lod=2.0)
    np.testing.assert_allclose(med_l, [0.0, 0.0, 3.5])


@needs_kramer2022
def test_tchla_matches_kramer_table2():
    prod = sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz'
    if not os.path.isfile(prod):
        pytest.skip('run scripts/sdp/reproduce_gsm.py first')
    z = np.load(prod)
    _, X = spectral.difference(z['wave'], z['dRrs'], order=2)
    y = sdpdata.load_kramer2022().pigments['Tchla'].to_numpy()
    s = models.train_pcr_kramer(X, y, seed=1).summary()
    # Kramer 2022 Table 2, Tchla: R² 0.72 ± 0.15; normalized MAD 0.498 ± 0.127
    assert abs(s['R2'][0] - 0.72) < 0.15
    assert abs(s['MAE_norm_pred'][0] - 0.498) < 0.127
