"""Tests for :mod:`epft_up.sdp.validate` (splits, baselines, scores)."""

import numpy as np
import pandas as pd
import pytest

from epft_up.sdp import validate as V


def test_replace_zeros_and_lod_proxy():
    np.testing.assert_allclose(V.replace_zeros([0.0, 0.2, -1.0], 0.004), [0.002, 0.2, 0.002])
    np.testing.assert_allclose(V.replace_zeros([0.0], 0.004, frac=0.1), [0.0004])
    pig = pd.DataFrame({'withzeros': [0.0, 0.003, 0.010], 'nozeros': [0.019, 0.5, 1.0]})
    lods = V.lod_proxy(pig)
    assert lods['withzeros'] == pytest.approx(0.003)
    assert lods['nozeros'] == V.PIGMENT_RESOLUTION


def test_oc4v6_value():
    wave = np.array([443.0, 490.0, 510.0, 555.0])
    Rrs = np.array([[0.002, 0.004, 0.003, 0.004]])      # max ratio = 1
    assert V.oc4v6(wave, Rrs)[0] == pytest.approx(10**0.3272)


def test_splits_are_partitions():
    sp = V.random_splits(20, n_perm=5, train_frac=0.75, seed=3)
    for tr, te in sp:
        assert len(tr) == 15 and len(te) == 5
        assert set(tr) | set(te) == set(range(20)) and not set(tr) & set(te)
    names, lo = V.loco_splits(np.array(['a', 'b', 'a', 'c', 'b']))
    assert names == ['a', 'b', 'c']
    np.testing.assert_array_equal(lo[0][1], [0, 2])
    res = [(te, np.full(len(te), float(i))) for i, (_, te) in enumerate(lo)]
    np.testing.assert_allclose(V.pooled_predictions(res, 5), [0, 1, 0, 2, 1])
    with pytest.raises(ValueError):
        V.pooled_predictions(res[:2], 5)


def test_loglinear_and_constant_fitters():
    rng = np.random.default_rng(0)
    x = 10**rng.uniform(-2, 0.5, 60)
    y = 0.3 * x**0.8
    f = V.LogLinearFitter(x, np.log10(y))
    tr, te = np.arange(40), np.arange(40, 60)
    np.testing.assert_allclose(f(tr, te), y[te], rtol=1e-10)
    flog = V.LogLinearFitter(x, np.log10(y), output='log')
    np.testing.assert_allclose(flog(tr, te), np.log10(y[te]), rtol=1e-10)
    with pytest.raises(ValueError):
        V.LogLinearFitter(np.r_[x[:-1], 0.0], np.log10(y))
    c = V.ConstantFitter(y)
    np.testing.assert_allclose(c(tr, te), y[tr].mean())


def test_pcr_fitter_clips_and_cross_validates():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(50, 12))
    y = X @ rng.normal(size=12)
    res = V.cross_validate(V.PCRFitter(X, y, 'pigment', seed=2, max_pcs=8),
                           V.random_splits(50, 3, seed=4))
    assert len(res) == 3 and all((pr >= 0).all() for _, pr in res)
    res_n = V.cross_validate(V.PCRFitter(X, y, None, seed=2, max_pcs=8),
                             V.random_splits(50, 3, seed=4))
    assert any((pr < 0).any() for _, pr in res_n)


def test_scores_and_paired():
    obs = np.array([1.0, 2.0, 3.0, 4.0])
    s = V.score_linear(obs, obs)
    assert s['R2'] == pytest.approx(1.0) and s['MAE'] == 0.0
    assert s['MADn'] == 0.0
    sl = V.score_log(np.log10(obs), np.log10(obs) + 0.1)
    assert sl['bias'] == pytest.approx(0.1) and sl['R2'] == pytest.approx(1.0)
    good = {'R2': np.array([0.8, 0.7]), 'MAE': np.array([1.0, 1.0])}
    bad = {'R2': np.array([0.5, 0.75]), 'MAE': np.array([2.0, 2.0])}
    p = V.paired(good, bad)
    assert p['R2']['frac_model_better'] == 0.5
    assert p['MAE']['frac_model_better'] == 1.0
    assert p['SS_MAE'] == pytest.approx(0.5)


def test_bootstrap_paired_detects_better_model():
    rng = np.random.default_rng(5)
    obs = rng.normal(size=200)
    good = obs + 0.1 * rng.normal(size=200)
    bad = obs + 1.0 * rng.normal(size=200)
    b = V.bootstrap_paired(obs, good, bad, V.score_linear, n_boot=300)
    assert b['dR2'] > 0 and b['dR2_ci'][0] > 0 and b['SS_MAE_ci'][0] > 0
    assert b['frac_boot_R2_better'] == 1.0


def _bench_fixture(seed=0, n=48):
    rng = np.random.default_rng(seed)
    camps = np.repeat(['a', 'b', 'c', 'd'], n // 4)
    tchla = 10**rng.uniform(-1.5, 0.5, n)
    pig = pd.DataFrame({'Tchla': tchla,
                        'Fuco': 0.2 * tchla**1.2 * 10**(0.05 * rng.normal(size=n)),
                        'Zea': np.where(rng.random(n) < 0.2, 0.0, 0.05 * 10**rng.normal(0, .2, n))})
    return pig, tchla, camps


def test_benchmark_targets_and_baselines():
    pig, tchla, camps = _bench_fixture()
    B = V.Benchmark(pig, tchla * 1.1, tchla * 0.9, camps, n_perm=5, seed=2)
    assert set(B.targets) == {'abs:Tchla', 'abs:Fuco', 'abs:Zea', 'ratio:Fuco', 'ratio:Zea'}
    assert B.loco_names == ['a', 'b', 'c', 'd']
    # a model that knows the truth beats the nulls on the Zea ratio (not Tchla-driven)
    truth = {t: v['y'] for t, v in B.targets.items()}

    def oracle_factory(y, constraint):
        return lambda tr, te: y[te]

    df = B.evaluate(oracle_factory, 'truth', targets=['ratio:Zea', 'abs:Fuco'])
    assert df.loc['ratio:Zea', 'R2_random'] == pytest.approx(1.0)
    assert bool(df.loc['ratio:Zea', 'beats_null_loco'])
    assert df.loc['ratio:Zea', 'RMS_loco'] == pytest.approx(0.0)
    assert {'R2_loco', 'dR2_loco_lo', 'frac_wins_random'} <= set(df.columns)
    del truth


def test_benchmark_null_model_scores_itself_as_no_gain():
    pig, tchla, camps = _bench_fixture(seed=3)
    B = V.Benchmark(pig, tchla, tchla, camps, n_perm=5, seed=2)

    def gsm_null_factory(y, constraint):
        return V.LogLinearFitter(tchla, np.log10(V.replace_zeros(y, B.lods['Fuco'])))

    df = B.evaluate(gsm_null_factory, 'nullcopy', targets=['abs:Fuco'])
    assert df.loc['abs:Fuco', 'dR2_random'] == pytest.approx(0.0, abs=1e-12)
    assert not bool(df.loc['abs:Fuco', 'beats_null_loco'])
