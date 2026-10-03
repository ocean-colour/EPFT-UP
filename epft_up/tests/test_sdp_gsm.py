"""Tests for :mod:`epft_up.sdp.gsm` and :mod:`epft_up.sdp.spectral`.

Tier 1 uses synthetic reference tables and self-consistent synthetic spectra
(the GSM shapes depend on the spectrum's own band ratios, so a synthetic
spectrum is generated as a fixed point of the forward model). Tier 2
(``needs_kramer2022``) checks the real tables, the Python port's b_sw, and the
Kramer 2022 Fig. 4B statistics.
"""

import numpy as np
import pytest

from epft_up.sdp import gsm, spectral
from epft_up.sdp import data as sdpdata
from epft_up.tests.conftest import needs_kramer2022

WAVE = np.arange(400.0, 701.0)


def _synthetic_tables(wave=WAVE):
    """Smooth stand-ins for A(λ), B(λ), a_w with the right magnitudes."""
    A = 0.02 + 0.03 * np.exp(-0.5 * ((wave - 440) / 30)**2) \
        + 0.015 * np.exp(-0.5 * ((wave - 675) / 10)**2)
    B = 0.75 + 0.2 * (wave - 400) / 300
    aw = 0.005 + 0.5 * np.exp((wave - 700) / 60)
    return gsm.RefTables(wave=wave, A=A, B=B, aw=aw, source={'synthetic': True})


def _synthetic_spectrum(params, T=18.0, S=35.0, tables=None, n_iter=50):
    """Above-water Rrs that is a fixed point of the GSM-like forward model."""
    tables = tables or _synthetic_tables()
    Rrs = np.full((1, WAVE.size), 0.004)
    for _ in range(n_iter):
        sh = gsm.build_shapes(WAVE, Rrs, [T], [S], tables=tables)
        rrs = gsm.forward_rrs(params, sh.aw, sh.A, sh.B, sh.bbw[0], sh.adg_star[0],
                              sh.bbp_star[0])
        new = gsm.rrs_to_Rrs(rrs)[None, :]
        if np.max(np.abs(new - Rrs)) < 1e-15:
            break
        Rrs = new
    return Rrs


def test_rrs_conversion_roundtrip():
    Rrs = np.linspace(1e-5, 0.03, 50)
    np.testing.assert_allclose(gsm.rrs_to_Rrs(gsm.Rrs_to_rrs(Rrs)), Rrs, rtol=1e-12)


def test_bsw_zhang2009_physics():
    w = np.array([400.0, 500.0, 600.0])
    b35 = gsm.bsw_zhang2009(w, 20.0, 35.0)
    b0 = gsm.bsw_zhang2009(w, 20.0, 0.0)
    assert b35.shape == (3,)
    assert 2.4e-3 < b35[1] < 2.7e-3                  # seawater b_sw(500) [m^-1]
    assert np.all((b35 / b0 > 1.25) & (b35 / b0 < 1.35))   # ~30% salt enhancement
    slope = np.log(b35[0] / b35[2]) / np.log(400 / 600)
    assert -4.35 < slope < -4.15
    # vectorised over samples
    bv = gsm.bsw_zhang2009(w, np.array([20.0, 20.0]), np.array([35.0, 0.0]))
    np.testing.assert_allclose(bv, np.vstack([b35, b0]))


def test_forward_model_limits():
    t = _synthetic_tables()
    one = np.ones_like(WAVE)
    # no particles / no dissolved: u = bbw / (aw + bbw)
    r = gsm.forward_rrs((1e-12, 0.0, 0.0), t.aw, t.A, t.B, 1e-3 * one, one, one)
    u = 1e-3 / (t.aw + 1e-3)
    np.testing.assert_allclose(r, gsm.G0 * u + gsm.G1 * u**2, rtol=1e-6)
    # more chlorophyll -> less blue reflectance
    lo = gsm.forward_rrs((0.05, 0.01, 0.002), t.aw, t.A, t.B, 1e-3 * one, one, one)
    hi = gsm.forward_rrs((2.0, 0.01, 0.002), t.aw, t.A, t.B, 1e-3 * one, one, one)
    assert hi[40] < lo[40]


def test_sdg_sign_follows_carder():
    # Carder et al. (1999): positive slope magnitude ~0.015 nm^-1, decaying a_dg
    s = gsm.sdg_slope(0.006, 0.002)
    assert s == pytest.approx(0.01447 + 0.00033 * 3)
    sh = gsm.build_shapes(WAVE, np.full((1, WAVE.size), 0.004), [18.0], [35.0],
                          tables=_synthetic_tables())
    assert sh.adg_star[0, 0] > 1 > sh.adg_star[0, -1]
    assert sh.adg_star[0, 43] == pytest.approx(1.0)        # 443 nm


@pytest.mark.parametrize('method', ['kramer', 'lsq'])
def test_fit_recovers_synthetic_truth(method):
    tables = _synthetic_tables()
    truth = np.array([0.3, 0.02, 0.0025])
    Rrs = _synthetic_spectrum(truth, tables=tables)
    fit = gsm.fit_gsm(WAVE, Rrs, 18.0, 35.0, tables=tables, method=method)
    assert fit.converged.all()
    np.testing.assert_allclose(fit.params[0], truth, rtol=1e-3)
    assert np.max(np.abs(fit.dRrs)) < 1e-7
    assert fit.settings['method'] == method


def test_fit_rejects_unknown_method():
    with pytest.raises(ValueError):
        gsm.fit_gsm(WAVE, np.full((1, WAVE.size), 0.004), 18.0, 35.0,
                    tables=_synthetic_tables(), method='bogus')


def test_moving_mean_and_trim():
    y = 3.0 + 0.5 * np.arange(20.0)
    m = spectral.moving_mean(y, 5)
    assert np.isnan(m[:2]).all() and np.isnan(m[-2:]).all()
    np.testing.assert_allclose(m[2:-2], y[2:-2])        # linear is preserved
    w, t = spectral.trim_edges(np.arange(20.0), m, 4)
    assert w[0] == 4 and w[-1] == 15 and np.isfinite(t).all()
    with pytest.raises(ValueError):
        spectral.moving_mean(y, 4)


def test_kramer_preprocess_grid():
    wave = np.arange(390.0, 712.0, 3.3)
    Rrs = np.vstack([0.004 + 1e-6 * (wave - 400.0)])
    out = spectral.kramer_preprocess(wave, Rrs)
    assert out.shape == (1, 301)
    np.testing.assert_allclose(out[0], 0.004 + 1e-6 * (WAVE - 400.0), atol=1e-9)
    with pytest.raises(ValueError):
        spectral.kramer_preprocess(wave, Rrs, width=11, n_trim=4)


@needs_kramer2022
def test_real_reference_tables():
    t = gsm.load_ref_tables()
    assert t.wave.size == 301
    assert t.A[43] == pytest.approx(0.0501146) and t.B[43] == pytest.approx(0.75803)
    assert t.aw[0] == pytest.approx(0.00222)
    assert t.source['repo_sha'] == sdpdata.REF_REPOS['Rrs_pigments']['sha']


@needs_kramer2022
def test_reproduces_kramer_fig4b():
    d = sdpdata.load_kramer2022()
    fit = gsm.fit_gsm(d.wave, d.Rrs, d.meta['temp'], d.meta['sal'])
    assert fit.converged.all()
    x = np.log10(d.pigments['Tchla'].to_numpy())
    y = np.log10(fit.params[:, 0])
    worst = np.argmax(np.abs(y - x))
    assert d.meta['campaign'][worst] == 'SABOR'
    keep = np.arange(d.n) != worst
    r2 = np.corrcoef(x[keep], y[keep])[0, 1]**2
    slope, icpt = np.polyfit(x[keep], y[keep], 1)
    # Kramer 2022 Fig. 4B: y = 0.96x - 0.093, R^2 = 0.86
    assert r2 == pytest.approx(0.86, abs=0.01)
    assert slope == pytest.approx(0.96, abs=0.01)
    assert icpt == pytest.approx(-0.093, abs=0.01)


def test_spline_residual_removes_smooth_keeps_bump():
    w = WAVE
    smooth = 0.004 * np.exp(-(w - 400) / 120)
    bump = 2e-5 * np.exp(-0.5 * ((w - 550) / 4)**2)
    r = spectral.spline_residual(w, np.vstack([smooth, smooth + bump]))
    # the smooth part is removed to < 0.5% of the spectrum (M1 is a heavy spline)
    assert np.max(np.abs(r[0])) < 0.005 * smooth.max()
    assert r[1][150] - r[0][150] > 0.5 * 2e-5        # narrow bump retained


def test_savgol_second_derivative_quadratic():
    w = WAVE
    d2 = spectral.savgol_second_derivative(w, 3e-6 * (w - 500.0)**2, 11)
    np.testing.assert_allclose(d2[0, 10:-10], 6e-6, rtol=1e-6)


def test_smoothed_tables_flatten_structure():
    t = _synthetic_tables()
    s = gsm.smoothed_tables(t, 80)
    assert s.source['smoothed_fwhm_nm'] == 80
    assert np.all(s.A > 0)
    # the narrow 675 nm feature is strongly reduced, the broad mean kept
    assert s.A[275] < t.A[275]
    assert abs(np.mean(np.log(s.A)) - np.mean(np.log(t.A))) < 0.05
