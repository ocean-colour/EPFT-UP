"""Tests for the MOANA track (`epft_up/moana/`; design doc §8).

Tier 1 (always run): the vendored LUTs and their pinned hashes, the retrieval
invariants, the PC-mapping flag, NASA-compat post-processing, every pipeline
stage on fabricated streams with planted answers, the training recipe on a
planted linear model, and the pandas time-resolution regression.

Tier 2 (`needs_amt24`): the real AMT24 delivery — .sav layout, the ES≢LT
guard, the BODC flow-cytometry code mapping, and one real day through the
full Level-2 chain. Earthdata tests (`needs_netrc`) arrive with validation
target (iii); none exist yet.
"""

import hashlib

import numpy as np
import pandas as pd
import pytest

from epft_up.moana import io as mio
from epft_up.moana import algorithm as malg
from epft_up.moana import pipeline as mpipe
from epft_up.moana import train as mtrain
from epft_up.tests.conftest import needs_amt24

# ---------------------------------------------------------------------------
# LUTs (design §8 "LUT tests")
# ---------------------------------------------------------------------------

#: sha256 pinned at vendoring time (epft_up/data/moana/README.md); a change
#: means NASA moved the algorithm from under us and every result needs rerun.
LUT_SHA256 = {
    'pca_picophyto.h5':
        '478ace8798ce525bc6a09567b3773cadd532056c13d20ce6daf0a529f46deb57',
    'picophyt.json':
        'cd6b5c1cca879bc97acf6c5a3d4f47bf658d5a839b99ff4b1f534eac881756d4',
}


@pytest.fixture(scope='module')
def lut():
    return mio.load_luts()


def test_lut_files_pinned():
    for name, want in LUT_SHA256.items():
        got = hashlib.sha256((mio.LUT_DIR / name).read_bytes()).hexdigest()
        assert got == want, f'{name} changed since vendoring'


def test_lut_shapes_and_grid(lut):
    assert lut['wave'].shape == (124,)
    assert lut['wave'][0] == 414 and lut['wave'][-1] == 660
    assert np.all(np.diff(lut['wave']) == 2)
    assert lut['V'].shape == (124, 45)
    assert lut['npc'] == 17
    assert lut['pro_coef'].size == lut['npc'] + 2   # intercept + SST + U1..17
    assert lut['syn_coef'].size == lut['npc'] + 1
    assert lut['apeuk_coef'].size == lut['npc'] + 1


def test_lut_orthonormal(lut):
    VtV = lut['V'].T @ lut['V']
    assert np.abs(VtV - np.eye(45)).max() < 1e-6  # float32 storage precision


# ---------------------------------------------------------------------------
# Retrieval invariants and steps (design §3.1)
# ---------------------------------------------------------------------------

def _synthetic_rrs(n=5, seed=0):
    """Smooth, positive, blue-peaked spectra on a 400-700 @ 2 nm grid."""
    rng = np.random.default_rng(seed)
    wave = np.arange(400.0, 701.0, 2.0)
    amp = rng.uniform(0.005, 0.02, size=(n, 1))
    rrs = amp * np.exp(-((wave - 440.0) / 80.0) ** 2) + 0.002
    return wave, rrs + rng.normal(scale=1e-5, size=(n, wave.size))


def test_standardize_invariants():
    _, rrs = _synthetic_rrs()
    std = malg.standardize(rrs[:, :124])
    # Every spectrum: mean exactly 0, L2 norm exactly sqrt(N-1) (design §3.1).
    assert np.allclose(std.mean(axis=1), 0, atol=1e-12)
    assert np.allclose(np.linalg.norm(std, axis=1), np.sqrt(123), atol=1e-9)


def test_standardize_matches_formula():
    # Independent hand computation of Lange Eq. 3 with the N-1 denominator.
    rng = np.random.default_rng(3)
    x = rng.uniform(0.001, 0.02, size=124)
    manual = (x - x.mean()) / x.std(ddof=1)
    assert np.allclose(malg.standardize(x)[0], manual)


def test_scores_bounded(lut):
    wave, rrs = _synthetic_rrs(n=20)
    out = malg.run_moana(wave, rrs, sst=25.0, lut=lut)
    # |U_i| <= sqrt(123) for an orthonormal basis (design §3.1).
    assert np.abs(out['scores']).max() <= np.sqrt(123) + 1e-9
    assert np.all((out['scores'] ** 2).sum(axis=1) <= 123 + 1e-6)


def test_interp_flags(lut):
    wave, rrs = _synthetic_rrs(n=1)
    # Full coverage: no flag.
    _, flags = malg.interp_to_moana(wave, rrs, lut['wave'])
    assert flags[0] == 0
    # Input stops at 600 nm: extrapolation flag.
    keep = wave <= 600
    _, flags = malg.interp_to_moana(wave[keep], rrs[:, keep], lut['wave'])
    assert flags[0] & malg.FLAG_EXTRAPOLATED
    # Nearly all bands invalid: too-few-bands flag, NaN output.
    sparse = np.full_like(rrs, np.nan)
    sparse[0, :5] = rrs[0, :5]
    out, flags = malg.interp_to_moana(wave, sparse, lut['wave'])
    assert flags[0] & malg.FLAG_TOO_FEW_BANDS
    assert np.isnan(out[0]).all()


def test_evaluate_flags_and_positivity(lut):
    wave, rrs = _synthetic_rrs(n=3)
    out = malg.run_moana(wave, rrs, sst=None, lut=lut)
    # No SST: Pro is NaN and flagged, the log10 taxa are finite and positive.
    assert np.isnan(out['pro']).all()
    assert np.all(out['flags'] & malg.FLAG_BAD_SST)
    assert np.all(out['syn'] > 0) and np.all(out['apeuk'] > 0)
    # Non-positive SST is flagged per element.
    out = malg.run_moana(wave, rrs, sst=np.array([25.0, -1.0, 20.0]), lut=lut)
    assert not out['flags'][0] & malg.FLAG_BAD_SST
    assert out['flags'][1] & malg.FLAG_BAD_SST


def test_pc_mapping_moves_only_disputed_slots(lut):
    op = malg._coef_matrix(lut, 'operational')
    at = malg._coef_matrix(lut, 'atbd')
    # apeuk agrees everywhere (report §7.1).
    assert np.array_equal(op['apeuk'][2], at['apeuk'][2])
    # Pro: the U7 weight moves to U17, all else equal.
    w_op, w_at = op['pro'][2], at['pro'][2]
    assert w_at[17 - 1] == w_op[7 - 1] and w_at[7 - 1] == 0
    keep = np.ones(45, bool)
    keep[[7 - 1, 17 - 1]] = False
    assert np.array_equal(w_op[keep], w_at[keep])
    # Syn: U16 -> U13 likewise.
    w_op, w_at = op['syn'][2], at['syn'][2]
    assert w_at[13 - 1] == w_op[16 - 1] and w_at[16 - 1] == 0
    with pytest.raises(ValueError):
        malg._coef_matrix(lut, 'nonsense')


def test_mappings_agree_when_disputed_scores_vanish(lut):
    # Build rrs' orthogonal to PC7/PC13/PC16/PC17: both mappings must agree.
    rng = np.random.default_rng(7)
    coefs = rng.normal(size=6)
    keep_pcs = lut['V'][:, [0, 1, 2, 3, 4, 5]]  # span of undisputed PCs only
    base = keep_pcs @ coefs
    rrs_std = (base - base.mean()) / base.std(ddof=1)
    scores = malg.project(rrs_std, lut['V'])
    a = malg.evaluate(scores, sst=22.0, lut=lut, pc_mapping='operational')
    b = malg.evaluate(scores, sst=22.0, lut=lut, pc_mapping='atbd')
    for k in ('pro', 'syn', 'apeuk'):
        assert np.allclose(a[k], b[k])


def test_nasa_compat_clamps_and_truncates(lut):
    wave, rrs = _synthetic_rrs(n=2)
    # SST ~ 0.1 C drives log10(SST) large-negative -> Pro deeply negative.
    raw = malg.run_moana(wave, rrs, sst=0.1, lut=lut)
    assert np.all(raw['pro'] < 0)                        # raw floats kept
    assert np.all(raw['flags'] & malg.FLAG_NEGATIVE_PRO)  # ... and flagged
    compat = malg.run_moana(wave, rrs, sst=0.1, lut=lut, nasa_compat=True)
    assert np.all(compat['pro'] == 0)                    # OCSSW clamp
    # Truncation, not rounding: fractional parts drop toward zero.
    assert np.array_equal(compat['syn'], np.trunc(raw['syn']))
    assert np.array_equal(compat['apeuk'], np.trunc(raw['apeuk']))


def test_reconstruction_residual_detects_out_of_basis(lut):
    # A spectrum built inside the basis has ~zero residual; a sharp line
    # (nothing like the training set) does not.
    inside = lut['V'][:, :3] @ np.array([2.0, -1.0, 0.5])
    inside = (inside - inside.mean()) / inside.std(ddof=1)
    spike = np.zeros(124)
    spike[60] = 1.0
    spike = (spike - spike.mean()) / spike.std(ddof=1)
    r = malg.reconstruction_residual(np.vstack([inside, spike]), lut['V'])
    assert r[0] < 1e-6
    assert r[1] > 1.0


# ---------------------------------------------------------------------------
# Pipeline stages on fabricated streams (design §8 "pipeline tests")
# ---------------------------------------------------------------------------

def _fake_day(n=8):
    """Minimal day dict for the screening/selection stages."""
    return {
        'wave': mio.HSAS_WAVE.copy(),
        'time': np.linspace(12.0, 12.0 + (n - 1) / 3600.0 * 4, n),
        'lat': np.zeros(n), 'lon': np.zeros(n),
        'pitch': np.zeros(n), 'roll': np.zeros(n),
        'solar_zenith': np.full(n, 40.0),
        'delta_azimuth': np.full(n, 90.0),
        'lt': np.ones((n, mio.HSAS_WAVE.size)),
        'li': np.ones((n, mio.HSAS_WAVE.size)),
        'es': np.ones((n, mio.HSAS_WAVE.size)),
        'sst': np.full(n, 20.0),
    }


def test_tilt_from_pitch_roll():
    assert mpipe.tilt_from_pitch_roll(0.0, 0.0) == 0.0
    assert np.isclose(mpipe.tilt_from_pitch_roll(3.0, 0.0), 3.0)
    assert np.isclose(mpipe.tilt_from_pitch_roll(0.0, 4.0), 4.0)
    # Combined tilt exceeds either axis alone.
    assert mpipe.tilt_from_pitch_roll(3.0, 4.0) > 4.0


def test_screen_geometry_thresholds_are_lange_exact():
    day = _fake_day(n=7)
    day['pitch'] = np.array([0.0, 5.1, 0.0, 0.0, 0.0, 0.0, 0.0])
    day['solar_zenith'] = np.array([40, 40, 9.0, 81.0, 40, 40, 40.0])
    day['delta_azimuth'] = np.array([90, 90, 90, 90, 49.0, 171.0, -90.0])
    keep = mpipe.screen_geometry(day)
    # Row 0 passes; 1 fails tilt; 2-3 fail zenith; 4-5 fail azimuth;
    # 6 (negative dphi, |.|=90) passes — the screen is on the magnitude.
    assert keep.tolist() == [True, False, False, False, False, False, True]


def test_select_glint_minima_takes_darkest_per_minute():
    day = _fake_day(n=6)
    # Two 1-minute intervals of three spectra each.
    day['time'] = np.array([12.0, 12.005, 12.010, 12.020, 12.025, 12.028])
    nir = day['wave'] >= 750
    for i, lvl in enumerate([3.0, 1.0, 2.0, 5.0, 4.0, 6.0]):
        day['lt'][i, nir] = lvl
    idx = mpipe.select_glint_minima(day, np.ones(6, bool))
    assert idx.tolist() == [1, 4]   # the NIR-darkest of each minute
    # Screened-out spectra never win.
    keep = np.ones(6, bool)
    keep[1] = False
    assert mpipe.select_glint_minima(day, keep).tolist() == [2, 4]


def test_fit_glint_recovers_planted_parameters():
    rng = np.random.default_rng(11)
    wave = mio.HSAS_WAVE
    n = 4
    rho_true = np.array([0.02, 0.05, 0.08, 0.03])
    lnir_true = np.array([0.01, 0.0, 0.02, 0.005])
    li = 5.0 + 3.0 * np.exp(-((wave - 450) / 120.0) ** 2)[None, :] \
        + rng.uniform(0, 0.1, size=(n, wave.size))
    water = 0.02 * np.exp(-((wave - 500) / 60.0) ** 2)   # ~0 in the NIR
    lt = rho_true[:, None] * li + lnir_true[:, None] + water[None, :]
    rho, lnir = mpipe.fit_glint(lt, li, wave)
    assert np.allclose(rho, rho_true, atol=2e-3)
    assert np.allclose(lnir, lnir_true, atol=2e-3)


def test_qc_spectra_each_screen(lut):
    wave_n = mio.HSAS_WAVE
    wave_m = lut['wave']
    n = 4
    rrs_n = np.full((n, wave_n.size), 0.005)
    rrs_m = np.full((n, wave_m.size), 0.005)
    time = np.full(n, 12.0)          # noon UTC at lon 0 -> local noon
    lon = np.zeros(n)
    # Row 1: outside the 09-17 local window (via longitude, not clock).
    lon[1] = -179.0                  # local ~00:04
    # Row 2: a negative dip in the visible (native grid).
    rrs_n[2, np.argmin(np.abs(wave_n - 550))] = -1e-4
    # Row 3: a second-derivative spike in 610-660 on the 2 nm grid.
    rrs_m[3, np.argmin(np.abs(wave_m - 640))] += 10 * 2e-4 * 2.0 ** 2
    ok = mpipe.qc_spectra(rrs_n, rrs_m, wave_n, wave_m, time, lon)
    assert ok.tolist() == [True, False, False, False]


def test_resample_matches_np_interp(lut):
    wave, rrs = _synthetic_rrs(n=2)
    out = mpipe.resample_to_moana(rrs, wave, lut['wave'])
    assert out.shape == (2, 124)
    assert np.allclose(out[0], np.interp(lut['wave'], wave, rrs[0]))


def test_time_ns_is_resolution_proof():
    # Regression for the pandas>=2 trap found in the prompt-12 smoke: the
    # same instant expressed at second- and nanosecond-resolution must give
    # identical ns values (a bare astype('int64') differs by 1e9).
    t_s = pd.Series(pd.to_datetime(['2014-10-07 10:00']).as_unit('s'))
    t_ns = pd.Series(pd.to_datetime(['2014-10-07 10:00']).as_unit('ns'))
    assert mpipe._time_ns(t_s)[0] == mpipe._time_ns(t_ns)[0]


def test_attach_sst_interpolates_and_masks():
    sst = pd.DataFrame({
        'datetime': pd.date_range('2014-10-01', periods=11, freq='1h'),
        'sst': np.linspace(20.0, 21.0, 11),
    })
    q = pd.Series(pd.to_datetime(
        ['2014-10-01 05:00', '2014-10-01 05:30', '2014-09-30 00:00']))
    out = mpipe.attach_sst(q, sst)
    assert np.isclose(out[0], 20.5)
    assert np.isclose(out[1], 20.55)
    assert np.isnan(out[2])          # outside the series span


def test_build_training_matrix_median_and_strict():
    # A tiny stream: 40 spectra one minute apart, constant except a glint
    # spike the median must reject; one FCM sample in the middle.
    n = 40
    times = pd.date_range('2014-10-07 10:00', periods=n, freq='1min')
    rrs = np.full((n, 124), 0.004)
    rrs[20] = 0.4                      # the spike, exactly at the sample time
    rec = pd.DataFrame({'datetime': times, 'sst_hsas': 20.0})
    stream = {'records': rec, 'rrs': rrs}
    fcm = pd.DataFrame({
        'station': ['S1'], 'datetime': [times[20]],
        'lat': [0.0], 'lon': [0.0], 'depth': [5.0],
        'pro': [1e5], 'syn': [1e3], 'peuk': [5e2],
    })
    sst = pd.DataFrame({'datetime': times, 'sst': np.full(n, 21.0)})
    out = mpipe.build_training_matrix(stream, fcm, sst)
    assert len(out['matchups']) == 1
    assert out['matchups']['n_in_bin'].iloc[0] == 31    # +/-15 min inclusive
    assert np.allclose(out['rrs'][0], 0.004)            # median beat the spike
    assert np.isclose(out['matchups']['sst'].iloc[0], 21.0)
    # Lange-strict mode takes the single nearest spectrum -> the spike.
    strict = mpipe.build_training_matrix(stream, fcm, sst,
                                         config={'match_window_min': 0})
    assert np.allclose(strict['rrs'][0], 0.4)
    assert strict['matchups']['n_in_bin'].iloc[0] == 1


# ---------------------------------------------------------------------------
# Training (design §6)
# ---------------------------------------------------------------------------

def test_fit_pca_orthonormal_and_ordered():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(50, 124))
    pca = mtrain.fit_pca(X, center=False)
    V = pca['V']
    assert np.abs(V.T @ V - np.eye(V.shape[1])).max() < 1e-9
    assert np.all(np.diff(pca['sdev']) <= 1e-12)         # descending
    # center=True subtracts the column mean it reports.
    pca_c = mtrain.fit_pca(X, center=True)
    assert np.allclose(pca_c['mean'], X.mean(axis=0))


def test_truncate_pcs():
    sdev = np.array([10.0, 5.0, 0.02, 0.009, 0.001])
    # cut at 0.1% of PC1 (=0.01): keeps the first three.
    assert mtrain.truncate_pcs(sdev, frac=1e-3) == 3


def test_backward_stepwise_keeps_planted_predictors():
    rng = np.random.default_rng(4)
    n = 80
    X = rng.normal(size=(n, 8))
    y = 2.0 + 1.5 * X[:, 0] - 2.0 * X[:, 3] + rng.normal(scale=0.05, size=n)
    names = [f'U{j+1}' for j in range(8)]
    sel = mtrain.backward_stepwise_aic(X, y, names)
    assert set(sel['names']) >= {'U1', 'U4'}             # truth retained
    assert len(sel['names']) <= 4                        # noise mostly gone
    assert sel['fit']['r2_adj'] > 0.99


def test_train_moana_recovers_planted_model(lut):
    rng = np.random.default_rng(5)
    n = 80
    base = 0.008 * np.exp(-((lut['wave'] - 450) / 90.0) ** 2) + 0.001
    amps = rng.normal(scale=[3.0, 1.0, 0.5], size=(n, 3))
    rrs = base + 0.001 * amps @ lut['V'][:, :3].T \
        + rng.normal(scale=2e-5, size=(n, 124))
    sst = rng.uniform(15, 28, n)
    U = malg.standardize(rrs) @ lut['V'][:, :3]
    counts = {
        'pro': 2e5 + 8e4 * np.log10(sst) + 1e4 * U[:, 0]
               + rng.normal(scale=5e3, size=n),
        'syn': 10 ** (3.5 + 0.15 * U[:, 0] - 0.10 * U[:, 1]
                      + rng.normal(scale=0.03, size=n)),
        'peuk': 10 ** (3.0 + 0.12 * U[:, 1] + 0.20 * U[:, 2]
                       + rng.normal(scale=0.03, size=n)),
    }
    model = mtrain.train_moana(rrs, counts, sst=sst)
    assert 'logSST' in model['models']['pro']['names']
    assert model['models']['pro']['response'] == 'linear'
    assert model['models']['syn']['response'] == 'log10'
    for tax in ('pro', 'syn', 'peuk'):
        assert model['models'][tax]['r2_adj'] > 0.5
    # The per-taxon SST switch (paper Eq. 7 variant, design §4.3): SST must
    # have entered peuk's candidate pool — retained or explicitly dropped —
    # and must be absent entirely under the default (operational) setting.
    m3 = mtrain.train_moana(rrs, counts, sst=sst, use_sst={'peuk': True})
    assert 'logSST' in (m3['models']['peuk']['names']
                        + m3['models']['peuk']['dropped'])
    assert 'logSST' not in (model['models']['peuk']['names']
                            + model['models']['peuk']['dropped'])
    # Saturation guard: center=True explodes the PC count but must not crash.
    m2 = mtrain.train_moana(rrs, counts, sst=sst, center=True)
    assert m2['n_pc'] <= n - 12


def test_train_moana_raises_on_starved_taxon(lut):
    wave, rrs = _synthetic_rrs(n=20)
    rrs124 = mpipe.resample_to_moana(rrs, wave, lut['wave'])
    counts = {'pro': np.full(20, np.nan),                # nothing usable
              'syn': np.full(20, 1e3), 'peuk': np.full(20, 5e2)}
    with pytest.raises(ValueError, match='pro'):
        mtrain.train_moana(rrs124, counts, sst=np.full(20, 20.0))


def test_compare_loadings_identity(lut):
    cmp = mtrain.compare_loadings(lut['V'], lut['V'], n_compare=10)
    assert np.array_equal(cmp['best_match'], np.arange(10))
    assert np.allclose(cmp['best_cos'], 1.0, atol=1e-6)


# ---------------------------------------------------------------------------
# Validation metrics (prompt 14)
# ---------------------------------------------------------------------------

def test_seegers_metrics_known_values():
    from epft_up.moana.validation import seegers_metrics
    obs = np.array([100.0, 1000.0, 10000.0])         # varied, so R2 defined
    m = seegers_metrics(2 * obs, obs)                # uniformly 2x high
    assert np.isclose(m['bias'], 2.0) and np.isclose(m['mae'], 2.0)
    assert m['frac_unphysical'] == 0.0
    # A perfect prediction: bias = MAE = 1, R2 = 1.
    m = seegers_metrics(obs, obs)
    assert np.isclose(m['bias'], 1.0) and np.isclose(m['mae'], 1.0)
    assert np.isclose(m['r2'], 1.0)


def test_seegers_metrics_censors_clipped():
    from epft_up.moana.validation import seegers_metrics
    obs = np.array([100.0, 300.0, 900.0, 100.0, 100.0, 100.0])
    pred = np.array([100.0, 300.0, 900.0, 100.0, 0.0, -5.0])  # 2 unphysical
    m = seegers_metrics(pred, obs)
    # Censored values leave the log-space scores untouched (Q&A #13) ...
    assert m['n'] == 4 and np.isclose(m['mae'], 1.0)
    # ... and surface as the headline unphysical fraction.
    assert np.isclose(m['frac_unphysical'], 2 / 6)
    assert m['n_censored'] == 2


# ---------------------------------------------------------------------------
# Tier 2 — the real AMT24 delivery
# ---------------------------------------------------------------------------

@needs_amt24
def test_read_hsas_day_layout_and_es_guard():
    day = mio.read_hsas_day(mio.hsas_day_paths()[0])
    n = day['time'].size
    assert day['wave'].shape == (141,)
    assert day['wave'][0] == 306.0 and day['wave'][-1] == 796.0
    for key in ('lt', 'li', 'es', 'rrs_provider'):
        assert day[key].shape == (n, 141)
    # The ES bug is a CSV-only defect: the .sav's ES must differ from LT ...
    assert not np.array_equal(day['es'], day['lt'])
    # ... and be irradiance-magnitude (much brighter than the radiances).
    assert np.nanmedian(day['es']) > 10 * np.nanmedian(day['lt'])
    # The provider Rrs is the fixed-rho combination (prompt-8 finding).
    rho = (day['lt'] - day['rrs_provider'] * day['es']) / day['li']
    assert np.isclose(np.nanmedian(rho), 0.0280, atol=1e-3)


@needs_amt24
def test_fcm_mapping_and_surface_rule():
    fcm = mio.load_fcm()
    assert {'pro', 'syn', 'peuk'} <= set(fcm.columns)
    assert len(fcm) == 814                    # the CTD-only deposit (Q&A #35)
    # AMT24 tropical surface waters: Pro >> Syn on average — the mapping
    # would fail this instantly if P700/P701 were swapped.
    surf = mio.surface_fcm(fcm)
    assert (surf['depth'] <= 10.0).all()
    assert surf['station'].is_unique
    tropics = surf[np.abs(surf['lat']) < 20]
    assert (tropics['pro'] > tropics['syn']).mean() > 0.8
    # Shallowest-bottle rule (Q&A #33): never deeper than any same-station
    # bottle that also passed the depth cut.
    merged = fcm[fcm['depth'] <= 10.0].groupby('station')['depth'].min()
    assert np.allclose(surf.set_index('station')['depth'], merged[surf['station']])


@needs_amt24
def test_process_day_real():
    out = mpipe.process_day(mio.hsas_day_paths()[14])   # a mid-cruise day
    a = out['attrition']
    assert a['raw'] >= a['geometry'] >= a['one_per_minute'] >= a['qc'] > 0
    rrs = out['rrs']
    assert rrs.shape[1] == 124
    assert np.isfinite(rrs).all()
    # QC guarantees no negative Rrs survived into 414-660.
    assert (rrs >= 0).all()
    # Glint parameters within bounds.
    assert (out['records']['rho_sky'] >= 0).all()
    assert (out['records']['rho_sky'] <= mpipe.DEFAULT_PIPELINE['rho_max']).all()


@needs_amt24
def test_brewin2023_loader_untangles_column_families():
    # The CSV interleaves three Rrs families (plain, BRDF-corrected,
    # uncertainty-%) — regression for the prompt-14 bug where all three were
    # swept into one 453-column "spectrum".
    b = mio.load_brewin2023()
    assert b['wave'].shape == (151,)
    assert b['wave'][0] == 400.0 and b['wave'][-1] == 700.0
    assert b['rrs'].shape == (len(b['stations']), 151)
    assert set(b['stations']['cruise'].unique()) == {23, 25, 28}
    # Physical Rrs magnitudes, and uncertainties are percentages not Rrs.
    assert np.nanmedian(b['rrs']) < 0.05
    assert np.nanmedian(b['rrs_unc_pct']) > 0.1
    # Plain vs BRDF-corrected must actually differ.
    b2 = mio.load_brewin2023(brdf=True)
    assert not np.allclose(np.nan_to_num(b['rrs']), np.nan_to_num(b2['rrs']))


@needs_amt24
def test_heldout_retrievals_physical():
    from epft_up.moana.validation import validate_heldout_cruises
    out = validate_heldout_cruises(verbose=False)
    r = out['retrievals']
    assert out['missing_counts'] == [23, 25, 28]     # counts not on disk yet
    assert len(r) == 90
    # Oligotrophic-Atlantic magnitudes; a scrambled spectrum matrix (the
    # prompt-14 bug) drove Pro to ~-5e6 and Syn to 0.
    assert 5e4 < r['pro'].median() < 1e6
    assert 1e3 < r['syn'].median() < 1e5
    assert (r['pro'] < 0).mean() < 0.2


@needs_amt24
def test_heldout_scoring_with_fcm_tables():
    # Target (ii) end-to-end with the AMT23/25/28 truth tables (prompt 15).
    from epft_up.moana.validation import validate_heldout_cruises
    tables = {c: mio.load_fcm(cruise=c) for c in (23, 25, 28)}
    for c, t in tables.items():
        assert {'pro', 'syn', 'peuk'} <= set(t.columns), c
        assert len(t) > 500, c
    out = validate_heldout_cruises(fcm_tables=tables, verbose=False)
    assert out['missing_counts'] == []
    assert set(out['metrics']) == {23, 25, 28, 'pooled'}
    pooled = out['metrics']['pooled']
    # Enough matchups to mean something, and the robust taxon behaves:
    # picoeukaryotes transferred with MAE < 2 in the prompt-15 run.
    assert all(pooled[t]['n'] >= 50 for t in ('pro', 'syn', 'peuk'))
    assert pooled['peuk']['mae'] < 2.0
    assert pooled['peuk']['r2'] > 0.5


@needs_amt24
def test_uway_sst_loads_and_covers_cruise():
    sst = mio.load_uway_sst()
    assert sst['sst'].between(5, 32).all()    # physical Atlantic range
    assert sst['datetime'].is_monotonic_increasing
    assert sst['datetime'].iloc[0].year == 2014
