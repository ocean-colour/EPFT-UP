"""AMT24 Level-2 → training matrix, and the AMT23/25/28 transfer set
(cf_design.md §3–4, §5.2, §8.5).

Stages
------
1  geometry screening (tilt, solar zenith, relative azimuth; Lange §2.3)
2  1-minute minimum-NIR selection and L1 glint fit (ρ_sky, L_NIR)
3  spectral QC (solar time, no negative Rrs, second-derivative noise filter
   evaluated on a temporary 2 nm resample)
4  linear resampling to the MOANA grid (QC and baseline only)
5  FCM matchup and ±15 min median binning with σ_Rrs, member spectra, SST,
   qualifier flags, censoring and counting noise

All thresholds and paths live in :data:`DEFAULT_PIPELINE`. Functions take and
return plain dicts / ``xarray.Dataset`` objects.
"""

import copy
import glob
import json
import os
import platform
import subprocess

import numpy as np
import pandas as pd
import xarray as xr

from epft_up.claude_fun import io as cfio
from epft_up.claude_fun.io import TAXA

DEFAULT_PIPELINE = {
    # ---- paths (§2): everything relative to AMT_ROOT = $OS_COLOR/AMT ----------
    'AMT_ROOT': os.path.join(os.environ.get('OS_COLOR', ''), 'AMT'),
    'hsas_glob': 'AMT24/Radiometry/level2/*/AMT24_HSAS_*.sav',
    'afc_files': {'AMT23': 'AMT23/AMT23_JR20131005_AFC_Dataset.csv',
                  'AMT24': 'AMT24/AMT24_JR20140922_AFC_Dataset.csv',
                  'AMT25': 'AMT25/AMT25_JR15001_AFC_Dataset.csv',
                  'AMT28': 'AMT28/AMT28_JR18001_AFC_Dataset.csv'},
    'brewin_file': 'AMT24/BODC_data_AMT_modern_and_historical_optical_observations.csv',
    'jordan_file': 'AMT24/amt24_final_with_debiased_chl.nc',
    'derived_dir': 'AMT24/derived/claude_fun',       # the repo stays data-free
    'train_cruise': 'AMT24',
    # ---- Stage 1: Lange et al. 2020 §2.3, steps 3 (tilt, θ0) and 4 (Δϕ) --------
    'tilt_max_deg': 5.0,
    'sza_range_deg': (10.0, 80.0),
    'dphi_range_deg': (50.0, 170.0),        # screened on |Δϕ|
    # ---- Stage 2: Lange's 1-min minimum-NIR selection and L1 ρ_sky/L_NIR fit ----
    'interval_min': 1.0,
    'glint_window_nm': (750.0, 800.0),
    'rho_sky_bounds': (0.0, 0.1),
    # ---- Stage 3: Lange's post-processing screens ------------------------------
    'solar_time_window_h': (9.0, 17.0),    # "local" time taken as solar time
    'nonneg_range_nm': (400.0, 700.0),
    'd2_threshold': 2e-4,                   # sr^-1 nm^-2 (paper prints nm^-1)
    'd2_range_nm': (610.0, 660.0),
    # ---- Stage 4: MOANA grid (QC, baseline, robustness only) --------------------
    'moana_grid_nm': (414.0, 660.0, 2.0),
    # ---- Stage 5 ---------------------------------------------------------------
    'max_depth_m': 10.0,                    # Q&A #33: shallowest bottle ≤ 10 m
    'sample_offset_min': 45.0,              # AFC time is cast start; +45 min
    'bin_half_window_min': 15.0,            # Q&A #32
    'strict_max_dt_min': 1.0,               # window=0 (Lange-strict): nearest within 1 min
    'afc_columns': {},                      # header overrides for io.read_afc
    'afc_qv_missing': ('N',),
    # one event = DL, V = 1/step (Stage 5 table). AMT25 fits 1.29–1.40: midpoint.
    'count_step': {'AMT23': 1.37, 'AMT24': 1.38, 'AMT25': 1.345, 'AMT28': 0.67},
    'count_step_fit_max': 20.0,
    'sst_var': 'uway_sst',
    # ---- σ floor (§5.2), fitted to Brewin et al. 2023 ---------------------------
    'sigma_f_rel': 0.027,
    'sigma_abs': 8.2e-5,                    # sr^-1
    'nir_diag_nm': (750.0, 796.0),          # σ_abs diagnostic band
    'sigma_abs_tolerance': 2.0,             # stop at checkpoint 1 if off by > 2x
    # ---- transfer test (§8.5) --------------------------------------------------
    'transfer_cruises': ('AMT23', 'AMT25', 'AMT28'),
    'transfer_max_dt_h': 1.5,
    'transfer_max_km': 10.0,
    # ---- SeaFlow + PACE (§8.6) -------------------------------------------------
    'pace_max_dt_h': 3.0,
    'pace_box_half': 2,                     # 5x5 box
    'pace_min_valid_frac': 0.5,
    'seaflow_abundance_factor': 1000.0,     # cells uL^-1 -> mL^-1; check units
}


def config(**overrides):
    """A deep copy of :data:`DEFAULT_PIPELINE` with overrides applied.

    Inputs: keyword overrides. Outputs: dict.
    """
    cfg = copy.deepcopy(DEFAULT_PIPELINE)
    unknown = set(overrides) - set(cfg)
    if unknown:
        raise KeyError(f"unknown pipeline keys {sorted(unknown)}")
    cfg.update(overrides)
    return cfg


def path(cfg, key_or_rel):
    """Absolute path under AMT_ROOT. Inputs: cfg, config key or relative path."""
    rel = cfg.get(key_or_rel, key_or_rel) if isinstance(key_or_rel, str) else key_or_rel
    return os.path.join(cfg['AMT_ROOT'], rel)


# ==========================================================================
# Stage 1 — geometry
# ==========================================================================

def wrap180(a):
    """Wrap angles [deg] to [−180, 180). Inputs/outputs: array."""
    return (np.asarray(a, float) + 180.0) % 360.0 - 180.0


def tilt_deg(pitch, roll):
    """Platform tilt from pitch and roll [deg]: acos(cos p · cos r), which is
    exact for any angle (not a small-angle approximation).

    Inputs: pitch, roll arrays [deg]. Outputs: tilt [deg].
    """
    c = np.cos(np.radians(pitch)) * np.cos(np.radians(roll))
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def geometry_masks(ds, cfg):
    """Stage 1 screens (Lange §2.3): tilt < 5°, 10° < θ0 < 80°,
    50° < |Δϕ| < 170°.

    Inputs: ds (one day from io.read_hsas_day); cfg.
    Outputs: dict of bool arrays ``tilt``, ``sza``, ``dphi`` and the float
    array ``tilt_deg``.
    """
    tilt = tilt_deg(ds['pitch'].values, ds['roll'].values)
    sza = ds['solar_zenith'].values
    dphi = np.abs(wrap180(ds['delta_azimuth'].values))
    lo, hi = cfg['sza_range_deg']
    plo, phi = cfg['dphi_range_deg']
    with np.errstate(invalid='ignore'):
        return {'tilt': tilt < cfg['tilt_max_deg'], 'sza': (sza > lo) & (sza < hi),
                'dphi': (dphi > plo) & (dphi < phi), 'tilt_deg': tilt}


def dphi_convention_check(days, max_per_day=5000, seed=0):
    """Unit check of the file's signed Δϕ against heading − solar azimuth.

    For each heading column (``ship_orientation``, ``ths_compass``) and sign
    s = ±1, computes the circular spread of ``Δϕ_file − s·wrap(heading − ϕ_sun)``
    (``azimuth`` taken as solar azimuth). A small spread means the file's Δϕ
    is a relative azimuth up to a fixed sensor-mounting offset; the screen uses
    |Δϕ|, which is unaffected by the sign convention.

    Inputs: days (list of datasets from io.read_hsas_day; a few clear days).
    Outputs: DataFrame, one row per (column, sign): offset_deg, spread_deg, n;
    sorted best first.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for col in ('ship_orientation', 'ths_compass'):
        for s in (1, -1):
            res = []
            for ds in days:
                d = ds['delta_azimuth'].values
                x = s * wrap180(ds[col].values - ds['azimuth'].values)
                r = np.radians(wrap180(d - x))
                r = r[np.isfinite(r)]
                if r.size > max_per_day:
                    r = rng.choice(r, max_per_day, replace=False)
                res.append(r)
            r = np.concatenate(res) if res else np.array([])
            if r.size == 0:
                continue
            z = np.mean(np.exp(1j * r))
            rows.append({'heading_column': col, 'sign': s,
                         'offset_deg': float(np.degrees(np.angle(z))),
                         'spread_deg': float(np.degrees(np.sqrt(-2 * np.log(max(abs(z), 1e-12))))),
                         'n': int(r.size)})
    return pd.DataFrame(rows).sort_values('spread_deg').reset_index(drop=True)


# ==========================================================================
# Stage 2 — 1-minute selection and glint correction
# ==========================================================================

def one_minute_select(time, lt_nir_mean, interval_min=1.0):
    """Index of the minimum mean-Lt(750–800) spectrum in each interval.

    Inputs: time (datetime64 array); lt_nir_mean (array); interval_min.
    Outputs: int array of selected positions (sorted by time).
    """
    t = pd.to_datetime(np.asarray(time))
    ok = np.isfinite(lt_nir_mean) & t.notna()
    if not ok.any():
        return np.array([], int)
    df = pd.DataFrame({'k': t[ok].floor(f'{interval_min}min'), 'v': np.asarray(lt_nir_mean)[ok],
                       'i': np.flatnonzero(ok)})
    pick = df.loc[df.groupby('k')['v'].idxmin(), 'i'].to_numpy()
    return np.sort(pick)


def fit_glint_l1(lt, li, bounds=(0.0, 0.1), chunk=512):
    """Exact L1 fit ``min Σ_λ |Lt − ρ·Lsky − L_NIR|`` with ρ ∈ bounds,
    L_NIR ≥ 0, per spectrum, over the NIR bands passed in.

    The LP optimum lies at a vertex with two active conditions among: a zero
    residual, a ρ bound, L_NIR = 0. The candidate ρ values are therefore the
    pairwise solutions, Lt/Lsky (L_NIR = 0) and the two bounds; for each the
    optimal L_NIR is max(0, median residual). The best candidate is the
    global optimum (no iteration, no tolerance).

    Inputs: lt, li (n, m) NIR radiances; bounds (ρ_min, ρ_max).
    Outputs: (rho, l_nir, objective), each (n,); NaN where the fit is undefined.
    """
    lt = np.atleast_2d(np.asarray(lt, float))
    li = np.atleast_2d(np.asarray(li, float))
    n, m = lt.shape
    lo, hi = bounds
    ii, jj = np.triu_indices(m, 1)
    rho = np.full(n, np.nan)
    lnir = np.full(n, np.nan)
    obj = np.full(n, np.nan)
    for a in range(0, n, chunk):
        L, S = lt[a:a + chunk], li[a:a + chunk]
        with np.errstate(divide='ignore', invalid='ignore'):
            cand = np.concatenate([(L[:, ii] - L[:, jj]) / (S[:, ii] - S[:, jj]), L / S,
                                   np.full((len(L), 1), lo), np.full((len(L), 1), hi)], axis=1)
        cand[~np.isfinite(cand) | (cand < lo) | (cand > hi)] = np.nan
        r = L[:, None, :] - cand[:, :, None] * S[:, None, :]
        c = np.maximum(np.median(r, axis=2), 0.0)
        f = np.abs(r - c[:, :, None]).sum(axis=2)
        f[~np.isfinite(f)] = np.inf
        k = np.argmin(f, axis=1)
        rows = np.arange(len(L))
        good = np.isfinite(f[rows, k]) & np.isfinite(L).all(1) & np.isfinite(S).all(1)
        rho[a:a + chunk] = np.where(good, cand[rows, k], np.nan)
        lnir[a:a + chunk] = np.where(good, c[rows, k], np.nan)
        obj[a:a + chunk] = np.where(good, f[rows, k], np.nan)
    return rho, lnir, obj


# ==========================================================================
# Stage 3 — spectral QC;  Stage 4 — resampling
# ==========================================================================

def solar_time_h(time, lon):
    """Solar time [h] = UTC + lon/15, mod 24 (the "local" time of Lange's
    09:00–17:00 screen; an assumption recorded in provenance).

    Inputs: time (datetime64 array), lon (deg E). Outputs: hours array.
    """
    t = pd.to_datetime(np.asarray(time))
    h = t.hour + t.minute / 60.0 + t.second / 3600.0
    return (np.asarray(h, float) + np.asarray(lon, float) / 15.0) % 24.0


def moana_grid(cfg=None):
    """414–660 nm @ 2 nm (124 bands). Inputs: cfg. Outputs: array."""
    lo, hi, step = (cfg or DEFAULT_PIPELINE)['moana_grid_nm']
    return np.arange(lo, hi + step / 2, step)


def interp_matrix(wave_in, wave_out):
    """Linear-interpolation weights A with ``rrs_out = rrs_in @ A.T``.

    Inputs: wave_in (n_in,) increasing; wave_out (n_out,).
    Outputs: (A, inside) — A (n_out, n_in); inside (n_out,) bool (outputs
    within the input range; others are NaN after :func:`resample_linear`).
    """
    wave_in = np.asarray(wave_in, float)
    wave_out = np.asarray(wave_out, float)
    A = np.stack([np.interp(wave_out, wave_in, e) for e in np.eye(wave_in.size)], axis=1)
    inside = (wave_out >= wave_in[0] - 1e-9) & (wave_out <= wave_in[-1] + 1e-9)
    return A, inside


def resample_linear(wave_in, rrs, wave_out, sigma=None):
    """Stage 4 resampler: linear interpolation of spectra to a new grid, with
    optional σ propagation ``σ_out = sqrt(Σ_j A_ij² σ_j²)`` (independent bands).

    Inputs: wave_in (n_in,); rrs (n, n_in) or (n_in,); wave_out (n_out,);
    sigma (same shape as rrs), optional.
    Outputs: rrs_out (n, n_out) [and sigma_out if sigma given]; NaN outside
    the input range or where any contributing band is NaN.
    """
    A, inside = interp_matrix(wave_in, wave_out)
    R = np.atleast_2d(np.asarray(rrs, float))
    out = np.full((R.shape[0], A.shape[0]), np.nan)
    nz = np.abs(A) > 0
    for i, r in enumerate(R):
        bad = ~np.isfinite(r)
        rr = np.where(bad, 0.0, r)
        o = rr @ A.T
        o[(nz & bad[None, :]).any(axis=1) | ~inside] = np.nan
        out[i] = o
    if sigma is None:
        return out
    S = np.atleast_2d(np.asarray(sigma, float))
    so = np.sqrt(np.where(np.isfinite(S), S, 0.0) ** 2 @ (A ** 2).T)
    so[~np.isfinite(out)] = np.nan
    return out, so


def second_derivative_max(wave, rrs, cfg):
    """Max |d²Rrs/dλ²| [sr⁻¹ nm⁻²] in the noise-filter range, by central
    differences on a temporary linear resample to the 2 nm MOANA grid (the
    threshold was set on a 2 nm grid and a central second difference scales as
    1/h², so it is evaluated there; only the flag is kept).

    Inputs: wave (n_in,), rrs (n, n_in), cfg. Outputs: (n,) max |d²|.
    """
    g = moana_grid(cfg)
    r = resample_linear(wave, rrs, g)
    h = np.diff(g)[0]
    d2 = (r[:, 2:] - 2 * r[:, 1:-1] + r[:, :-2]) / h ** 2
    centre = g[1:-1]
    lo, hi = cfg['d2_range_nm']
    sel = (centre >= lo - 1e-9) & (centre <= hi + 1e-9)
    with np.errstate(invalid='ignore'):
        return np.max(np.abs(d2[:, sel]), axis=1)


def qc_flags(wave, rrs, time, lon, cfg):
    """Stage 3 flags (True = pass).

    Inputs: wave (n_in,), rrs (n, n_in), time, lon (n,), cfg.
    Outputs: dict ``time``, ``nonneg``, ``d2`` (bool arrays), ``d2_max``,
    ``solar_time_h``.
    """
    st = solar_time_h(time, lon)
    t0, t1 = cfg['solar_time_window_h']
    lo, hi = cfg['nonneg_range_nm']
    band = (wave >= lo) & (wave <= hi)
    rr = np.asarray(rrs, float)[:, band]
    d2 = second_derivative_max(wave, rrs, cfg)
    with np.errstate(invalid='ignore'):
        return {'time': (st >= t0) & (st <= t1),
                'nonneg': np.isfinite(rr).all(1) & (rr >= 0).all(1),
                'd2': np.isfinite(d2) & (d2 <= cfg['d2_threshold']),
                'd2_max': d2, 'solar_time_h': st}


def sigma_floor(rrs, cfg):
    """σ_floor = sqrt((f_rel·Rrs)² + σ_abs²) (§5.2). Inputs: rrs, cfg."""
    return np.sqrt((cfg['sigma_f_rel'] * np.asarray(rrs, float)) ** 2 + cfg['sigma_abs'] ** 2)


# ==========================================================================
# Stages 0–3 over the cruise
# ==========================================================================

def process_day(file, cfg):
    """Stages 0–2 for one day, plus the Stage 3 flags.

    Inputs: file (.sav path), cfg.
    Outputs: (ds, counts) — ds dims (obs, wavelength) with the selected 1-min
    spectra: ``rrs`` (fitted), ``rrs_provider``, ``lt``, ``li``, ``es``,
    ``rho_sky``, ``l_nir``, ``glint_obj``, geometry and ancillary variables and
    Stage 3 flags ``qc_time``, ``qc_nonneg``, ``qc_d2``, ``qc_pass``, ``d2_max``;
    counts (dict) raw / tilt / geometry / glint_selected.
    """
    day = cfio.read_hsas_day(file)
    wave = day['wavelength'].values
    g = geometry_masks(day, cfg)
    t_ok = pd.to_datetime(day['time'].values).notna()
    m_geom = g['tilt'] & g['sza'] & g['dphi'] & t_ok
    counts = {'raw': int(day.sizes['spectrum']), 'tilt': int(g['tilt'].sum()),
              'geometry': int(m_geom.sum())}
    lo, hi = cfg['glint_window_nm']
    nir = (wave >= lo) & (wave <= hi)
    idx = np.flatnonzero(m_geom)
    lt = day['lt'].values
    pick = idx[one_minute_select(day['time'].values[idx], lt[idx][:, nir].mean(axis=1),
                                 cfg['interval_min'])]
    counts['glint_selected'] = int(pick.size)
    sub = day.isel(spectrum=pick)
    LT, LI, ES = sub['lt'].values, sub['li'].values, sub['es'].values
    rho, lnir, obj = fit_glint_l1(LT[:, nir], LI[:, nir], cfg['rho_sky_bounds'])
    with np.errstate(invalid='ignore', divide='ignore'):
        rrs = (LT - rho[:, None] * LI - lnir[:, None]) / ES
    q = qc_flags(wave, rrs, sub['time'].values, sub['lon'].values, cfg)
    sp = ('obs', 'wavelength')
    ds = xr.Dataset({
        'rrs': (sp, rrs), 'rrs_provider': (sp, sub['rrs_provider'].values),
        'lt': (sp, LT.astype(np.float32)), 'li': (sp, LI.astype(np.float32)),
        'es': (sp, ES.astype(np.float32)),
        'rho_sky': ('obs', rho), 'l_nir': ('obs', lnir), 'glint_obj': ('obs', obj),
        'tilt_deg': ('obs', g['tilt_deg'][pick]),
        'qc_time': ('obs', q['time']), 'qc_nonneg': ('obs', q['nonneg']),
        'qc_d2': ('obs', q['d2']), 'd2_max': ('obs', q['d2_max']),
        'solar_time_h': ('obs', q['solar_time_h']),
        'doy': ('obs', np.full(pick.size, day.attrs['doy'])),
    }, coords={'wavelength': wave, 'time': ('obs', sub['time'].values)})
    for v in ('lat', 'lon', 'solar_zenith', 'delta_azimuth', 'azimuth', 'ship_orientation',
              'ths_compass', 'sst_hsas', 'salinity', 'wind_speed', 'pitch', 'roll'):
        ds[v] = ('obs', sub[v].values)
    ds['qc_pass'] = ds['qc_time'] & ds['qc_nonneg'] & ds['qc_d2']
    return ds, counts


def process_cruise(cfg=None, files=None, verbose=True):
    """Stages 0–3 over all AMT24 .sav days.

    Inputs: cfg (default :data:`DEFAULT_PIPELINE`); files (list of .sav paths,
    default the ``hsas_glob`` match); verbose.
    Outputs: stream (xarray.Dataset, dims obs × wavelength) with
    ``attrs['attrition']`` (JSON: raw → tilt → geometry → glint_selected →
    qc_time → qc_nonneg → qc_d2), ``attrs['days']`` (JSON list of DOYs) and
    per-day counts in ``attrs['per_day']``.
    """
    cfg = cfg or DEFAULT_PIPELINE
    files = sorted(files or glob.glob(path(cfg, 'hsas_glob')))
    if not files:
        raise FileNotFoundError(f"no .sav files at {path(cfg, 'hsas_glob')}")
    parts, per_day = [], []
    for f in files:
        ds, c = process_day(f, cfg)
        c['doy'] = cfio.parse_hsas_name(f)[1]
        c['qc_pass'] = int(ds['qc_pass'].sum())
        per_day.append(c)
        parts.append(ds)
        if verbose:
            print(f"  DOY {c['doy']}: raw {c['raw']:6d}  geometry {c['geometry']:6d}  "
                  f"1-min {c['glint_selected']:4d}  QC {c['qc_pass']:4d}")
    stream = xr.concat(parts, dim='obs')
    att = {k: int(sum(c[k] for c in per_day)) for k in ('raw', 'tilt', 'geometry', 'glint_selected')}
    qt = stream['qc_time'].values
    qn = qt & stream['qc_nonneg'].values
    att.update(qc_time=int(qt.sum()), qc_nonneg=int(qn.sum()),
               qc_d2=int((qn & stream['qc_d2'].values).sum()))
    stream.attrs.update(attrition=json.dumps(att), per_day=json.dumps(per_day),
                        days=json.dumps(sorted(int(c['doy']) for c in per_day)),
                        files=json.dumps(files))
    return stream


# ==========================================================================
# Diagnostics for checkpoint 1
# ==========================================================================

def nir_sigma_abs(stream, cfg):
    """Robust spread (1.4826·MAD) of glint-corrected Rrs over 750–796 nm on
    QC-passing spectra, where Lw ≈ 0 (§5.2 diagnostic).

    Inputs: stream, cfg. Outputs: dict ``sigma_abs_nir``, ``ratio`` (to
    cfg['sigma_abs']), ``within_tolerance`` (bool), ``n``.
    """
    lo, hi = cfg['nir_diag_nm']
    w = stream['wavelength'].values
    r = stream['rrs'].values[stream['qc_pass'].values][:, (w >= lo) & (w <= hi)]
    r = r[np.isfinite(r)]
    s = 1.4826 * np.median(np.abs(r - np.median(r))) if r.size else np.nan
    ratio = s / cfg['sigma_abs']
    tol = cfg['sigma_abs_tolerance']
    return {'sigma_abs_nir': float(s), 'ratio': float(ratio),
            'within_tolerance': bool(1 / tol <= ratio <= tol), 'n': int(r.size)}


def provider_crosscheck(stream, cfg):
    """Per-band spread of provider fixed-ρ Rrs minus our fitted Rrs on
    QC-passing spectra, beside the median σ_floor (§5.2 sanity check).

    Inputs: stream, cfg. Outputs: DataFrame indexed by wavelength: median,
    half_p16_p84, robust_sd of the difference, sigma_floor_median.
    """
    ok = stream['qc_pass'].values
    d = stream['rrs_provider'].values[ok] - stream['rrs'].values[ok]
    p16, p50, p84 = np.nanpercentile(d, [16, 50, 84], axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(d - p50), axis=0)
    fl = np.nanmedian(sigma_floor(stream['rrs'].values[ok], cfg), axis=0)
    return pd.DataFrame({'median': p50, 'half_p16_p84': (p84 - p16) / 2, 'robust_sd': mad,
                         'sigma_floor_median': fl},
                        index=pd.Index(stream['wavelength'].values, name='wavelength'))


def rho_bound_fraction(stream, cfg):
    """Fraction of selected spectra with ρ_sky at each bound (pile-up = bad
    intervals). Inputs: stream, cfg. Outputs: dict lower, upper."""
    lo, hi = cfg['rho_sky_bounds']
    r = stream['rho_sky'].values
    r = r[np.isfinite(r)]
    return {'lower': float(np.mean(np.isclose(r, lo))), 'upper': float(np.mean(np.isclose(r, hi)))}


def count_step_table(cfg=None):
    """Re-run the detection-limit fit on all four AFC files (Stage 5 table).

    Inputs: cfg. Outputs: DataFrame per cruise: step, volume_ml, method,
    delta_loglik, second_step, configured step.
    """
    cfg = cfg or DEFAULT_PIPELINE
    rows = []
    for cr, rel in cfg['afc_files'].items():
        df = cfio.read_afc(path(cfg, rel), cr, cfg['afc_columns'])
        vals = np.concatenate([df[t].to_numpy(float) for t in TAXA])
        r = cfio.infer_count_step(vals, cfg['count_step_fit_max'])
        r.update(cruise=cr, configured_step=cfg['count_step'][cr])
        rows.append(r)
    return pd.DataFrame(rows).set_index('cruise')


# ==========================================================================
# Stage 5 — labels and matchup
# ==========================================================================

def surface_samples(afc, cruise, cfg):
    """Shallowest ≤ 10 m bottle per station with bottle flag 0, with labels.

    Inputs: afc (from io.read_afc), cruise (str), cfg.
    Outputs: (df, counts). df columns: station, time_cast, sample_time
    (= cast start + 45 min), lat, lon, depth, solar_time_h, doy, and per taxon
    ``<t>`` (abundance), ``<t>_observed``, ``<t>_censored``, ``<t>_missing``,
    ``<t>_dl`` (= one event), ``<t>_n_events``, ``<t>_y`` (log10; log10 DL if
    censored; NaN if missing), ``<t>_sigma_count`` (0.434/√n; 0 unless
    observed). counts: casts, surface, flagged_only.
    """
    step = cfg['count_step'][cruise]
    ok_depth = np.isfinite(afc['depth'].to_numpy()) & (afc['depth'].to_numpy() <= cfg['max_depth_m'])
    ok_flag = afc['bottle_flag'].to_numpy() == 0
    s = (afc[ok_depth & ok_flag].sort_values(['station', 'depth'], kind='stable')
         .drop_duplicates('station', keep='first').reset_index(drop=True))
    flagged_only = set(afc.loc[ok_depth, 'station']) - set(s['station'])
    st = cfio.afc_label_status(s, cfg['afc_qv_missing'])
    s['sample_time'] = s['time_cast'] + pd.Timedelta(minutes=cfg['sample_offset_min'])
    s['solar_time_h'] = solar_time_h(s['sample_time'].values, s['lon'].values)
    s['doy'] = s['sample_time'].dt.dayofyear
    for t in TAXA:
        v = s[t].to_numpy(float)
        o, c, m = st[t]['observed'], st[t]['censored'], st[t]['missing']
        n = np.where(o, np.floor(v / step + 0.5), 0.0)
        n = np.where(o, np.maximum(n, 1.0), n)
        with np.errstate(divide='ignore', invalid='ignore'):
            s[t + '_y'] = np.where(o, np.log10(v), np.where(c, np.log10(step), np.nan))
            s[t + '_sigma_count'] = np.where(o, np.log10(np.e) / np.sqrt(np.maximum(n, 1)), 0.0)
        s[t + '_observed'], s[t + '_censored'], s[t + '_missing'] = o, c, m
        s[t + '_dl'] = step
        s[t + '_n_events'] = n
    counts = {'casts': int(afc['station'].nunique()), 'surface': int(len(s)),
              'flagged_only': len(flagged_only)}
    return s, counts


def _label_vars(s):
    """(sample, taxon) label arrays from a surface-sample DataFrame."""
    def stack(suf, dtype=float):
        return np.stack([s[t + suf].to_numpy(dtype) for t in TAXA], axis=1)
    st = ('sample', 'taxon')
    return {'abundance': (st, stack('')), 'y_log10': (st, stack('_y')),
            'observed': (st, stack('_observed', bool)), 'censored': (st, stack('_censored', bool)),
            'missing': (st, stack('_missing', bool)), 'dl': (st, stack('_dl')),
            'n_events': (st, stack('_n_events')), 'sigma_count': (st, stack('_sigma_count'))}


def match_amt24(stream, surface, sst_series, cfg, window_min=None):
    """Stage 5 matchup and binning for AMT24.

    For each surface FCM sample (time = cast start + 45 min): the per-band
    median of QC-passing 1-min spectra within ±``window_min`` (default 15),
    σ_bin = (p84 − p16)/2, σ = max(σ_bin, σ_floor), member spectra kept for
    replicate augmentation. ``window_min=0`` is Lange-strict: the nearest
    1-min spectrum within ``strict_max_dt_min``; σ is then the floor alone.

    Inputs: stream (process_cruise), surface (surface_samples df),
    sst_series (io.read_jordan_sst), cfg, window_min.
    Outputs: xarray.Dataset dims (sample, wavelength, member, taxon); see
    module docs and cf_design.md §4. ``attrs['attrition']`` holds the matchup
    lines (surface → daylight → radiometry day → in-window → ≤1 censored).
    """
    window = cfg['bin_half_window_min'] if window_min is None else window_min
    qc_idx = np.flatnonzero(stream['qc_pass'].values)
    qc = stream.isel(obs=qc_idx)
    wave = stream['wavelength'].values
    t_q = pd.to_datetime(qc['time'].values).values.astype('datetime64[ns]')
    days = set(json.loads(stream.attrs.get('days', '[]'))) or set(qc['doy'].values.tolist())
    t0, t1 = cfg['solar_time_window_h']
    att = {'surface': int(len(surface))}
    day_ok = (surface['solar_time_h'] >= t0) & (surface['solar_time_h'] <= t1)
    att['daylight'] = int(day_ok.sum())
    rad_ok = day_ok & surface['doy'].isin(days)
    att['radiometry_day'] = int(rad_ok.sum())
    keep, members = [], []
    for i in np.flatnonzero(rad_ok.to_numpy()):
        dt = (t_q - np.datetime64(surface['sample_time'].iloc[i], 'ns')) / np.timedelta64(1, 'm')
        if window == 0:
            j = np.argmin(np.abs(dt)) if dt.size else -1
            mem = np.array([j]) if j >= 0 and abs(dt[j]) <= cfg['strict_max_dt_min'] else np.array([], int)
        else:
            mem = np.flatnonzero(np.abs(dt) <= window)
        if mem.size:
            keep.append(i)
            members.append((mem, dt[mem]))
    att['in_window'] = len(keep)
    cens = np.stack([surface[t + '_censored'].to_numpy(bool) for t in TAXA], axis=1)
    ok2 = [k for k in range(len(keep)) if cens[keep[k]].sum() <= 1]
    att['max_one_censored'] = len(ok2)
    keep = [keep[k] for k in ok2]
    members = [members[k] for k in ok2]
    att['matched'] = len(keep)
    s = surface.iloc[keep].reset_index(drop=True)
    n, w = len(keep), wave.size
    mmax = max([len(m) for m, _ in members] + [1])
    R = qc['rrs'].values
    med, p16, p84 = (np.full((n, w), np.nan) for _ in range(3))
    mem_rrs = np.full((n, mmax, w), np.nan)
    mem_idx = np.full((n, mmax), -1)
    mem_dt = np.full((n, mmax), np.nan)
    nbin = np.zeros(n, int)
    sst_h = np.full(n, np.nan)
    for k, (m, dt) in enumerate(members):
        x = R[m]
        med[k] = np.median(x, axis=0)
        p16[k], p84[k] = np.percentile(x, [16, 84], axis=0)
        mem_rrs[k, :len(m)] = x
        mem_idx[k, :len(m)] = qc_idx[m]          # positions in the stream
        mem_dt[k, :len(m)] = dt
        nbin[k] = len(m)
        sst_h[k] = np.nanmedian(qc['sst_hsas'].values[m])
    sig_bin = (p84 - p16) / 2.0
    floor = sigma_floor(med, cfg)
    sst, gap = (cfio.interp_series(sst_series, s['sample_time'].values) if len(s)
                else (np.array([]), np.array([])))
    ds = xr.Dataset({
        'rrs': (('sample', 'wavelength'), med), 'rrs_p16': (('sample', 'wavelength'), p16),
        'rrs_p84': (('sample', 'wavelength'), p84),
        'sigma_bin': (('sample', 'wavelength'), sig_bin),
        'sigma_floor': (('sample', 'wavelength'), floor),
        'sigma': (('sample', 'wavelength'), np.fmax(sig_bin, floor)),
        'member_rrs': (('sample', 'member', 'wavelength'), mem_rrs),
        'member_index': (('sample', 'member'), mem_idx),
        'member_dt_min': (('sample', 'member'), mem_dt),
        'n_in_bin': ('sample', nbin),
        'sst': ('sample', sst), 'sst_gap_min': ('sample', gap),
        'sst_hsas': ('sample', sst_h), 'sst_diff': ('sample', sst - sst_h),
        'lat': ('sample', s['lat'].to_numpy(float)), 'lon': ('sample', s['lon'].to_numpy(float)),
        'depth': ('sample', s['depth'].to_numpy(float)),
        'solar_time_h': ('sample', s['solar_time_h'].to_numpy(float)),
        'station': ('sample', s['station'].astype(str).to_numpy()),
        'cruise': ('sample', np.full(n, cfg['train_cruise'])),
        'time': ('sample', s['sample_time'].values.astype('datetime64[ns]')),
        'time_cast': ('sample', s['time_cast'].values.astype('datetime64[ns]')),
        **_label_vars(s),
    }, coords={'wavelength': wave, 'taxon': list(TAXA)})
    ds.attrs.update(attrition=json.dumps(att), window_min=float(window),
                    sst_source='jordan_uway_sst', grid='hsas_native_3.5nm')
    return ds


def build_training_matrix(cfg=None, stream=None, window_min=None, verbose=True):
    """Stages 0–5 for AMT24: the matched training matrix.

    Inputs: cfg; stream (reuse a processed stream); window_min (None → 15;
    0 → Lange-strict); verbose.
    Outputs: (matrix, stream, extras) — extras: dict with surface counts,
    the full attrition table (DataFrame), and the checkpoint-1 diagnostics.
    """
    cfg = cfg or DEFAULT_PIPELINE
    if stream is None:
        stream = process_cruise(cfg, verbose=verbose)
    afc = cfio.read_afc(path(cfg, cfg['afc_files'][cfg['train_cruise']]),
                        cfg['train_cruise'], cfg['afc_columns'])
    surf, scounts = surface_samples(afc, cfg['train_cruise'], cfg)
    sst = cfio.read_jordan_sst(path(cfg, 'jordan_file'), cfg['sst_var'])
    mat = match_amt24(stream, surf, sst, cfg, window_min)
    a1 = json.loads(stream.attrs['attrition'])
    a2 = json.loads(mat.attrs['attrition'])
    table = pd.DataFrame(
        [('raw spectra', a1['raw']), ('tilt < 5°', a1['tilt']),
         ('θ0 and |Δϕ| screens', a1['geometry']), ('1-min glint selection', a1['glint_selected']),
         ('solar time 09–17', a1['qc_time']), ('no negative Rrs 400–700', a1['qc_nonneg']),
         ('second-derivative filter', a1['qc_d2']),
         ('FCM casts in deposit', scounts['casts']),
         ('surface samples (≤10 m, flag 0)', a2['surface']),
         ('  no daylight radiometry (dropped)', a2['surface'] - a2['daylight']),
         ('daylight samples', a2['daylight']),
         ('  no radiometry that day (dropped)', a2['daylight'] - a2['radiometry_day']),
         ('on radiometry days', a2['radiometry_day']),
         ('≥1 QC spectrum in window', a2['in_window']),
         ('  ≥2 censored taxa (dropped)', a2['in_window'] - a2['max_one_censored']),
         ('matched labels', a2['matched'])], columns=['stage', 'n'])
    extras = {'surface_counts': scounts, 'attrition': table,
              'nir_sigma_abs': nir_sigma_abs(stream, cfg),
              'rho_bounds': rho_bound_fraction(stream, cfg)}
    return mat, stream, extras


def build_transfer_set(cfg=None):
    """AMT23/25/28 transfer set from the Brewin compilation and AFC files (§8.5).

    Brewin rows with full spectra → Stage 3 QC (no negative Rrs 400–700;
    second-derivative filter on the native 2 nm grid; no time screen, the
    stations are near solar noon) → matched to the shallowest ≤ 10 m AFC
    sample (cast + 45 min) within ±1.5 h and 10 km, nearest in time.
    σ = max(σ_file, σ_floor); σ_floor alone where the file has none.

    Inputs: cfg.
    Outputs: (ds, attrition) — ds in the training-matrix layout on the Brewin
    grid (plus dt_h, dist_km, shared_cast, ctd_method and qc:* columns);
    attrition DataFrame per cruise.
    """
    cfg = cfg or DEFAULT_PIPELINE
    B = cfio.read_brewin(path(cfg, 'brewin_file'))
    wave = B['wave']
    meta = B['meta']
    parts, att = [], []
    for cr in cfg['transfer_cruises']:
        rows = np.flatnonzero((meta['cruise'] == cr).to_numpy() & B['usable'])
        q = qc_flags(wave, B['rrs'][rows], meta['time'].values[rows], meta['lon'].values[rows], cfg)
        qc_ok = q['nonneg'] & q['d2']
        afc = cfio.read_afc(path(cfg, cfg['afc_files'][cr]), cr, cfg['afc_columns'])
        surf, _ = surface_samples(afc, cr, cfg)
        match = np.full(rows.size, -1)
        dt_h = np.full(rows.size, np.nan)
        dist = np.full(rows.size, np.nan)
        for k, r in enumerate(rows):
            d = (surf['sample_time'] - meta['time'].iloc[r]).dt.total_seconds().to_numpy() / 3600
            km = cfio.haversine_km(meta['lat'].iloc[r], meta['lon'].iloc[r],
                                   surf['lat'].to_numpy(), surf['lon'].to_numpy())
            ok = (np.abs(d) <= cfg['transfer_max_dt_h']) & (km <= cfg['transfer_max_km'])
            if ok.any():
                j = np.flatnonzero(ok)[np.argmin(np.abs(d[ok]))]
                match[k], dt_h[k], dist[k] = j, d[j], km[j]
        cens = np.stack([surf[t + '_censored'].to_numpy(bool) for t in TAXA], 1)
        two_c = np.array([m >= 0 and cens[m].sum() > 1 for m in match])
        final = qc_ok & (match >= 0) & ~two_c
        att.append({'cruise': cr, 'usable_rows': int(rows.size), 'matched_before_qc': int((match >= 0).sum()),
                    'fail_nonneg': int((~q['nonneg']).sum()), 'fail_d2': int((q['nonneg'] & ~q['d2']).sum()),
                    'qc_pass': int(qc_ok.sum()), 'two_censored': int(two_c.sum()), 'matched': int(final.sum())})
        if not final.any():
            continue
        r_sel, m_sel = rows[final], match[final]
        s = surf.iloc[m_sel].reset_index(drop=True)
        R = B['rrs'][r_sel]
        fl = sigma_floor(R, cfg)
        sf = B['sigma_file'][r_sel]
        sig = np.where(np.isfinite(sf), np.fmax(sf, fl), fl)
        ctd = meta['ctd_method'].values[r_sel] if 'ctd_method' in meta else np.full(len(r_sel), np.nan)
        ds = xr.Dataset({
            'rrs': (('sample', 'wavelength'), R), 'sigma_file': (('sample', 'wavelength'), sf),
            'sigma_floor': (('sample', 'wavelength'), fl), 'sigma': (('sample', 'wavelength'), sig),
            'n_in_bin': ('sample', np.ones(len(r_sel), int)),
            'sst': ('sample', meta['sst'].values[r_sel] if 'sst' in meta else np.full(len(r_sel), np.nan)),
            'ctd_method': ('sample', ctd.astype(float)),
            'lat': ('sample', s['lat'].to_numpy(float)), 'lon': ('sample', s['lon'].to_numpy(float)),
            'optics_lat': ('sample', meta['lat'].values[r_sel]),
            'optics_lon': ('sample', meta['lon'].values[r_sel]),
            'depth': ('sample', s['depth'].to_numpy(float)),
            'station': ('sample', s['station'].astype(str).to_numpy()),
            'cruise': ('sample', np.full(len(r_sel), cr)),
            'time': ('sample', meta['time'].values[r_sel].astype('datetime64[ns]')),
            'time_cast': ('sample', s['time_cast'].values.astype('datetime64[ns]')),
            'dt_h': ('sample', dt_h[final]), 'dist_km': ('sample', dist[final]),
            'shared_cast': ('sample', pd.Series(m_sel).duplicated(keep=False).to_numpy()),
            **_label_vars(s),
        }, coords={'wavelength': wave, 'taxon': list(TAXA)})
        for c in meta.columns:
            if c.startswith('qc:'):
                ds[c.replace(':', '_').replace(' ', '_')] = ('sample', meta[c].values[r_sel])
        parts.append(ds)
    out = xr.concat(parts, dim='sample') if parts else None
    att = pd.DataFrame(att).set_index('cruise')
    att.loc['total'] = att.sum()
    if out is not None:
        out.attrs.update(sst_source='brewin_temperature_above_secchi', grid='brewin_2nm',
                         attrition=att.to_json())
    return out, att


# ==========================================================================
# Products and provenance
# ==========================================================================

def git_hash():
    """Commit of the package's repository, or 'unknown'. Outputs: str."""
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                       cwd=os.path.dirname(os.path.abspath(__file__)),
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return 'unknown'


def provenance(pipeline_cfg=None, model_cfg=None, train_cfg=None, files=(), seeds=(),
               hash_files=True, **extra):
    """Serialisable provenance: the three configs, git hash, input-file hashes,
    seeds and library versions.

    Inputs: configs (dicts), files (paths), seeds, hash_files, extra items.
    Outputs: dict.
    """
    import scipy
    out = {'pipeline': pipeline_cfg, 'model': model_cfg, 'train': train_cfg,
           'git': git_hash(), 'seeds': list(seeds),
           'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'pandas': pd.__version__,
                        'xarray': xr.__version__},
           'assumptions': {'local_time': 'solar time = UTC + lon/15',
                           'dphi': 'screened on |delta_azimuth| from the file'}}
    try:
        import tensorflow as tf
        out['versions']['tensorflow'] = tf.__version__
    except Exception:
        pass
    out['files'] = {f: (cfio.sha256_file(f) if hash_files and os.path.exists(f) else None)
                    for f in files}
    out.update(extra)
    return out


def _encodable(ds):
    """Copy of ds with bool → int8 and object strings → str (netCDF-safe)."""
    ds = ds.copy()
    for v in list(ds.variables):
        if ds[v].dtype == bool:
            ds[v] = ds[v].astype(np.int8)
            ds[v].attrs['bool'] = 1
        elif ds[v].dtype == object:
            ds[v] = ds[v].astype(str)
    return ds


def save_dataset(ds, file):
    os.makedirs(os.path.dirname(os.path.abspath(file)), exist_ok=True)
    # Replace forward slashes in variable names with underscores to satisfy netcdf
    rename_dict = {v: v.replace('/', '_') for v in ds.data_vars if '/' in v}
    if rename_dict:
        ds = ds.rename(rename_dict)
    _encodable(ds).to_netcdf(file)
    return file


def load_dataset(file):
    """Read a dataset written by :func:`save_dataset`, restoring bools."""
    ds = xr.load_dataset(file)
    for v in list(ds.variables):
        if ds[v].attrs.get('bool') == 1:
            ds[v] = ds[v].astype(bool)
    return ds
