"""AMT24 HyperSAS Level-2 → MOANA training matrix (design doc §5).

The delivery is calibrated but otherwise raw — fixed-ρ Rrs, no screening, no
binning (prompt-8 report) — so this module re-implements the Lange et al.
(2020) §2.3 processing chain from the LT / Lsky / Ed component streams
(Q&A #29). Stages, each a plain function taking a config dict:

1. :func:`screen_geometry` — tilt / solar-zenith / relative-azimuth cuts;
2. :func:`select_glint_minima` — one spectrum per minute, minimum NIR Lt;
3. :func:`fit_glint` — per-spectrum (ρ_sky, L_NIR) by L1 minimisation over
   750–800 nm, then ``Rrs = (Lt − ρ_sky·Lsky − L_NIR)/Ed``;
4. :func:`qc_spectra` — local-solar-time window, no visible negatives,
   second-derivative noise filter;
5. :func:`resample_to_moana` — 141-band native grid → 414–660 nm @ 2 nm;
6. :func:`build_training_matrix` — median-bin the screened stream onto the
   surface flow-cytometry samples (±15 min, Q&A #32/#33) with the Jordan
   ``uway_sst`` attached (Q&A #34).

Every Lange threshold lives in :data:`DEFAULT_PIPELINE` with a citation, and
:func:`process_day` records per-stage attrition counts so a run's
raw → matched table (the first thing to compare against Lange's n) is a
standard output.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from epft_up.moana.io import (HSAS_WAVE, hsas_day_paths, load_luts,
                              read_hsas_day)

#: Every tunable in the chain. Values are Lange et al. (2020) §2.3 verbatim
#: unless marked as ours; a run's provenance is this dict, serialised.
DEFAULT_PIPELINE = {
    # -- stage 1: geometry (Lange step 3-4) --
    'max_tilt': 5.0,            # deg; discard tilt >= 5
    'min_solar_zenith': 10.0,   # deg; discard theta0 <= 10
    'max_solar_zenith': 80.0,   # deg; discard theta0 >= 80
    'min_delta_azimuth': 50.0,  # deg; discard |dphi| <= 50
    'max_delta_azimuth': 170.0,  # deg; discard |dphi| >= 170
    # -- stage 2: glint (Lange step 5, Eqs. 1-2) --
    'glint_interval_s': 60.0,   # 1-minute selection intervals
    'nir_lo': 750.0,            # nm; NIR window for selection + rho fit
    'nir_hi': 800.0,            # nm (delivery grid tops out at 796)
    'rho_max': 0.1,             # ours: physical upper bound on rho_sky
    # -- stage 3: spectral QC (Lange §2.3 last para) --
    'local_hour_min': 9.0,      # keep 09:00-17:00 local solar time
    'local_hour_max': 17.0,     # ("local" undefined in the paper; we use
                                #  solar time = UTC + lon/15, design §5.4)
    'neg_lo': 400.0,            # nm; no negative Rrs anywhere in 400-700
    'neg_hi': 700.0,
    'd2_max': 2e-4,             # sr^-1 nm^-2; |d2Rrs/dlambda2| cut ...
    'd2_lo': 610.0,             # ... evaluated in 610-660 nm on the 2 nm grid
    'd2_hi': 660.0,
    # -- stage 5: matchup (Q&A #32/#33) --
    'match_window_min': 15.0,   # +/- minutes around each FCM sample; 0 =>
                                # Lange-strict nearest-minute matching
    'fcm_max_depth': 10.0,      # m; shallowest bottle per station (Q&A #33)
    'min_in_bin': 3,            # ours: fewest spectra for a robust median
}


def tilt_from_pitch_roll(pitch, roll):
    """Platform tilt from pitch and roll [deg] (design §5.2).

    ``tilt = arccos(cos(pitch)·cos(roll))`` — exact for small angles and the
    standard combination for a two-axis attitude sensor.

    Parameters
    ----------
    pitch, roll : array — degrees.

    Returns
    -------
    array — tilt in degrees.
    """
    p, r = np.radians(pitch), np.radians(roll)
    return np.degrees(np.arccos(np.clip(np.cos(p) * np.cos(r), -1.0, 1.0)))


def screen_geometry(day, config=None):
    """Stage 1: Lange geometry screens on one day's stream.

    Parameters
    ----------
    day : dict — from :func:`epft_up.moana.io.read_hsas_day`.
    config : dict, optional — :data:`DEFAULT_PIPELINE` overrides.

    Returns
    -------
    (n,) bool — True where the spectrum survives all three screens
    (tilt < 5°, 10° < θ₀ < 80°, 50° < |Δϕ| < 170°).
    """
    cfg = {**DEFAULT_PIPELINE, **(config or {})}
    tilt = tilt_from_pitch_roll(day['pitch'], day['roll'])
    sza = day['solar_zenith']
    # The delivery's delta-azimuth is signed (sensor − sun, −180..180);
    # Lange's screen is on the magnitude of the separation.
    dphi = np.abs(day['delta_azimuth'])
    return ((tilt < cfg['max_tilt'])
            & (sza > cfg['min_solar_zenith']) & (sza < cfg['max_solar_zenith'])
            & (dphi > cfg['min_delta_azimuth'])
            & (dphi < cfg['max_delta_azimuth']))


def select_glint_minima(day, keep, config=None):
    """Stage 2a: keep one spectrum per interval — the darkest in the NIR.

    Lange step 5: partition the day into 1-minute intervals; within each,
    retain only the spectrum whose mean ``Lt(750–800 nm)`` is minimal (the
    least glint-affected sample of that minute).

    Parameters
    ----------
    day : dict — from :func:`read_hsas_day`.
    keep : (n,) bool — stage-1 survivor mask.
    config : dict, optional.

    Returns
    -------
    (m,) int — indices into the day's arrays, time-ordered, one per
    surviving interval.
    """
    cfg = {**DEFAULT_PIPELINE, **(config or {})}
    idx = np.flatnonzero(keep)
    if idx.size == 0:
        return idx
    nir = (day['wave'] >= cfg['nir_lo']) & (day['wave'] <= cfg['nir_hi'])
    nir_mean = np.nanmean(day['lt'][idx][:, nir], axis=1)
    interval = np.floor(day['time'][idx] * 3600.0
                        / cfg['glint_interval_s']).astype(int)
    # Within each interval label, pick the index with minimal NIR Lt.
    order = np.lexsort((nir_mean, interval))
    interval_sorted = interval[order]
    first = np.ones(order.size, dtype=bool)
    first[1:] = interval_sorted[1:] != interval_sorted[:-1]
    return np.sort(idx[order[first]])


def fit_glint(lt, li, wave, config=None):
    """Stage 2b: per-spectrum sky-glint parameters (Lange Eqs. 1–2).

    Solves ``min_{ρ,L_NIR} Σ_{750–800} |Lt − ρ·Lsky − L_NIR|``. For a fixed
    ρ the optimal L1 offset is the median of ``Lt − ρ·Lsky``, so the problem
    reduces to a bounded 1-D minimisation over ρ.

    Parameters
    ----------
    lt, li : (n, w) arrays — total and sky radiance.
    wave : (w,) array — wavelengths [nm].
    config : dict, optional.

    Returns
    -------
    (rho, l_nir) : ((n,), (n,)) float64 — per-spectrum ρ_sky ∈ [0, rho_max]
    and L_NIR (same units as Lt). L_NIR is *not* forced positive: a negative
    offset is diagnostic of over-corrected sky radiance and we want to see it.
    """
    cfg = {**DEFAULT_PIPELINE, **(config or {})}
    nir = (wave >= cfg['nir_lo']) & (wave <= cfg['nir_hi'])
    lt_n, li_n = lt[:, nir], li[:, nir]
    n = lt_n.shape[0]
    rho = np.full(n, np.nan)
    l_nir = np.full(n, np.nan)
    for i in range(n):
        t, s = lt_n[i], li_n[i]
        ok = np.isfinite(t) & np.isfinite(s)
        if ok.sum() < 4:
            continue
        t, s = t[ok], s[ok]
        cost = lambda r: np.abs(t - r * s - np.median(t - r * s)).sum()
        res = minimize_scalar(cost, bounds=(0.0, cfg['rho_max']),
                              method='bounded')
        rho[i] = res.x
        l_nir[i] = np.median(t - res.x * s)
    return rho, l_nir


def resample_to_moana(rrs, wave_in, wave_out=None):
    """Stage 4: linear resampling onto the MOANA grid (design §5.5).

    One shared implementation with the retrieval's interpolation step —
    plain ``np.interp`` per spectrum. The double-interpolation smoothing
    (native ~3.3 nm → delivered 3.5 nm → 2 nm) is accepted per Q&A #31.

    Parameters
    ----------
    rrs : (n, w) array — spectra on the native grid.
    wave_in : (w,) array — native wavelengths [nm].
    wave_out : (124,) array, optional — defaults to the LUT grid.

    Returns
    -------
    (n, 124) float64.
    """
    if wave_out is None:
        wave_out = load_luts()['wave']
    rrs = np.atleast_2d(rrs)
    out = np.empty((rrs.shape[0], wave_out.size))
    for i in range(rrs.shape[0]):
        out[i] = np.interp(wave_out, wave_in, rrs[i])
    return out


def qc_spectra(rrs_native, rrs_moana, wave_native, wave_moana, time_utc, lon,
               config=None):
    """Stage 3: Lange's three post-processing screens (design §5.4).

    Parameters
    ----------
    rrs_native : (n, w) array — Rrs on the native 3.5 nm grid (negatives
        check runs here, over 400–700 nm).
    rrs_moana : (n, 124) array — resampled Rrs (second-derivative filter runs
        here, over 610–660 nm, matching Lange's interpolate-then-QC order).
    wave_native, wave_moana : arrays — the two grids [nm].
    time_utc : (n,) array — decimal hour UTC.
    lon : (n,) array — longitude [deg E].
    config : dict, optional.

    Returns
    -------
    (n,) bool — True where the spectrum passes all three screens.
    """
    cfg = {**DEFAULT_PIPELINE, **(config or {})}
    # 1. Local solar time window ("local" per design §5.4 assumption).
    local = (time_utc + lon / 15.0) % 24.0
    ok = (local >= cfg['local_hour_min']) & (local <= cfg['local_hour_max'])
    # 2. No negative Rrs in the visible, native grid.
    vis = (wave_native >= cfg['neg_lo']) & (wave_native <= cfg['neg_hi'])
    ok &= ~(rrs_native[:, vis] < 0).any(axis=1)
    # 3. Second-derivative noise filter on the 2 nm grid, 610-660 nm.
    step = np.median(np.diff(wave_moana))
    d2 = np.diff(rrs_moana, n=2, axis=1) / step ** 2
    mid = wave_moana[1:-1]  # central-difference wavelengths
    red = (mid >= cfg['d2_lo']) & (mid <= cfg['d2_hi'])
    ok &= ~(np.abs(d2[:, red]) > cfg['d2_max']).any(axis=1)
    # Non-finite spectra fail QC outright.
    ok &= np.isfinite(rrs_moana).all(axis=1)
    return ok


def process_day(sav_path, config=None, wave_moana=None):
    """Stages 0–4 for one delivery day: .sav → screened 1-minute Rrs records.

    Parameters
    ----------
    sav_path : path-like — the day's ``.sav`` file.
    config : dict, optional — :data:`DEFAULT_PIPELINE` overrides.
    wave_moana : (124,) array, optional — target grid; defaults to the LUT's.

    Returns
    -------
    dict with
        ``records`` : pandas.DataFrame — one row per surviving 1-minute
            spectrum: ``doy``, ``time`` (decimal hr UTC), ``datetime``,
            ``lat``, ``lon``, ``rho_sky``, ``l_nir``, ``sst_hsas``
            (cross-check channel), and ``rrs`` holding the (124,) spectrum
            (object column of arrays).
        ``rrs`` : (m, 124) float64 — the same spectra as a matrix.
        ``attrition`` : dict — counts after each stage
            (raw / geometry / one_per_minute / qc).
    """
    cfg = {**DEFAULT_PIPELINE, **(config or {})}
    if wave_moana is None:
        wave_moana = load_luts()['wave']
    day = read_hsas_day(sav_path)

    keep = screen_geometry(day, cfg)                      # stage 1
    idx = select_glint_minima(day, keep, cfg)             # stage 2a
    lt, li, es = day['lt'][idx], day['li'][idx], day['es'][idx]
    rho, l_nir = fit_glint(lt, li, day['wave'], cfg)      # stage 2b
    with np.errstate(invalid='ignore', divide='ignore'):
        rrs_native = (lt - rho[:, None] * li - l_nir[:, None]) / es
    rrs124 = resample_to_moana(rrs_native, day['wave'], wave_moana)
    ok = qc_spectra(rrs_native, rrs124, day['wave'], wave_moana,  # stage 3
                    day['time'][idx], day['lon'][idx], cfg)

    attrition = {'raw': int(day['time'].size),
                 'geometry': int(keep.sum()),
                 'one_per_minute': int(idx.size),
                 'qc': int(ok.sum())}

    sel = idx[ok]
    # Decimal DOY + hour -> real timestamps (year is fixed by the cruise).
    base = pd.Timestamp('2014-01-01')
    dt = (base + pd.to_timedelta(day['doy'] - 1, unit='D')
          + pd.to_timedelta(day['time'][sel], unit='h'))
    records = pd.DataFrame({
        'doy': day['doy'],
        'time': day['time'][sel],
        'datetime': dt,
        'lat': day['lat'][sel], 'lon': day['lon'][sel],
        'rho_sky': rho[ok], 'l_nir': l_nir[ok],
        'sst_hsas': day['sst'][sel],
    })
    rrs_ok = rrs124[ok]
    records['rrs'] = list(rrs_ok)
    return {'records': records, 'rrs': rrs_ok, 'attrition': attrition}


def process_cruise(os_color=None, config=None, verbose=True):
    """Run :func:`process_day` over every available delivery day.

    Parameters
    ----------
    os_color : path-like, optional — data-tree root (default ``$OS_COLOR``).
    config : dict, optional.
    verbose : bool, optional — print one attrition line per day.

    Returns
    -------
    dict with ``records`` (concatenated DataFrame), ``rrs`` ((M, 124)
    matrix), and ``attrition`` (per-day dict of dicts).
    """
    frames, mats, attr = [], [], {}
    for path in hsas_day_paths(os_color):
        out = process_day(path, config)
        attr[out['records']['doy'].iloc[0] if len(out['records']) else
             int(path.stem.rsplit('-', 1)[1])] = out['attrition']
        if verbose:
            a = out['attrition']
            print(f"{path.stem}: raw {a['raw']:6d} -> geom {a['geometry']:6d}"
                  f" -> 1/min {a['one_per_minute']:4d} -> qc {a['qc']:4d}")
        if len(out['records']):
            frames.append(out['records'])
            mats.append(out['rrs'])
    records = pd.concat(frames, ignore_index=True)
    return {'records': records,
            'rrs': np.vstack(mats),
            'attrition': attr}


def _time_ns(times):
    """Datetimes → float64 nanoseconds since epoch, resolution-proof.

    pandas ≥ 2 keeps whatever resolution the parser produced (s/us/ns), and a
    bare ``astype('int64')`` silently returns *that* unit — mixing sources
    then shifts time axes by orders of magnitude. Force ns first.

    Parameters
    ----------
    times : pandas.Series/DatetimeIndex/array of datetimes.

    Returns
    -------
    (n,) float64 — ns since epoch.
    """
    return (pd.to_datetime(times).astype('datetime64[ns]')
            .astype('int64').to_numpy(dtype=np.float64))


def attach_sst(match_times, uway_sst):
    """Interpolate the Jordan ``uway_sst`` series to given timestamps.

    Parameters
    ----------
    match_times : pandas.Series of Timestamps.
    uway_sst : pandas.DataFrame — from
        :func:`epft_up.moana.io.load_uway_sst` (columns datetime, sst).

    Returns
    -------
    (n,) float64 — SST [°C]; NaN outside the series' time span.
    """
    t = _time_ns(uway_sst['datetime'])
    y = uway_sst['sst'].to_numpy()
    q = _time_ns(match_times)
    out = np.interp(q, t, y)
    out[(q < t[0]) | (q > t[-1])] = np.nan
    return out


def build_training_matrix(stream, fcm_surface, uway_sst, config=None):
    """Stage 5: median-bin the screened Rrs stream onto the FCM samples.

    Parameters
    ----------
    stream : dict — from :func:`process_cruise` (needs ``records`` + ``rrs``).
    fcm_surface : pandas.DataFrame — from
        :func:`epft_up.moana.io.surface_fcm` (one row per station).
    uway_sst : pandas.DataFrame — from :func:`load_uway_sst`; the primary SST
        (Q&A #34). The HSAS ancillary SST median rides along as cross-check.
    config : dict, optional — ``match_window_min = 0`` selects Lange-strict
        nearest-minute matching (design §5.6); otherwise the ±window median
        with the 16–84 % spread kept per band.

    Returns
    -------
    dict with
        ``matchups`` : pandas.DataFrame — one row per matched station:
            station, datetime, lat, lon, depth, pro/syn/peuk counts,
            ``n_in_bin``, ``dt_min`` (bin-centre offset of nearest spectrum),
            ``sst`` (primary), ``sst_hsas`` (cross-check), ``sst_diff``.
        ``rrs`` : (k, 124) float64 — the matched (median) spectra.
        ``rrs_spread`` : (k, 124) float64 — half the 16–84 % range per band
            (NaN in Lange-strict mode).
    """
    cfg = {**DEFAULT_PIPELINE, **(config or {})}
    rec, rrs = stream['records'], stream['rrs']
    t_rrs = _time_ns(rec['datetime'])

    rows, spectra, spreads = [], [], []
    for _, sample in fcm_surface.iterrows():
        t0 = float(pd.Timestamp(sample['datetime']).as_unit('ns').value)
        dt_min = (t_rrs - t0) / 60e9  # minutes
        if cfg['match_window_min'] > 0:
            inside = np.abs(dt_min) <= cfg['match_window_min']
            if inside.sum() < cfg['min_in_bin']:
                continue
            block = rrs[inside]
            med = np.median(block, axis=0)
            lo, hi = np.percentile(block, [16, 84], axis=0)
            spread = 0.5 * (hi - lo)
            n_in = int(inside.sum())
            nearest = float(np.abs(dt_min[inside]).min())
            sst_hsas = float(np.median(rec['sst_hsas'].to_numpy()[inside]))
        else:
            # Lange-strict: the single nearest 1-minute spectrum, within 1 min.
            j = int(np.argmin(np.abs(dt_min)))
            nearest = float(abs(dt_min[j]))
            if nearest > 1.0:
                continue
            med = rrs[j]
            spread = np.full(rrs.shape[1], np.nan)
            n_in = 1
            sst_hsas = float(rec['sst_hsas'].iloc[j])
        row = sample.to_dict()
        row.update({'n_in_bin': n_in, 'dt_min': nearest,
                    'sst_hsas': sst_hsas})
        rows.append(row)
        spectra.append(med)
        spreads.append(spread)

    matchups = pd.DataFrame(rows)
    if len(matchups):
        matchups['sst'] = attach_sst(matchups['datetime'], uway_sst)
        matchups['sst_diff'] = matchups['sst'] - matchups['sst_hsas']
    return {'matchups': matchups,
            'rrs': np.array(spectra) if spectra else np.empty((0, 124)),
            'rrs_spread': np.array(spreads) if spreads else np.empty((0, 124))}
