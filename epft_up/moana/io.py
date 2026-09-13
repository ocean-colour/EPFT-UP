"""Readers for everything the MOANA track consumes (design doc §5.1, §2).

Four data sources, four loaders, no processing here:

- :func:`load_luts` — the two vendored NASA OCSSW files in
  ``epft_up/data/moana/`` (PCA loadings + regression coefficients). These are
  the complete set of algorithm constants; see the README next to them.
- :func:`read_hsas_day` — one day of the AMT24 HyperSAS Level-2 delivery.
  **Reads the IDL ``.sav`` only.** The sibling ``*_ES.dat`` CSV is a
  byte-identical copy of ``*_LT.dat`` on all 37 days (Q&A #28), so the CSVs
  must never be used; the ``.sav``'s ``matrix_es`` is correct.
- :func:`load_fcm` — the BODC flow-cytometry bottle file (training truth).
- :func:`load_uway_sst` — the QC'd 1-minute underway SST from the Jordan
  et al. (2025) netCDF, the primary SST source (Q&A #34).

All loaders return plain dicts / :class:`pandas.DataFrame`; the wavelength
axis is always explicit so no downstream code indexes raw ``.sav`` columns.
"""

from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.io import readsav

#: Directory holding the vendored NASA LUTs (sha256-pinned; see its README).
LUT_DIR = Path(__file__).resolve().parent.parent / 'data' / 'moana'

# --- AMT24 HyperSAS Level-2 .sav column layout (prompt-8 inspection) --------
#
# Every matrix is (n_spectra, 158): 14 leading ancillary columns, 141
# spectral columns on 306.0..796.0 nm at exactly 3.5 nm, then 3 trailing
# ancillary columns. The CSV headers carry the same order; we verified the
# .sav matches the CSVs column-for-column (RRS/LT/LI) on real days.
_ANC_LEAD = [
    'ship_start_time',   # decimal hour, ship datalogger start
    'time',              # decimal hour, UTC
    'lat', 'lon',
    'pressure_atmos',    # hPa
    'par',
    'wind_dir',          # deg
    'wind_speed',        # m/s
    'ship_orientation',  # deg
    'salinity',          # ppt (corrupt on DOY 266)
    'sst',               # degC (corrupt on DOY 266; cross-check only, Q&A #34)
    'delta_azimuth',     # deg, signed (sensor - sun); screened on |value|
    'azimuth',           # deg, sensor azimuth
    'solar_zenith',      # deg
]
_ANC_TRAIL = ['ths_compass', 'pitch', 'roll']  # deg
_N_LEAD = len(_ANC_LEAD)           # 14
_N_WAVE = 141
#: Native wavelength grid of the Level-2 delivery [nm].
HSAS_WAVE = 306.0 + 3.5 * np.arange(_N_WAVE)

# --- BODC flow-cytometry parameter codes -------------------------------------
#
# Mapping transcribed from the deposit's own metadata document
# (`AMT24_AFC_Document_468667.htm`, BODC document 468667) — NOT guessed from
# the code strings, which are misleading (P700* is Synechococcus, P701* is
# Prochlorococcus). load_fcm() asserts all expected codes are present.
#: BODC parameter code -> our short taxon name. Derived from each cruise's own
#: metadata document, and **confirmed directly by PML** (Tom, 2026-09-10: "the
#: cell counts are the final set of columns ... P700A90Z is Abundance of
#: Synechococcus"). Verified cell-for-cell against the 2026-09 AMT24 re-delivery.
FCM_CODES = {
    'P701A90Z': 'pro',     # Prochlorococcus       [cells/ml]
    'P700A90Z': 'syn',     # Synechococcus         [cells/ml]
    'PYEUA00A': 'peuk',    # picoeukaryotes        [cells/ml]
    'J79A0596': 'crypto',  # Cryptophyceae         [cells/ml] (not used by MOANA)
    'P490A00Z': 'cocco',   # coccolithophores      [cells/ml] (not used by MOANA)
    'X726A86B': 'nanoeuk',  # nanoeukaryotes 2-12um [cells/ml] (not used by MOANA)
}


def _os_color():
    """Return the ``$OS_COLOR`` data-tree root as a Path.

    Returns
    -------
    Path
    Raises
    ------
    RuntimeError : if the environment variable is not set (Tier-2 data needed).
    """
    root = os.environ.get('OS_COLOR')
    if root is None:
        raise RuntimeError(
            "MOANA data loaders need the $OS_COLOR data tree; "
            "the variable is not set.")
    return Path(root)


def _amt24_dir():
    """The AMT24 data directory, wherever the tree currently keeps it.

    The layout changed once already (2026-08-16: ``$OS_COLOR/AMT24/`` →
    ``$OS_COLOR/AMT/AMT24/``), so every loader resolves through here.

    Returns
    -------
    Path
    """
    root = _os_color()
    for cand in (root / 'AMT' / 'AMT24', root / 'AMT24'):
        if cand.is_dir():
            return cand
    raise FileNotFoundError(
        f'AMT24 data not found under {root}/AMT/AMT24 or {root}/AMT24')


def load_luts(lut_dir=None):
    """Load the vendored NASA MOANA constants.

    Parameters
    ----------
    lut_dir : path-like, optional
        Directory holding ``pca_picophyto.h5`` and ``picophyt.json``.
        Defaults to the vendored copies in ``epft_up/data/moana/``.

    Returns
    -------
    dict with keys
        ``wave`` : (124,) float64 — 414..660 nm at 2 nm (from the HDF5).
        ``V`` : (124, 45) float64 — PCA loading matrix, columns orthonormal.
        ``npc`` : int — number of scores the regressions consume (17).
        ``pro_coef`` : (npc+2,) float64 — [intercept, log10(SST), U1..U17].
        ``syn_coef`` : (npc+1,) float64 — [intercept, U1..U17].
        ``apeuk_coef`` : (npc+1,) float64 — [intercept, U1..U17].
    """
    lut_dir = Path(lut_dir) if lut_dir is not None else LUT_DIR
    with h5py.File(lut_dir / 'pca_picophyto.h5', 'r') as f:
        V = f['component'][()].astype(np.float64)
        wave = f['wavelength'][()].astype(np.float64)
    # The HDF5 may store components as (45, 124) or (124, 45) depending on
    # the writer's row/column convention; normalise to (n_wave, n_pc).
    if V.shape[0] != wave.size:
        V = V.T
    with open(lut_dir / 'picophyt.json') as f:
        raw = json.load(f)
    npc = int(raw['npc'])
    parse = lambda s: np.array([float(x) for x in s.split(',')])
    return {
        'wave': wave,
        'V': V,
        'npc': npc,
        'pro_coef': parse(raw['pro_coef']),
        'syn_coef': parse(raw['syn_coef']),
        'apeuk_coef': parse(raw['apeuk_coef']),
    }


def hsas_day_paths(os_color=None):
    """List the available AMT24 HyperSAS ``.sav`` files, sorted by day.

    Parameters
    ----------
    os_color : path-like, optional
        Data-tree root; defaults to ``$OS_COLOR``.

    Returns
    -------
    list of Path — one ``.sav`` per available day (37 for the 2026-08 delivery;
    DOY 267, 273, 298 are absent upstream, Q&A #30).
    """
    if os_color is not None:
        base = Path(os_color) / 'AMT24'
    else:
        base = _amt24_dir()
    pat = str(base / 'Radiometry' / 'level2' / '*' / '*.sav')
    return sorted(Path(p) for p in glob.glob(pat))


def read_hsas_day(sav_path):
    """Read one day of the AMT24 HyperSAS Level-2 delivery from its ``.sav``.

    Parameters
    ----------
    sav_path : path-like
        e.g. ``.../level2/2014280/AMT24_HSAS_2014-280.sav``.

    Returns
    -------
    dict with keys
        ``doy`` : int — day of year (from the filename).
        ``wave`` : (141,) float64 — native grid, 306..796 nm @ 3.5 nm.
        ``lt``, ``li``, ``es``, ``rrs_provider`` : (n, 141) float64 —
            total radiance, sky radiance, downwelling irradiance
            [µW cm⁻² nm⁻¹ (sr⁻¹)], and the provider's fixed-ρ Rrs [sr⁻¹]
            (cross-check channel only, Q&A #29).
        one (n,) float64 array per ancillary name in the module-level
        ``_ANC_LEAD`` / ``_ANC_TRAIL`` lists (``time`` is decimal hour UTC).

    Notes
    -----
    ES comes from ``matrix_es`` — correct in the ``.sav``, broken in the CSV
    export (Q&A #28). A guard raises if ES and LT are identical, so the known
    bug can never silently re-enter through a "fixed" future delivery that
    breaks the ``.sav`` instead.
    """
    sav_path = Path(sav_path)
    s = readsav(str(sav_path))
    lt = np.asarray(s['matrix_lt'], dtype=np.float64)
    li = np.asarray(s['matrix_li'], dtype=np.float64)
    es = np.asarray(s['matrix_es'], dtype=np.float64)
    rrs = np.asarray(s['matrix_rrs'], dtype=np.float64)
    if np.array_equal(np.nan_to_num(es), np.nan_to_num(lt)):
        raise ValueError(
            f"{sav_path.name}: matrix_es is identical to matrix_lt — "
            "the ES-duplication bug (Q&A #28) has reached the .sav files.")
    wsl = slice(_N_LEAD, _N_LEAD + _N_WAVE)
    out = {
        'doy': int(sav_path.stem.rsplit('-', 1)[1]),
        'wave': HSAS_WAVE.copy(),
        'lt': lt[:, wsl], 'li': li[:, wsl], 'es': es[:, wsl],
        'rrs_provider': rrs[:, wsl],
    }
    # Ancillaries are identical across the four matrices; take them from RRS.
    for j, name in enumerate(_ANC_LEAD):
        out[name] = rrs[:, j]
    for j, name in enumerate(_ANC_TRAIL):
        out[name] = rrs[:, _N_LEAD + _N_WAVE + j]
    return out


def load_fcm(csv_path=None, cruise=24):
    """Load a BODC AMT flow-cytometry bottle dataset (cell-count truth).

    All four deposits on disk (AMT23/24/25/28, Tarran & Zubkov 2020 / Tarran
    2020) share one layout and one BODC parameter vocabulary — the
    code-to-taxon mapping was verified against each cruise's own metadata
    document (P700A90Z = *Synechococcus*, P701A90Z = *Prochlorococcus*).

    Parameters
    ----------
    csv_path : path-like, optional
        Explicit dataset path; overrides ``cruise``.
    cruise : int, optional
        AMT cruise number (23, 24, 25 or 28); resolves
        ``$OS_COLOR/AMT/AMT<cruise>/AMT<cruise>_*_AFC_Dataset.csv``, falling
        back to the legacy ``$OS_COLOR/AMT24/`` location for cruise 24.

    Returns
    -------
    pandas.DataFrame — one row per bottle, columns:
        ``station`` (str), ``datetime`` (UTC), ``lat``, ``lon``,
        ``depth`` [m], and one cells mL⁻¹ column per taxon short name in
        :data:`FCM_CODES` (``pro``, ``syn``, ``peuk``, ...). Missing counts
        are NaN.

    Notes
    -----
    These deposits are CTD-bottle-only; the underway FCM samples Lange+2020
    also trained on are in none of them (Q&A #35 — PML follow-up list).
    """
    if csv_path is None:
        hits = sorted(glob.glob(str(
            _os_color() / 'AMT' / f'AMT{cruise}' / f'AMT{cruise}_*_AFC_Dataset.csv')))
        if not hits and cruise == 24:
            hits = sorted(glob.glob(str(
                _amt24_dir() / 'AMT24_*_AFC_Dataset.csv')))
        if not hits:
            raise FileNotFoundError(
                f'no AFC dataset for AMT{cruise} under $OS_COLOR/AMT/')
        csv_path = hits[0]
    df = pd.read_csv(csv_path)
    missing = [c for c in FCM_CODES if f'{c}[#/ml]' not in df.columns]
    if missing:
        raise ValueError(
            f"{csv_path}: expected BODC parameter codes absent: {missing} — "
            "re-check the code mapping against the deposit's metadata "
            "document before proceeding.")
    out = pd.DataFrame({
        'station': df['Orig_stn'],
        # BODC writes dd/mm/yyyy hh:mm in the timestamp-labelled column.
        'datetime': pd.to_datetime(
            df['yyyy-mm-ddThh24:mi:ss[GMT]'], format='%d/%m/%Y %H:%M'),
        'lat': df['Latitude[deg+veN]'],
        'lon': df['Longitude[deg+veE]'],
        'depth': df['Bot_depth[metres]'],
    })
    for code, short in FCM_CODES.items():
        out[short] = pd.to_numeric(df[f'{code}[#/ml]'], errors='coerce')
    return out


def surface_fcm(fcm, max_depth=10.0):
    """Reduce bottle samples to one surface sample per station (Q&A #33).

    Keeps bottles with ``depth <= max_depth`` and, per station, the single
    *shallowest* one — not a 0–10 m mean, to avoid mixing across the
    near-surface Prochlorococcus gradient (Lange et al. 2018).

    Parameters
    ----------
    fcm : pandas.DataFrame — from :func:`load_fcm`.
    max_depth : float, optional — depth cut [m]; confirmed default 10 m.

    Returns
    -------
    pandas.DataFrame — one row per station, same columns as the input.
    """
    shallow = fcm[fcm['depth'] <= max_depth]
    idx = shallow.groupby('station')['depth'].idxmin()
    return shallow.loc[idx].sort_values('datetime').reset_index(drop=True)


def load_brewin2023(csv_path=None, cruises=(23, 25, 28), brdf=False):
    """Load the Brewin et al. (2023) AMT station optical dataset.

    In-situ hyperspectral Rrs (BRDF-corrected after Lee et al. 2011, 2 nm
    from 400 nm) plus temperature and ancillaries at 127 stations on AMT
    23/25/26/28 — the held-out-cruise Rrs for validation target (ii).
    BODC DOI 10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6.

    Parameters
    ----------
    csv_path : path-like, optional
        Defaults to ``$OS_COLOR/AMT24/BODC_data_AMT_modern_and_historical_
        optical_observations.csv`` (where the deposit was unzipped).
    cruises : tuple of int, optional
        AMT cruise numbers to keep. Default (23, 25, 28) — the three that are
        Lange+2020 held-out validation cruises. AMT26 is available but has no
        Lange comparison line.
    brdf : bool, optional
        The file carries **three** interleaved 151-column Rrs families on the
        same 400–700 @ 2 nm grid: plain ``Rrs(λ) [1/sr]``, BRDF-corrected
        ``Rrs(λ) (BDRF-corrected Lee et al. 2011) [1/sr]``, and
        ``Rrs(λ) Uncertainty [%]``. Default False selects the *plain* family
        — the above-water geometry MOANA's training Rrs had; True selects the
        BRDF-corrected one (closer to satellite normalised Rrs). The percent
        uncertainties return alongside either way.

    Returns
    -------
    dict with
        ``stations`` : pandas.DataFrame — cruise, datetime (UTC), lat, lon,
            ``sst`` [°C] (temperature above Secchi depth — the in-situ
            near-surface temperature), chl, solar zenith angle.
        ``wave`` : (151,) float64 — Rrs wavelengths [nm], 400–700 @ 2 nm.
        ``rrs`` : (n, 151) float64 — Rrs [sr⁻¹]; BODC −999 fills as NaN.
        ``rrs_unc_pct`` : (n, 151) float64 — per-band uncertainty [%].
    """
    if csv_path is None:
        csv_path = (_amt24_dir() /
                    'BODC_data_AMT_modern_and_historical_optical_observations.csv')
    df = pd.read_csv(csv_path)
    df = df[df['AMT Cruise'].isin(cruises)].reset_index(drop=True)

    def family(predicate):
        # One (wavelength-sorted) column family; -999 is the BODC fill.
        cols = [c for c in df.columns if c.startswith('Rrs(') and predicate(c)]
        w = np.array([float(c.split('(')[1].split(')')[0]) for c in cols])
        order = np.argsort(w)
        m = df[cols].to_numpy(dtype=np.float64)[:, order]
        m[m <= -998] = np.nan
        return w[order], m

    want = ('BDRF' if brdf else None)
    wave, rrs = family(lambda c: ('Uncertainty' not in c) and
                       (('BDRF' in c) == (want == 'BDRF')))
    wave_u, rrs_unc = family(lambda c: 'Uncertainty' in c)
    assert np.array_equal(wave, wave_u), 'Rrs and uncertainty grids differ'
    stations = pd.DataFrame({
        'cruise': df['AMT Cruise'],
        'datetime': pd.to_datetime(
            df['Year-Month-Day'] + ' ' + df['Time [GMT]']),
        'lat': df['Latitude [deg+veN]'],
        'lon': df['Longitude [deg+veE]'],
        'sst': pd.to_numeric(
            df['Temperature above Secchi depth [degC]'], errors='coerce'
        ).replace(-999.0, np.nan),
        'chl': pd.to_numeric(
            df['Total chlorophyll-a [mg/m^3]'], errors='coerce'
        ).replace(-999.0, np.nan),
        'solar_zenith': pd.to_numeric(
            df['Solar Zenith Angle [deg]'], errors='coerce'),
    })
    return {'stations': stations, 'wave': wave, 'rrs': rrs,
            'rrs_unc_pct': rrs_unc}


def load_uway_sst(nc_path=None):
    """Load the QC'd 1-minute underway SST (primary SST source, Q&A #34).

    Parameters
    ----------
    nc_path : path-like, optional
        Defaults to ``$OS_COLOR/AMT24/amt24_final_with_debiased_chl.nc``
        (Jordan et al. 2025, DOI 10.5281/zenodo.12527954).

    Returns
    -------
    pandas.DataFrame — columns ``datetime`` (UTC) and ``sst`` [°C], NaNs
    dropped, sorted by time.
    """
    import xarray as xr  # local import: only this loader needs it
    if nc_path is None:
        nc_path = _amt24_dir() / 'amt24_final_with_debiased_chl.nc'
    ds = xr.open_dataset(nc_path)
    df = pd.DataFrame({
        'datetime': pd.to_datetime(ds['time'].values),
        'sst': np.asarray(ds['uway_sst'].values, dtype=np.float64),
    })
    ds.close()
    return df.dropna(subset=['sst']).sort_values('datetime').reset_index(drop=True)
