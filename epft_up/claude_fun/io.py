"""Readers for Claude Fun (cf_design.md §2, §3 Stage 0 and Stage 5, §5.2, §8.5, §8.6).

Every reader returns plain containers (``xarray.Dataset``, ``pandas.DataFrame``
or dict) with explicit wavelength axes and named columns, so no downstream code
indexes raw file columns.

Readers
-------
read_hsas_day      AMT24 HyperSAS Level-2 IDL save (.sav only; the CSVs are
                   never read: ``*_ES.dat`` is a copy of ``*_LT.dat``, Q&A #28).
read_afc           BODC flow-cytometry CSV (AMT23/24/25/28), header-mapped.
infer_count_step   detection limit (one event) from count quantisation.
read_jordan_sst    ``uway_sst`` from the Jordan et al. 2025 netCDF.
read_brewin        Brewin et al. 2023 optical compilation (uncorrected Rrs).
fit_sigma_floor    two-term σ floor fit to the Brewin uncertainties (§5.2).
read_seaflow       SeaFlow v2 abundance table (§8.6, test only).
read_pace_box      5×5 PACE OCI L2 box around a position (§8.6, test only).
"""

import csv
import hashlib
import os
import re
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import readsav
from scipy.optimize import least_squares

#: Taxon order used everywhere in Claude Fun.
TAXA = ('pro', 'syn', 'peuk')

# --------------------------------------------------------------------------
# HyperSAS Level-2 (.sav)
# --------------------------------------------------------------------------

#: 141 native bands, 306.0–796.0 nm at exactly 3.5 nm (columns 14–154).
HSAS_WAVE = 306.0 + 3.5 * np.arange(141)
#: Columns 0–13 of every matrix (identical in all four matrices).
HSAS_ANCILLARY = ('ship_start_time', 'time_utc_h', 'lat', 'lon', 'pressure',
                  'par', 'wind_dir', 'wind_speed', 'ship_orientation',
                  'salinity', 'sst_hsas', 'delta_azimuth', 'azimuth',
                  'solar_zenith')
#: Columns 155–157.
HSAS_TAIL = ('ths_compass', 'pitch', 'roll')
#: Output variable → IDL save variable.
HSAS_MATRICES = {'lt': 'matrix_lt', 'li': 'matrix_li', 'es': 'matrix_es',
                 'rrs_provider': 'matrix_rrs'}
_N_COLS = 158


def parse_hsas_name(path):
    """Year and day of year from ``AMT24_HSAS_<yyyy-ddd>.sav``.

    Inputs: path (str). Outputs: (year, doy) ints.
    """
    m = re.search(r'HSAS_(\d{4})-(\d{3})\.sav$', os.path.basename(path))
    if m is None:
        raise ValueError(f"not an HSAS Level-2 save name: {path}")
    return int(m.group(1)), int(m.group(2))


def read_hsas_day(path):
    """Read one per-day HyperSAS Level-2 IDL save.

    Inputs
    ------
    path : str — ``.../AMT24_HSAS_<yyyy-ddd>.sav``.

    Outputs
    -------
    xarray.Dataset, dims (spectrum, wavelength):
      ``lt``, ``li``, ``es``, ``rrs_provider`` (spectrum, wavelength) — the four
      matrices' spectral columns on :data:`HSAS_WAVE`;
      one (spectrum,) variable per :data:`HSAS_ANCILLARY` / :data:`HSAS_TAIL`
      name; coordinate ``time`` (UTC, from the file's day of year plus the
      decimal-hour column). Attributes: file, year, doy,
      ``ancillary_consistent`` (1 if the ancillary columns agree across the
      four matrices).
    """
    raw = readsav(path, python_dict=True, verbose=False)
    keys = {k.lower(): k for k in raw}
    mats = {}
    for name, key in HSAS_MATRICES.items():
        if key not in keys:
            raise KeyError(f"{path}: missing {key}; has {sorted(keys)}")
        a = np.asarray(raw[keys[key]], dtype=np.float64)
        # readsav reverses IDL's column-major order; accept either orientation.
        if a.ndim == 2 and a.shape[1] != _N_COLS and a.shape[0] == _N_COLS:
            a = a.T
        if a.ndim != 2 or a.shape[1] != _N_COLS:
            raise ValueError(f"{path}: {key} has shape {a.shape}, expected (n, 158)")
        mats[name] = a
    n = {m.shape[0] for m in mats.values()}
    if len(n) != 1:
        raise ValueError(f"{path}: matrices disagree in length {n}")
    ref = mats['lt']
    consistent = all(
        np.allclose(m[:, :14], ref[:, :14], equal_nan=True)
        and np.allclose(m[:, 155:], ref[:, 155:], equal_nan=True)
        for m in mats.values())
    year, doy = parse_hsas_name(path)
    hours = ref[:, 1]
    time = (pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy - 1)
            + pd.to_timedelta(hours, unit='h'))
    data = {k: (('spectrum', 'wavelength'), m[:, 14:155]) for k, m in mats.items()}
    for j, nm in enumerate(HSAS_ANCILLARY):
        data[nm] = ('spectrum', ref[:, j])
    for j, nm in enumerate(HSAS_TAIL):
        data[nm] = ('spectrum', ref[:, 155 + j])
    ds = xr.Dataset(data, coords={'wavelength': HSAS_WAVE,
                                  'time': ('spectrum', np.asarray(time, 'datetime64[ns]'))})
    ds.attrs.update(file=os.path.abspath(path), year=year, doy=doy,
                    ancillary_consistent=int(consistent))
    return ds


# --------------------------------------------------------------------------
# BODC flow cytometry (AFC) files
# --------------------------------------------------------------------------

#: Expected mapping, BODC series 2009640 (AMT24 AFC). Note 700 = Syn, 701 = Pro.
BODC_TAXON_CODES = {'P700A90Z': 'syn', 'P701A90Z': 'pro', 'PYEUA00A': 'peuk'}
_TAXON_WORDS = {'syn': ('synechococcus',), 'pro': ('prochlorococcus',),
                'peuk': ('eukaryot', 'picoeuk')}
_CODE_RE = re.compile('|'.join(BODC_TAXON_CODES))
#: Time formats tried in order; the files hold dd/mm/yyyy hh:mm GMT even
#: though the header reads 'yyyy-mm-ddThh24:mi:ss[GMT]'.
_AFC_TIME_FORMATS = ('%d/%m/%Y %H:%M', '%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M', '%Y-%m-%d %H:%M:%S',
                     '%d/%m/%Y %H:%M:%S')


def resolve_taxon_codes(headers):
    """Map BODC parameter-code columns to taxa from the deposit's own header.

    The parameter-code meaning is taken from the header text when it names the
    organism; otherwise the BODC series documentation mapping
    (:data:`BODC_TAXON_CODES`) is used. A header that names a different
    organism than the documented mapping raises.

    Inputs: headers (list of str).
    Outputs: dict code → {'taxon', 'column', 'flag_column', 'source'}.
    """
    out = {}
    for i, h in enumerate(headers):
        m = _CODE_RE.search(h)
        if m is None:
            continue
        code = m.group(0)
        if code in out:
            raise ValueError(f"parameter code {code} appears twice in the header")
        if i + 1 >= len(headers) or not headers[i + 1].strip().upper().startswith('QV'):
            raise ValueError(f"column {h!r} is not followed by its QV:BODC flag column")
        low = h.lower()
        named = [t for t, words in _TAXON_WORDS.items() if any(w in low for w in words)]
        taxon, source = BODC_TAXON_CODES[code], 'bodc_series_2009640'
        if named:
            if named != [taxon]:
                raise ValueError(f"header {h!r} names {named}, BODC docs say {taxon}")
            source = 'header'
        out[code] = {'taxon': taxon, 'column': i, 'flag_column': i + 1, 'source': source}
    missing = set(BODC_TAXON_CODES) - set(out)
    if missing:
        raise ValueError(f"AFC header lacks parameter codes {sorted(missing)}")
    return out


def _pick(headers, used, patterns, prefer=None):
    """Indices of headers matching any regex in ``patterns`` (case-insensitive),
    skipping used columns and QV columns; ``prefer`` regexes sort first."""
    hits = [i for i, h in enumerate(headers) if i not in used
            and not h.strip().upper().startswith('QV')
            and any(re.search(p, h, re.I) for p in patterns)]
    if prefer:   # ranked: a header matching an earlier preference sorts first
        rank = lambda i: next((r for r, p in enumerate(prefer) if re.search(p, headers[i], re.I)), len(prefer))
        hits.sort(key=rank)
    return hits


def _parse_number(s):
    s = s.strip()
    if s == '' or s.upper() == 'N':
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def read_afc(path, cruise, columns=None):
    """Read a BODC AFC CSV, mapping every column by header, never by position.

    Inputs
    ------
    path : str — ``AMTxx_<ship cruise>_AFC_Dataset.csv``.
    cruise : str — label stored in the output, e.g. ``'AMT24'``.
    columns : dict, optional — exact header text for any of ``station``,
        ``time``, ``lat``, ``lon``, ``depth``, ``bottle_flag`` that the
        automatic header matching gets wrong.

    Outputs
    -------
    pandas.DataFrame, one row per bottle sample: ``cruise``, ``station``,
    ``time_cast`` (UTC cast start), ``lat``, ``lon``, ``depth`` [m],
    ``bottle_flag`` (int; blank read as 0), and per taxon ``<t>`` (cells mL⁻¹,
    NaN for ``N``/blank) and ``<t>_qv`` (str QV:BODC flag). ``df.attrs`` holds
    ``taxon_map`` (:func:`resolve_taxon_codes`), ``columns`` (header used per
    field) and ``qv_counts``.
    """
    with open(path, newline='', encoding='utf-8-sig', errors='replace') as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    headers, data = [h.strip() for h in rows[0]], rows[1:]
    tmap = resolve_taxon_codes(headers)
    used = {v['column'] for v in tmap.values()} | {v['flag_column'] for v in tmap.values()}
    col = lambda i: [r[i] if i < len(r) else '' for r in data]
    columns = dict(columns or {})
    chosen = {}

    def take(field, patterns, prefer=None, test=None):
        if field in columns:
            i = headers.index(columns[field])
        else:
            hits = _pick(headers, used, patterns, prefer)
            if test is not None:
                hits = [i for i in hits if test(col(i))]
            if not hits:
                raise ValueError(f"{path}: no header for {field!r}; headers are {headers}. "
                                 f"Pass columns={{'{field}': '<header>'}}.")
            i = hits[0]
        used.add(i)
        chosen[field] = headers[i]
        return i

    def time_format(vals):
        v = pd.Series([x.strip() for x in vals if x.strip()])
        for f in _AFC_TIME_FORMATS:
            if len(v) and pd.to_datetime(v, format=f, errors='coerce').notna().mean() > 0.9:
                return f
        return None

    i_time = take('time', [r'date', r'time', r'yyyy', r'gmt'], test=lambda v: time_format(v) is not None)
    t_fmt = time_format(col(i_time))
    i_lat = take('lat', [r'lat'], prefer=[r'^lat'])
    i_lon = take('lon', [r'lon'], prefer=[r'^lon'])
    # BODC files carry 'Orig_stn' (originator's CTD station, used for the
    # verified 68-cast count) and 'BODC_stn'; prefer the originator's.
    i_stn = take('station', [r'stn', r'stat', r'cast', r'ctd'], prefer=[r'^orig', r'stn', r'stat'])
    i_flag = take('bottle_flag', [r'flag'], prefer=[r'bot'])
    # Sample depth varies within a station; sea-floor depth does not.
    stn = col(i_stn)

    def varies_within_station(vals):
        d = pd.DataFrame({'s': stn, 'd': [_parse_number(x) for x in vals]})
        return d.groupby('s')['d'].nunique().max() > 1
    if 'depth' in columns:
        i_depth = take('depth', [])
    else:
        cands = _pick(headers, used, [r'depth', r'pres'], prefer=[r'bot', r'samp', r'depth'])
        varying = [i for i in cands if varies_within_station(col(i))]
        if not varying:
            raise ValueError(f"{path}: no sample-depth header; headers are {headers}")
        i_depth = varying[0]
        used.add(i_depth)
        chosen['depth'] = headers[i_depth]

    flag = [x.strip() for x in col(i_flag)]
    df = pd.DataFrame({
        'cruise': cruise,
        'station': [x.strip() for x in stn],
        'time_cast': pd.to_datetime(pd.Series(col(i_time)).str.strip(),
                                    format=t_fmt, errors='coerce'),
        'lat': [_parse_number(x) for x in col(i_lat)],
        'lon': [_parse_number(x) for x in col(i_lon)],
        'depth': [_parse_number(x) for x in col(i_depth)],
        'bottle_flag': [int(float(x)) if x not in ('', 'N') else 0 for x in flag],
    })
    qv_counts = {}
    for code, info in tmap.items():
        t = info['taxon']
        df[t] = [_parse_number(x) for x in col(info['column'])]
        df[t + '_qv'] = [x.strip() for x in col(info['flag_column'])]
        qv_counts[t] = df[t + '_qv'].value_counts().to_dict()
    df.attrs.update(taxon_map=tmap, columns=chosen, qv_counts=qv_counts, time_format=t_fmt,
                    file=os.path.abspath(path))
    return df


def afc_label_status(df, qv_missing=('N',)):
    """Classify each taxon value as observed, censored or missing.

    Rules (Stage 5): ``N`` (or blank) → missing; a ``<`` qualifier or a zero
    abundance → left-censored at the detection limit; any other non-blank QV
    flag listed in ``qv_missing`` → missing. Other QV flags are kept as
    observed and counted in ``df.attrs['qv_counts']``.

    Inputs: df from :func:`read_afc`; qv_missing (tuple of str).
    Outputs: dict taxon → dict of bool arrays ``observed``, ``censored``,
    ``missing``.
    """
    out = {}
    for t in TAXA:
        v = df[t].to_numpy(float)
        q = df[t + '_qv'].astype(str).str.strip().to_numpy()
        missing = ~np.isfinite(v) | np.isin(q, qv_missing)
        censored = ~missing & ((v == 0) | (q == '<'))
        out[t] = {'observed': ~missing & ~censored, 'censored': censored,
                  'missing': missing}
    return out


def _lattice_prob(step, n_int, smear, powers):
    """P(v | step) for v = 1..n_int: v = round(k·step·(1+δ)), δ ~ U(−smear, smear)
    (per-sample analysed-volume spread), k ≥ 1 with p(k) ∝ k^−a, conditional
    on v ≤ n_int. Returns (len(powers), n_int)."""
    kmax = int(np.ceil((n_int + 0.5) / (step * (1 - smear)))) + 1
    k = np.arange(1, kmax + 1)[:, None]
    lo, hi = k * step * (1 - smear), k * step * (1 + smear)
    v = np.arange(1, n_int + 1)[None, :]
    if smear > 0:
        M = np.clip(np.minimum(hi, v + 0.5) - np.maximum(lo, v - 0.5), 0, None) / (hi - lo)
    else:
        M = (np.floor(k * step + 0.5) == v).astype(float)
    W = k[:, 0][None, :] ** -np.asarray(powers, float)[:, None]          # (A, K)
    P = W @ M
    return P / np.maximum(P.sum(1, keepdims=True), 1e-300)


def _continuous_loglik(vi, n_int, eps, powers, amins=np.arange(0.5, 1.51, 0.05)):
    """Best log-likelihood of the *unquantised* null: abundance continuous with
    density ∝ a^−b on [a_min, ∞), rounded to whole cells, conditional on
    v ≤ n_int; b and a_min profiled."""
    v = np.arange(1, n_int + 1, dtype=float)
    best = -np.inf
    for b in powers:
        F = (lambda a: np.log(a)) if abs(b - 1) < 1e-9 else (lambda a, b=b: a ** (1 - b) / (1 - b))
        for am in amins:
            lo, hi = np.maximum(v - 0.5, am), v + 0.5
            P = np.where(hi > lo, F(hi) - F(np.where(hi > lo, lo, hi)), 0.0)
            if P.sum() <= 0:
                continue
            P = (1 - eps) * P / P.sum() + eps / n_int
            best = max(best, float(np.log(P[vi - 1]).sum()))
    return best


def infer_count_step(values, max_value=20.0, step_grid=None, eps=0.01, min_gain=6.0,
                     on_fraction=0.8, smear=0.035, powers=(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)):
    """Smallest non-zero abundance step (one event = 1/V) from quantisation.

    Abundances are event counts divided by the analysed volume V, so one event
    sets the detection limit. Fractional files (AMT23) show the step directly
    as the cluster of single-event values. Files rounded to whole cells are fit
    by maximum likelihood with a *smeared* lattice: ``v = round(k·step·(1+δ))``
    with δ ~ U(±smear) for the per-sample volume spread (AMT23's single-event
    values span 1.31–1.40, about ±3.5 %), and p(k) ∝ k^−a with the exponent
    profiled over ``powers`` (event counts pile up at small k). A rigid
    lattice with uniform k aliases the step to about half its value, because
    the spread fills the predicted gaps at high counts and a half step absorbs
    the pile-up at one event. An ``eps`` uniform outlier term keeps the
    likelihood finite. The best step is compared with an unquantised null
    (continuous power-law abundance with a fitted lower cut-off, rounded to
    whole cells); if quantisation is not detected (gain < ``min_gain`` nats,
    calibrated on simulated null data, or fewer than
    ``on_fraction`` of fractional values on the lattice below three events),
    the fallback is the smallest positive value. The configured steps in
    DEFAULT_PIPELINE are what the pipeline uses; this fit is the check.

    Inputs
    ------
    values : array — abundances (any taxa; zeros and NaN ignored).

    Outputs
    -------
    dict — ``step``, ``volume_ml`` (1/step), ``method`` ('direct', 'fit' or
    'fallback_min'), ``delta_loglik`` (fit vs the unquantised null), ``second_step``
    (best local maximum ≥ 0.05 away; NaN if none), ``power`` (fitted a),
    ``n_values``, ``hist`` (counts of the values 1..19, for inspection).
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v) & (v > 0) & (v < max_value)]
    out = {'n_values': int(v.size), 'delta_loglik': np.nan, 'second_step': np.nan, 'power': np.nan}
    if v.size == 0:
        raise ValueError("no positive values below max_value")
    if np.any(np.abs(v - np.round(v)) > 1e-6):
        single = v[v < 1.5 * v.min()]
        step = float(np.median(single))
        # per-sample volumes differ by a few %, so the lattice is only checked
        # below three events, where it cannot smear into the next multiple
        low = v[v < 3 * step]
        q = np.abs(low / step - np.floor(low / step + 0.5))
        if np.mean(q < 0.15) >= on_fraction:
            out.update(step=step, method='direct')
        else:
            out.update(step=float(v.min()), method='fallback_min')
        out['volume_ml'] = 1.0 / out['step']
        return out
    vi = np.floor(v + 0.5).astype(int)
    n_int = int(np.ceil(max_value)) - 1            # attainable integers 1..n_int
    out['hist'] = np.bincount(vi, minlength=n_int + 1)[1:].tolist()
    grid = np.arange(0.5, 2.5, 0.001) if step_grid is None else np.asarray(step_grid)

    def loglik(step):
        P = (1 - eps) * _lattice_prob(step, n_int, smear, powers) + eps / n_int
        L = np.log(P[:, vi - 1]).sum(1)
        j = int(np.argmax(L))
        return float(L[j]), powers[j]
    res = [loglik(st) for st in grid]
    ll = np.array([r[0] for r in res])
    ll1 = _continuous_loglik(vi, n_int, eps, powers)
    i = int(np.argmax(ll))
    gain = ll[i] - ll1
    peaks = [j for j in range(1, len(ll) - 1) if ll[j] >= ll[j - 1] and ll[j] >= ll[j + 1]
             and abs(grid[j] - grid[i]) >= 0.05]
    if peaks:
        out['second_step'] = float(grid[max(peaks, key=lambda j: ll[j])])
    out['delta_loglik'] = float(gain)
    out['power'] = float(res[i][1])
    if gain < min_gain:
        out.update(step=float(v.min()), method='fallback_min')
    else:
        out.update(step=float(grid[i]), method='fit')
    out['volume_ml'] = 1.0 / out['step']
    return out


# --------------------------------------------------------------------------
# Jordan et al. 2025 underway SST
# --------------------------------------------------------------------------

def read_jordan_sst(path, var='uway_sst'):
    """Underway SST time series from the Jordan et al. 2025 netCDF.

    Inputs: path (str); var (str, default ``uway_sst``).
    Outputs: pandas.Series of SST [°C] indexed by UTC time, sorted, finite only.
    """
    with xr.open_dataset(path) as ds:
        if var not in ds:
            raise KeyError(f"{path}: no variable {var!r}; has {list(ds.data_vars)}")
        da = ds[var].squeeze()
        if da.ndim != 1:
            raise ValueError(f"{var} has dims {da.dims}; expected one time-like dim")
        dim = da.dims[0]
        names = [c for c in list(ds.coords) + list(ds.data_vars)
                 if ds[c].dims == (dim,) and np.issubdtype(ds[c].dtype, np.datetime64)]
        if not names:
            raise ValueError(f"{path}: no decoded datetime variable along {dim!r}")
        names.sort(key=lambda c: c.lower() != 'time')
        t = pd.to_datetime(ds[names[0]].values)
        s = pd.Series(da.values.astype(float), index=t)
    s = s[np.isfinite(s.values) & s.index.notna()].sort_index()
    return s[~s.index.duplicated()]


def interp_series(series, times):
    """Linear interpolation of a time series; NaN outside its range.

    Inputs: series (pandas.Series, datetime index); times (array-like datetimes).
    Outputs: (values, gap_min) — interpolated values and the distance in
    minutes to the nearest valid sample (a diagnostic, not a filter).
    """
    x = series.index.values.astype('datetime64[ns]').astype(np.int64).astype(float)
    t = pd.to_datetime(np.asarray(times)).values.astype('datetime64[ns]').astype(np.int64).astype(float)
    val = np.interp(t, x, series.values, left=np.nan, right=np.nan)
    j = np.clip(np.searchsorted(x, t), 1, len(x) - 1)
    gap = np.minimum(np.abs(t - x[j - 1]), np.abs(x[j] - t)) / 6e10
    return val, gap


# --------------------------------------------------------------------------
# Brewin et al. 2023 compilation
# --------------------------------------------------------------------------

_RRS_RE = re.compile(r'^Rrs\((\d+(?:\.\d+)?)\)\s*\[1/sr\]$')
_UNC_RE = re.compile(r'^Rrs\((\d+(?:\.\d+)?)\)\s*Uncertainty')
_BREWIN_META = {'cruise': 'AMT Cruise', 'date': 'Year-Month-Day', 'time': 'Time [GMT]',
                'lat': 'Latitude [deg+veN]', 'lon': 'Longitude [deg+veE]',
                'sza': 'Solar Zenith Angle [deg]',
                'sst': 'Temperature above Secchi depth [degC]',
                'ctd_method': 'CTD Method (1=Profiler; 2=Underway)'}
_BREWIN_QC = (r'cloud', r'glint', r'pitch', r'roll', r'Li\s*\(?750', r'750.*Es')


def read_brewin(path, sentinel=-999.0):
    """Read the Brewin et al. 2023 AMT optical compilation.

    Uses the uncorrected ``Rrs(λ) [1/sr]`` columns, not the BRDF-corrected
    ones (the AMT24 pipeline applies no BRDF correction, §8.5). The signed
    percentage uncertainty is converted to absolute,
    ``σ_file = |U/100 · Rrs|``.

    Inputs: path (str); sentinel (missing-value code, −999).
    Outputs: dict — ``meta`` (DataFrame: cruise label 'AMTnn', time [UTC],
    lat, lon, sza, sst, ctd_method, plus ancillary QC columns carried for
    stratified reporting), ``wave`` (151,), ``rrs``, ``unc_pct``,
    ``sigma_file`` (n, 151; NaN for missing), ``usable`` (n,) bool (all bands
    valid), ``rrs_columns``, ``unc_columns`` (selected headers).
    """
    df = pd.read_csv(path)
    rcols, ucols = {}, {}
    for c in df.columns:
        if 'BDRF' in c.upper() or 'BRDF' in c.upper():
            continue
        m = _RRS_RE.match(c.strip())
        if m:
            rcols[float(m.group(1))] = c
        m = _UNC_RE.match(c.strip())
        if m:
            ucols[float(m.group(1))] = c
    wave = np.array(sorted(rcols))
    if wave.size == 0:
        raise ValueError(f"{path}: no uncorrected Rrs columns found")
    R = df[[rcols[w] for w in wave]].to_numpy(float)
    R[R <= sentinel + 1e-6] = np.nan
    U = np.full_like(R, np.nan)
    for j, w in enumerate(wave):
        if w in ucols:
            U[:, j] = df[ucols[w]].to_numpy(float)
    U[U <= sentinel + 1e-6] = np.nan
    meta = pd.DataFrame({
        'cruise': ['AMT%d' % int(c) for c in df[_BREWIN_META['cruise']]],
        'time': pd.to_datetime(df[_BREWIN_META['date']].astype(str).str.strip() + ' '
                               + df[_BREWIN_META['time']].astype(str).str.strip(),
                               errors='coerce'),
        'lat': df[_BREWIN_META['lat']].to_numpy(float),
        'lon': df[_BREWIN_META['lon']].to_numpy(float),
    })
    for k in ('sza', 'sst', 'ctd_method'):
        if _BREWIN_META[k] in df:
            v = df[_BREWIN_META[k]].to_numpy(float)
            meta[k] = np.where(v <= sentinel + 1e-6, np.nan, v)
    for c in df.columns:
        if c not in rcols.values() and c not in ucols.values() and not c.startswith('Rrs(') \
                and any(re.search(p, c, re.I) for p in _BREWIN_QC):
            meta['qc:' + c] = df[c].to_numpy()
    return {'meta': meta, 'wave': wave, 'rrs': R, 'unc_pct': U,
            'sigma_file': np.abs(U / 100.0 * R), 'usable': np.isfinite(R).all(axis=1),
            'rrs_columns': [rcols[w] for w in wave],
            'unc_columns': [ucols.get(w) for w in wave]}


def fit_sigma_floor(rrs, unc_pct):
    """Robust fit of ``σ = sqrt((f_rel·Rrs)² + σ_abs²)`` to absolute
    uncertainties (§5.2; reproduces f_rel ≈ 0.027, σ_abs ≈ 8.2×10⁻⁵ sr⁻¹ on the
    Brewin file).

    Inputs: rrs, unc_pct — arrays of Rrs [sr⁻¹] and signed % uncertainty.
    Outputs: dict ``f_rel``, ``sigma_abs``, ``log_resid_p16``, ``_p50``,
    ``_p84``, ``n``.
    """
    R = np.asarray(rrs, float)
    U = np.asarray(unc_pct, float)
    A = np.abs(U / 100.0 * R)
    ok = np.isfinite(R) & np.isfinite(U) & (np.abs(R) > 0) & (A > 0)
    A, Rv = A[ok], np.abs(R[ok])

    def res(p):
        f, s = np.exp(p)
        return np.log(np.sqrt((f * Rv) ** 2 + s ** 2)) - np.log(A)
    fit = least_squares(res, np.log([0.03, 7e-5]), loss='soft_l1')
    f, s = np.exp(fit.x)
    r = res(fit.x)
    p16, p50, p84 = np.percentile(r, [16, 50, 84])
    return {'f_rel': float(f), 'sigma_abs': float(s), 'log_resid_p16': float(p16),
            'log_resid_p50': float(p50), 'log_resid_p84': float(p84), 'n': int(ok.sum())}


# --------------------------------------------------------------------------
# SeaFlow v2 and PACE OCI L2 (§8.6, test only this stage)
# --------------------------------------------------------------------------

_SEAFLOW_GUESS = {'time': [r'^time$', r'date'], 'lat': [r'^lat'], 'lon': [r'^lon'],
                  'sst': [r'^temp', r'sst'],
                  'pro': [r'prochloro'], 'syn': [r'synecho'], 'peuk': [r'picoeuk']}


def read_seaflow(path, columns=None, abundance_factor=1000.0):
    """Read a SeaFlow v2 abundance table (CSV, TSV or parquet).

    Inputs
    ------
    path : str.
    columns : dict, optional — exact column names for ``time``, ``lat``,
        ``lon``, ``sst``, ``pro``, ``syn``, ``peuk`` (guessed otherwise; an
        abundance column containing 'abundance' is preferred).
    abundance_factor : float — multiplies the file's abundances to give cells
        mL⁻¹. Default 1000 assumes cells µL⁻¹; check the file's units.

    Outputs: DataFrame with ``time`` (UTC), ``lat``, ``lon``, ``sst``, ``pro``,
    ``syn``, ``peuk`` [cells mL⁻¹]; ``attrs['columns']`` records the mapping.
    """
    if str(path).endswith('.parquet'):
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep=None, engine='python')
    columns = dict(columns or {})
    chosen = {}
    for k, pats in _SEAFLOW_GUESS.items():
        if k in columns:
            chosen[k] = columns[k]
            continue
        hits = [c for c in df.columns if any(re.search(p, c, re.I) for p in pats)]
        hits.sort(key=lambda c: 'abundance' not in c.lower())
        if not hits:
            raise ValueError(f"{path}: no column for {k!r}; columns: {list(df.columns)}")
        chosen[k] = hits[0]
    out = pd.DataFrame({'time': pd.to_datetime(df[chosen['time']], utc=True).dt.tz_localize(None),
                        'lat': df[chosen['lat']].astype(float),
                        'lon': df[chosen['lon']].astype(float),
                        'sst': df[chosen['sst']].astype(float)})
    for t in TAXA:
        out[t] = df[chosen[t]].astype(float) * abundance_factor
    out.attrs['columns'] = chosen
    out.attrs['abundance_factor'] = abundance_factor
    return out


#: Standard L2 flags excluded in the PACE box (§8.6).
PACE_L2_FLAGS = ('ATMFAIL', 'LAND', 'HIGLINT', 'HILT', 'HISATZEN', 'STRAYLIGHT',
                 'CLDICE', 'COCCOLITH', 'HISOLZEN', 'LOWLW', 'CHLFAIL', 'NAVWARN',
                 'MAXAERITER', 'CHLWARN', 'ATMWARN', 'NAVFAIL')


def read_pace_box(path, lat, lon, half=2, flags=PACE_L2_FLAGS):
    """Extract the (2·half+1)² pixel box nearest a position from a PACE OCI L2
    AOP granule.

    Inputs: path (netCDF L2 granule); lat, lon (deg); half (int, 2 → 5×5);
    flags (L2 flag names that invalidate a pixel).

    Outputs: dict or None (position not in granule) — ``wave`` (nλ,),
    ``rrs`` and ``rrs_unc`` (n_pix, nλ; NaN where flagged or absent),
    ``valid`` (n_pix,) bool, ``valid_fraction``, ``dist_km`` to the nearest
    pixel, ``time`` (granule scan time of the centre line, UTC).
    """
    with xr.open_dataset(path, group='navigation_data') as nav, \
            xr.open_dataset(path, group='geophysical_data') as geo, \
            xr.open_dataset(path, group='sensor_band_parameters') as sbp:
        la, lo = nav['latitude'].values, nav['longitude'].values
        d = _haversine_km(lat, lon, la, lo)
        if not np.isfinite(d).any():
            return None
        i, j = np.unravel_index(np.nanargmin(d), d.shape)
        if i < half or j < half or i + half >= d.shape[0] or j + half >= d.shape[1]:
            return None
        sl = (slice(i - half, i + half + 1), slice(j - half, j + half + 1))
        rrs = geo['Rrs'].values[sl].reshape(-1, geo['Rrs'].shape[-1]).astype(float)
        unc = (geo['Rrs_unc'].values[sl].reshape(rrs.shape).astype(float)
               if 'Rrs_unc' in geo else np.full_like(rrs, np.nan))
        fl = geo['l2_flags']
        bits = dict(zip(fl.attrs.get('flag_meanings', '').split(),
                        np.atleast_1d(fl.attrs.get('flag_masks', []))))
        mask = 0
        for f in flags:
            mask |= int(bits.get(f, 0))
        bad = (fl.values[sl].reshape(-1).astype(np.int64) & mask) != 0
        wave = np.asarray(sbp['wavelength_3d'].values, float)
    valid = ~bad & np.isfinite(rrs).any(axis=1)
    rrs[~valid] = np.nan
    unc[~valid] = np.nan
    with xr.open_dataset(path, group='scan_line_attributes') as sla:
        t = None
        if 'time' in sla:
            tv = sla['time'].values[i]
            t = pd.Timestamp(tv) if np.issubdtype(np.asarray(tv).dtype, np.datetime64) else None
    if t is None:
        with xr.open_dataset(path) as root:
            t = pd.Timestamp(root.attrs.get('time_coverage_start', 'NaT')).tz_localize(None)
    return {'wave': wave, 'rrs': rrs, 'rrs_unc': unc, 'valid': valid,
            'valid_fraction': float(valid.mean()), 'dist_km': float(d[i, j]), 'time': t}


def _haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance [km]; arrays broadcast."""
    a1, o1, a2, o2 = map(np.radians, (lat1, lon1, lat2, lon2))
    h = np.sin((a2 - a1) / 2) ** 2 + np.cos(a1) * np.cos(a2) * np.sin((o2 - o1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))


haversine_km = _haversine_km


def sha256_file(path, chunk=1 << 22):
    """SHA-256 of a file (provenance). Inputs: path. Outputs: hex str."""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(chunk), b''):
            h.update(b)
    return h.hexdigest()
