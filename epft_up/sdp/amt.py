"""AMT above-water Rrs + HPLC matchups (Tom Jordan's compilation): loader and harmonisation.

The file ``AMT_Rrs_HPLC.csv`` (received 2026-10-24; staged by
``scripts/sdp/stage_amt.py`` in ``$OS_COLOR/AMT/Rrs_HPLC/``) holds 210 HPLC
samples from 8 AMT cruises, each matched to a 400-700 nm, 1 nm Rrs spectrum.
:func:`load_amt` checks the pinned checksum and, unless ``raw=True``, applies
the harmonisation rules of ``reports/provenance/AMT_benchmark_protocol.md``
§3, in this order:

1. drop Sample_ID 107/108 (one spectrum, Tchla 0.073 vs 0.165: a likely
   mismatch; :data:`DROP_SAMPLE_IDS`);
2. exclude spectra with any Rrs ≤ 0;
3. collapse spectra shared by several HPLC rows to one row (pigments, T and
   S averaged);
4. replace implausible open-ocean salinity (< :data:`SAL_MIN_OPEN_OCEAN`)
   with WOA23 climatology (``woa23_salinity_fix.csv``, written by
   ``scripts/sdp/amt_harmonise.py``);
5. fill missing MVchlb with Tchlb (MVchlb + DVchlb) and flag it; MVchlb is
   to be *scored* only where measured (AMT29/30);
6. zero the 13 modelled pigments below the deposit's LOD proxy
   (:data:`LOD_PROXY_DEPOSIT`), as Kramer et al. (2022) zeroed below-LOD
   values. Ratios then use ½·LOD via ``validate.replace_zeros``.

The 5 nm moving mean is **not** applied here. It belongs to the modelling
pipeline (``spectral.moving_mean(..., edges='shrink')``), so the loader
returns the spectra as delivered.
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from epft_up import __version__
from epft_up.sdp import data as sdpdata

#: SHA-256 of the 2026-10-24 delivery (``scripts/sdp/inspect_amt.py``).
AMT_SHA256 = '0290c9f49fae79dd52dbf920e77ff580c15af9115cd26ed734659d2c7827bca6'
AMT_FILE = 'AMT_Rrs_HPLC.csv'
SALINITY_FIX_FILE = 'woa23_salinity_fix.csv'

#: AMT column -> our pigment name (data.PIGMENTS_13 plus extras)
PIG_MAP = {'Tchla': 'Tchla', 'Hex': 'HexFuco', 'But': 'ButFuco', 'Allo': 'Allo',
           'Fuco': 'Fuco', 'Perid': 'Perid', 'Zea': 'Zea', 'DVchla': 'DVchla',
           'MVchlb': 'MVchlb', 'Chlc12': 'Chlc12', 'Chlc3': 'Chlc3', 'Neo': 'Neo',
           'Viola': 'Viola', 'Tchlb': 'Tchlb', 'Tchlc': 'Tchlc', 'ABcaro': 'ABcaro',
           'Diadino': 'Diadino', 'Diato': 'Diato', 'MVchla': 'MVchla', 'DVchlb': 'DVchlb',
           'Lut': 'Lut', 'Pras': 'Pras', 'Psc': 'Psc', 'Ppc': 'Ppc'}
META_COLS = ('Sample_ID', 'Original_ID', 'Cruise', 'Data_link', 'HPLC_lab', 'Rad_method',
             'DateTime', 'Lat', 'Lon', 'Temp', 'Sal')

#: Shared spectrum whose two HPLC samples disagree ×2.3 in Tchla (asked of Tom Jordan).
DROP_SAMPLE_IDS = (107, 108)
#: Below this, a surface salinity on the AMT track is treated as a bad value.
SAL_MIN_OPEN_OCEAN = 30.0
#: ``validate.lod_proxy`` of the Kramer deposit [mg m^-3] (protocol §3; a test
#: checks it against the deposit). Fixed here so both sets use one rule.
LOD_PROXY_DEPOSIT = {'Tchla': 0.001, 'HexFuco': 0.003, 'ButFuco': 0.002, 'Allo': 0.001,
                     'Fuco': 0.002, 'Perid': 0.004, 'Zea': 0.001, 'DVchla': 0.002,
                     'MVchlb': 0.003, 'Chlc12': 0.001, 'Chlc3': 0.001, 'Neo': 0.002,
                     'Viola': 0.001}
#: Cruises whose HPLC lab (NASA GSFC) reported MVchlb; score MVchlb only there.
MVCHLB_MEASURED_CRUISES = ('AMT29', 'AMT30')


def amt_dir():
    """``$OS_COLOR/AMT/Rrs_HPLC``."""
    return Path(os.environ['OS_COLOR']) / 'AMT' / 'Rrs_HPLC'


@dataclass
class AMTData:
    """AMT matchups, raw or harmonised.

    Attributes
    ----------
    wave : ndarray, shape (n_wave,)
        400-700 nm, 1 nm.
    Rrs : ndarray, shape (n, n_wave)
        Above-water Rrs [sr^-1], as delivered (not smoothed).
    meta : DataFrame
        ``sample_id`` (``'AMT:<Sample_ID>'``, the smallest ID of a collapsed
        group), ``sample_ids`` (all IDs, ``;``-joined), ``original_ids``,
        ``cruise`` (the CV group), ``hplc_lab``, ``time`` (UTC), ``lat``,
        ``lon``, ``temp``, ``sal``, ``sal_source``, ``n_hplc``.
    pigments : DataFrame
        Our pigment names [mg m^-3]; zeros are below the LOD proxy.
    flags : DataFrame
        ``MVchlb_proxy`` (filled from Tchlb), ``MVchlb_scorable``,
        ``sal_replaced``, and ``lod_zeroed:<pigment>``.
    steps : DataFrame
        Rows per cruise after each harmonisation rule.
    provenance : dict
    """
    wave: np.ndarray
    Rrs: np.ndarray
    meta: pd.DataFrame
    pigments: pd.DataFrame
    flags: pd.DataFrame
    steps: pd.DataFrame
    provenance: dict = field(default_factory=dict)

    @property
    def n(self):
        return self.Rrs.shape[0]

    @property
    def pigments13(self):
        return self.pigments[list(sdpdata.PIGMENTS_13)]


def spectrum_key(Rrs):
    """Deterministic per-row key identifying identical spectra."""
    R = np.round(np.asarray(Rrs, float), 10)
    return np.array([hashlib.sha1(r.tobytes()).hexdigest() for r in R])


def read_csv(path=None, verify=True):
    """The delivered CSV as a DataFrame, after checking its pinned SHA-256."""
    path = Path(path) if path is not None else amt_dir() / AMT_FILE
    sha = sdpdata.sha256_of(path)
    if verify and sha != AMT_SHA256:
        raise ValueError(f'{path}: SHA-256 {sha} is not the pinned {AMT_SHA256}; '
                         'stage a new delivery under a new name and re-pin it')
    return pd.read_csv(path, na_values=['nan', 'NaN', '']), path, sha


def read_salinity_fix(path=None):
    """``{Sample_ID: S_woa}`` from ``woa23_salinity_fix.csv``."""
    path = Path(path) if path is not None else amt_dir() / SALINITY_FIX_FILE
    f = pd.read_csv(path)
    return dict(zip(f['Sample_ID'].astype(int), f['S_woa'].astype(float)))


def _count(df, step, rows):
    c = df.groupby('Cruise').size()
    c['total'] = int(len(df))
    rows.append(c.rename(step))


def harmonise(df, sal_fix=None, lods=None):
    """Apply protocol §3 to the raw table (see module docstring for the order).

    Parameters
    ----------
    df : DataFrame
        As read by :func:`read_csv`.
    sal_fix : dict, optional
        ``{Sample_ID: salinity}`` for rows with ``Sal < SAL_MIN_OPEN_OCEAN``.
        Required when such rows survive the earlier rules.
    lods : dict, optional
        Per-pigment LOD; default :data:`LOD_PROXY_DEPOSIT`.

    Returns
    -------
    out : DataFrame
        One row per spectrum with meta, pigments (our names), Rrs columns and
        flag columns.
    steps : DataFrame
        Rows per cruise (columns) after each rule (rows).
    """
    lods = LOD_PROXY_DEPOSIT if lods is None else lods
    rcols = [c for c in df.columns if c.startswith('Rrs_')]
    pig_in = [c for c in PIG_MAP if c in df.columns]
    rows = []
    d = df.copy()
    _count(d, '0 raw rows', rows)
    d = d[~d['Sample_ID'].isin(DROP_SAMPLE_IDS)]
    _count(d, '1 drop 107/108 (mismatched pair)', rows)
    d = d[(d[rcols].to_numpy(float) > 0).all(axis=1)]
    _count(d, '2 exclude non-positive Rrs', rows)

    d = d.assign(_key=spectrum_key(d[rcols].to_numpy(float)))
    agg = {c: 'mean' for c in pig_in + ['Temp', 'Sal', 'Lat', 'Lon']}
    agg.update({c: 'first' for c in rcols + ['Cruise', 'HPLC_lab', 'Data_link',
                                              'Rad_method', 'DateTime']})
    g = d.groupby('_key', sort=False)
    out = g.agg(agg)
    out['Sample_ID'] = g['Sample_ID'].min()
    out['sample_ids'] = g['Sample_ID'].apply(lambda s: ';'.join(map(str, sorted(s))))
    out['original_ids'] = g['Original_ID'].apply(
        lambda s: ';'.join('' if pd.isna(x) else str(x) for x in s))
    out['n_hplc'] = g.size()
    out = out.reset_index(drop=True).sort_values('Sample_ID').reset_index(drop=True)
    _count(out, '3 one row per spectrum (HPLC averaged)', rows)

    bad = out['Sal'] < SAL_MIN_OPEN_OCEAN
    out['sal_source'] = 'in situ'
    out['sal_replaced'] = False
    if bad.any():
        if sal_fix is None:
            raise ValueError('implausible salinity present; pass sal_fix '
                             '(see scripts/sdp/amt_harmonise.py)')
        for i in np.where(bad)[0]:
            ids = [int(x) for x in out.loc[i, 'sample_ids'].split(';')]
            vals = [sal_fix[x] for x in ids if x in sal_fix]
            if not vals:
                raise KeyError(f'no salinity fix for Sample_ID {ids}')
            out.loc[i, 'Sal'] = float(np.mean(vals))
            out.loc[i, 'sal_source'] = 'WOA23'
            out.loc[i, 'sal_replaced'] = True
    _count(out, '4 salinity fixed (count unchanged)', rows)

    out = out.rename(columns=PIG_MAP)
    proxy = out['MVchlb'].isna()
    out['MVchlb_proxy'] = proxy
    out.loc[proxy, 'MVchlb'] = out.loc[proxy, 'Tchlb']
    out['MVchlb_scorable'] = out['Cruise'].isin(MVCHLB_MEASURED_CRUISES) & ~proxy
    _count(out, '5 MVchlb proxy filled (count unchanged)', rows)

    for p, lod in lods.items():
        z = (out[p] < lod) & (out[p] > 0)
        out[f'lod_zeroed:{p}'] = z | (out[p] == 0)
        out.loc[z, p] = 0.0
    _count(out, '6 LOD-proxy zeroing (count unchanged)', rows)
    steps = pd.DataFrame(rows).fillna(0).astype(int)
    return out, steps


def _to_struct(out, rcols, provenance, steps):
    wave = np.array([float(c[4:]) for c in rcols])
    pig = [p for p in PIG_MAP.values() if p in out.columns]
    meta = pd.DataFrame({
        'sample_id': 'AMT:' + out['Sample_ID'].astype(str),
        'sample_ids': out.get('sample_ids', out['Sample_ID'].astype(str)),
        'original_ids': out.get('original_ids', out.get('Original_ID')),
        'cruise': out['Cruise'], 'hplc_lab': out['HPLC_lab'],
        'time': pd.to_datetime(out['DateTime'], format='%m/%d/%Y %H:%M', utc=True),
        'lat': out['Lat'], 'lon': out['Lon'], 'temp': out['Temp'], 'sal': out['Sal'],
        'sal_source': out.get('sal_source', 'in situ'),
        'n_hplc': out.get('n_hplc', 1)})
    fcols = [c for c in out.columns
             if c in ('MVchlb_proxy', 'MVchlb_scorable', 'sal_replaced')
             or c.startswith('lod_zeroed:')]
    return AMTData(wave=wave, Rrs=out[rcols].to_numpy(float), meta=meta.reset_index(drop=True),
                   pigments=out[pig].reset_index(drop=True),
                   flags=out[fcols].reset_index(drop=True), steps=steps,
                   provenance=provenance)


def load_amt(raw=False, path=None, verify=True, sal_fix=None, lods=None):
    """Load the AMT matchups.

    Parameters
    ----------
    raw : bool, optional
        True: all delivered rows, renamed but otherwise untouched.
    path : str or Path, optional
        Default ``amt_dir() / AMT_FILE``.
    verify : bool, optional
        Refuse a file whose SHA-256 differs from :data:`AMT_SHA256`.
    sal_fix : dict, optional
        Default: read ``woa23_salinity_fix.csv`` next to the CSV when needed.
    lods : dict, optional
        Default :data:`LOD_PROXY_DEPOSIT`.

    Returns
    -------
    AMTData
    """
    df, p, sha = read_csv(path, verify=verify)
    rcols = [c for c in df.columns if c.startswith('Rrs_')]
    prov = {'source': 'AMT_Rrs_HPLC.csv, compiled by Tom Jordan (received 2026-10-24)',
            'file': str(p), 'sha256': sha, 'sha256_pinned': AMT_SHA256,
            'sha256_verified': sha == AMT_SHA256, 'raw': raw,
            'loader': 'epft_up.sdp.amt.load_amt', 'epft_up_version': __version__,
            'protocol': 'reports/provenance/AMT_benchmark_protocol.md §3'}
    if raw:
        steps = pd.DataFrame([df.groupby('Cruise').size().rename('0 raw rows')])
        return _to_struct(df.rename(columns=PIG_MAP), rcols, prov, steps)
    if sal_fix is None and (df['Sal'] < SAL_MIN_OPEN_OCEAN).any():
        fp = Path(p).parent / SALINITY_FIX_FILE
        if not fp.exists():
            raise FileNotFoundError(f'{fp} missing; run scripts/sdp/amt_harmonise.py')
        sal_fix = read_salinity_fix(fp)
        prov['salinity_fix'] = {'file': str(fp), 'sha256': sdpdata.sha256_of(fp)}
    out, steps = harmonise(df, sal_fix=sal_fix, lods=lods)
    prov['rules'] = list(steps.index)
    prov['lod_proxy'] = dict(LOD_PROXY_DEPOSIT if lods is None else lods)
    prov['dropped_sample_ids'] = list(DROP_SAMPLE_IDS)
    return _to_struct(out, rcols, prov, steps)
