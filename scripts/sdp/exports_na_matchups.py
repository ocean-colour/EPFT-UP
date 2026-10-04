"""Execution #9: build the EXPORTS-NA (May 2021) HPLC matchups for Kramer's 17 test spectra.

Inputs (both pinned in the output provenance):

* the 17 spectra of ``Kramer_rrs_testdata.mat`` (``sashajane19/Rrs_pigments``
  @ pinned commit): Rrs 400-700 nm, T, S, lat/lon and HPLC chlorophyll-a,
  with **no sample times**;
* the SeaBASS files the user downloaded by hand (``docs/HOWTO_SeaBASS_EXPORTS_NA.md``)
  into ``$OS_COLOR/SeaBASS/EXPORTS``. The three UCSB/CRSEO rosette HPLC files
  (JC214, DY131, DY130) are the ones used here.

Findings that shaped the method (report §9):

* The test spectra are **not** NASA GSFC's DY131 HyperSAS Rrs. No HyperSAS
  spectrum matches any of them (best relative RMS 5-13%). Their source
  radiometer is not in this SeaBASS download.
* The test **positions are HPLC sample positions**, and the test ``chl`` is
  exactly the **mean of the replicate surface HPLC Tchla** at that position
  (e.g. 0.998 = (0.960 + 1.036)/2).

Matching rule: average the replicates of each (ship, cast time) HPLC sample
at the shallowest bottle (≤ 12 m). For each test spectrum, keep the
candidates within 1 km whose replicate-mean Tchla equals the test ``chl`` to
≤ 0.0006 mg m⁻³ (rounding). Exactly one candidate is required; anything else
raises (stop and ask). Below-detection values (-8888) become 0, as Kramer
does; missing values (-9999) stay NaN.

Output: ``$OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA/exports_na_matchups.csv``
plus a JSON provenance sidecar.

Run with::

    conda run -n ocean14 python scripts/sdp/exports_na_matchups.py
"""
import datetime
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))

from epft_up import __version__  # noqa: E402
from epft_up.sdp import data as sdpdata, seabass  # noqa: E402

SEABASS_DIR = Path(os.environ['OS_COLOR']) / 'SeaBASS' / 'EXPORTS'
OUT_DIR = sdpdata.kramer2022_dir() / 'EXPORTS_NA'
MAX_KM = 1.0
CHL_TOL = 6e-4
MAX_DEPTH = 12.0
#: SeaBASS HPLC field -> our pigment abbreviation (data.PIGMENTS_13 + extras).
FIELD_MAP = {'Tot_Chl_a': 'Tchla', 'Hex-fuco': 'HexFuco', 'But-fuco': 'ButFuco',
             'Allo': 'Allo', 'Fuco': 'Fuco', 'Perid': 'Perid', 'Zea': 'Zea',
             'DV_Chl_a': 'DVchla', 'MV_Chl_b': 'MVchlb', 'Chl_c1c2': 'Chlc12',
             'Chl_c3': 'Chlc3', 'Neo': 'Neo', 'Viola': 'Viola',
             'Tot_Chl_b': 'Tchlb', 'Tot_Chl_c': 'Tchlc', 'alpha-beta-Car': 'ABcaro',
             'Diadino': 'Diadino', 'Diato': 'Diato', 'MV_Chl_a': 'MVchla',
             'Chlide_a': 'Chllide', 'DV_Chl_b': 'DVchlb', 'Lut': 'Lut',
             'Phytin_a': 'Phytin', 'Phide_a': 'Phide', 'Pras': 'Pras'}


def km(lat1, lon1, lat2, lon2):
    """Equirectangular distance [km] (fine at these scales)."""
    return np.hypot((lat2 - lat1) * 111.2,
                    (lon2 - lon1) * 111.2 * np.cos(np.radians(0.5 * (lat1 + lat2))))


def load_hplc(index):
    frames = []
    for _, r in index[index['name'].str.contains('rosette_HPLC')].iterrows():
        _, df, bdl = seabass.read_sb(r['path'])
        for f in FIELD_MAP:
            if f in df.columns:
                df.loc[bdl[f], f] = 0.0
        df['ship'] = r['name'].split('_')[2]
        df['file'] = r['name']
        frames.append(df)
    H = pd.concat(frames, ignore_index=True)
    H = H[H['depth'] <= MAX_DEPTH]
    # shallowest bottle of each (ship, cast time), replicate-averaged
    keep = H.groupby(['ship', 'datetime'])['depth'].transform('min') == H['depth']
    H = H[keep]
    agg = {f: 'mean' for f in FIELD_MAP if f in H.columns}
    agg.update({'lat': 'first', 'lon': 'first', 'depth': 'first', 'file': 'first',
                'quality': lambda s: ';'.join(sorted(set(map(str, s))))})
    G = H.groupby(['ship', 'datetime']).agg(agg)
    G['n_replicates'] = H.groupby(['ship', 'datetime']).size()
    return G.reset_index().rename(columns=FIELD_MAP)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tf = sdpdata.kramer2022_dir() / 'ref' / 'Rrs_pigments' / 'Kramer_rrs_testdata.mat'
    t = loadmat(tf)
    lat, lon, chl = t['latlon'][:, 0], t['latlon'][:, 1], t['chl'].ravel()
    index = seabass.index_dir(SEABASS_DIR)
    G = load_hplc(index)

    rows, problems = [], []
    for i in range(len(chl)):
        d = km(lat[i], lon[i], G['lat'].to_numpy(), G['lon'].to_numpy())
        cand = G[(d <= MAX_KM) & (np.abs(G['Tchla'] - chl[i]) <= CHL_TOL)]
        if len(cand) != 1:
            problems.append({'i': i, 'n_candidates': int(len(cand)),
                             'nearest_km': float(np.min(d))})
            continue
        c = cand.iloc[0]
        rows.append({'test_index': i, 'test_lat': lat[i], 'test_lon': lon[i],
                     'test_chl': chl[i], 'd_km': float(d[cand.index[0]]),
                     'chl_diff': float(c['Tchla'] - chl[i]), **c.to_dict()})
    if problems:
        raise SystemExit(f'STOP: ambiguous or missing matchups: {problems}')
    M = pd.DataFrame(rows).sort_values('test_index')
    M['T'] = t['T'].ravel()
    M['S'] = t['S'].ravel()
    out = OUT_DIR / 'exports_na_matchups.csv'
    M.to_csv(out, index=False)

    used = index[index['name'].str.contains('rosette_HPLC')]
    prov = {'script': 'scripts/sdp/exports_na_matchups.py', 'epft_up_version': __version__,
            'test_spectra': {'file': str(tf), 'sha256': sdpdata.sha256_of(tf),
                             'repo': sdpdata.REF_REPOS['Rrs_pigments']},
            'seabass_dir': str(SEABASS_DIR),
            'seabass_files_used': used[['name', 'sha256', 'n_copies']].to_dict(orient='records'),
            'seabass_files_total': int(len(index)),
            'rule': {'max_km': MAX_KM, 'chl_tol': CHL_TOL, 'max_depth_m': MAX_DEPTH,
                     'bdl': 'set to 0 (Kramer)', 'replicates': 'mean'},
            'summary': {'n_matched': int(len(M)), 'max_d_km': float(M['d_km'].max()),
                        'max_abs_chl_diff': float(M['chl_diff'].abs().max()),
                        'ships': M['ship'].value_counts().to_dict(),
                        'depth_range_m': [float(M['depth'].min()), float(M['depth'].max())],
                        'quality_flags': M['quality'].value_counts().to_dict()},
            'output': {'file': str(out), 'sha256': sdpdata.sha256_of(out)},
            'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(out.with_suffix('.json'), 'w', encoding='utf-8') as fh:
        json.dump(prov, fh, indent=1, default=str)
    pd.set_option('display.width', 220)
    print(M[['test_index', 'ship', 'datetime', 'depth', 'n_replicates', 'd_km', 'test_chl',
             'Tchla', 'Fuco', 'HexFuco', 'Zea', 'DVchla', 'quality']].round(4).to_string())
    print(json.dumps(prov['summary'], indent=1, default=str))
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
