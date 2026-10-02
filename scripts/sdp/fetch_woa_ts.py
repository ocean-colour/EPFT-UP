"""Extract WOA climatological surface T and S at the Kramer 2022 samples (Q&A #22).

Kramer 2022 computed b_bw with temperature and salinity from the World Ocean
Atlas 2013 1/4-degree climatology; the PANGAEA deposit instead carries
in-situ T/S, which we use by default. To quantify the difference, this script
pulls the **WOA23** (the current release; WOA13 is no longer served) 1/4°
*monthly* "decav" objectively analysed fields (``t_an``, ``s_an``) at the
surface level, for each sample's calendar month, by OPeNDAP from NOAA NCEI.
The nearest ocean grid point with data is used (searching outward up to
``MAX_RING`` cells when the nearest cell is land or missing).

Output: ``$OS_COLOR/PANGAEA/Kramer2022/woa23_surface_ts.csv`` with one row per
sample (same order as the deposit) and a JSON sidecar with provenance.

Run with::

    conda run -n ocean14 python scripts/sdp/fetch_woa_ts.py
"""
import datetime
import json
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))

from epft_up.sdp import data as sdpdata  # noqa: E402

URL = ('https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/{var}/netcdf/'
       'decav/0.25/woa23_decav_{v}{mm:02d}_04.nc')
VARS = {'temperature': 't', 'salinity': 's'}
#: Outward search radius [grid cells] when the nearest cell has no data.
MAX_RING = 4


def nearest_valid(field, lats, lons, lat0, lon0):
    """Value of the nearest non-NaN cell to (lat0, lon0) and its location."""
    i0 = int(np.argmin(np.abs(lats - lat0)))
    j0 = int(np.argmin(np.abs(lons - lon0)))
    for ring in range(MAX_RING + 1):
        best = None
        for di in range(-ring, ring + 1):
            for dj in range(-ring, ring + 1):
                if max(abs(di), abs(dj)) != ring:
                    continue
                i, j = i0 + di, (j0 + dj) % len(lons)
                if 0 <= i < len(lats) and np.isfinite(field[i, j]):
                    d2 = (lats[i] - lat0)**2 + ((lons[j] - lon0) * np.cos(np.deg2rad(lat0)))**2
                    if best is None or d2 < best[0]:
                        best = (d2, field[i, j], lats[i], lons[j])
        if best is not None:
            return best[1], best[2], best[3], ring
    return np.nan, np.nan, np.nan, -1


def main():
    data = sdpdata.load_kramer2022()
    meta = data.meta
    months = meta['time'].dt.month.to_numpy()
    out = pd.DataFrame({'month': months, 'lat': meta['lat'], 'lon': meta['lon']})
    urls = {}
    for var, v in VARS.items():
        vals = np.full(data.n, np.nan)
        glat = np.full(data.n, np.nan)
        glon = np.full(data.n, np.nan)
        ring = np.full(data.n, -1)
        for mm in sorted(set(months)):
            url = URL.format(var=var, v=v, mm=mm)
            urls[f'{var}_{mm:02d}'] = url
            print(f'GET {url}')
            ds = xr.open_dataset(url, decode_times=False)
            field = ds[f'{v}_an'].isel(time=0, depth=0).load().to_numpy()
            lats, lons = ds['lat'].to_numpy(), ds['lon'].to_numpy()
            for i in np.where(months == mm)[0]:
                vals[i], glat[i], glon[i], ring[i] = nearest_valid(
                    field, lats, lons, meta['lat'].iloc[i], meta['lon'].iloc[i])
            ds.close()
        out[f'{v.upper()}_woa'] = vals
        out[f'{v}_grid_lat'] = glat
        out[f'{v}_grid_lon'] = glon
        out[f'{v}_ring'] = ring

    dest = sdpdata.kramer2022_dir() / 'woa23_surface_ts.csv'
    out.to_csv(dest, index=False)
    side = {'source': 'NOAA NCEI World Ocean Atlas 2023, decav, 0.25 deg, monthly, '
                      'objectively analysed (t_an, s_an), surface level',
            'urls': urls, 'max_ring': MAX_RING,
            'n_missing_T': int(np.isnan(out['T_woa']).sum()),
            'n_missing_S': int(np.isnan(out['S_woa']).sum()),
            'retrieved_utc': datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec='seconds'),
            'kramer_input_sha256': data.provenance['sha256'],
            'output_sha256': sdpdata.sha256_of(dest)}
    with open(dest.with_suffix('.json'), 'w', encoding='utf-8') as fh:
        json.dump(side, fh, indent=1)
    print(out.describe().round(3).to_string())
    print(f'wrote {dest}')


if __name__ == '__main__':
    main()
