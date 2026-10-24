"""Test on new AMT data #1: harmonised AMT table, salinity fix and fold file.

Applies ``reports/provenance/AMT_benchmark_protocol.md`` §3–§4 through
:func:`epft_up.sdp.amt.load_amt` and writes, in ``$OS_COLOR/AMT/Rrs_HPLC/``:

* ``woa23_salinity_fix.csv`` (+ ``.json``): WOA23 monthly surface salinity
  (0.25°, objectively analysed, by OPeNDAP from NOAA NCEI) at the rows with
  implausible salinity. Fetched once and reused afterwards;
* ``AMT_Rrs_HPLC_harmonised.csv``: one row per spectrum with meta,
  pigments (our names), flags, the as-delivered Rrs, and Tchla from our
  GSM-like inversion and OC4v6, for both the 5 nm-smoothed spectra (the
  pipeline's primary input; MATLAB-style shrinking window at 400/700 nm) and
  the as-delivered spectra;
* ``AMT_Kramer_folds.csv``: ``sample_id, set, group, fold`` for the Kramer
  deposit (8 campaigns, folds 0-7) and AMT (8 cruises, folds 8-15);
* ``AMT_Rrs_HPLC_harmonised.json``: provenance, rows per cruise after each
  rule, GSM/OC4 diagnostics and output checksums.

Run with::

    conda run -n ocean14 python scripts/sdp/amt_harmonise.py
"""
import datetime
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from epft_up import __version__  # noqa: E402
from epft_up.sdp import amt, data as sdpdata, gsm, spectral  # noqa: E402
from epft_up.sdp import validate as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
D = amt.amt_dir()
SMOOTH = 5


def salinity_fix():
    """WOA23 surface S at each implausible-salinity row (cached on disk)."""
    dest = D / amt.SALINITY_FIX_FILE
    if dest.exists():
        print(f'reuse {dest}')
        return dest
    import xarray as xr
    from fetch_woa_ts import URL, nearest_valid
    raw = amt.load_amt(raw=True)
    bad = raw.meta['sal'] < amt.SAL_MIN_OPEN_OCEAN
    rows, urls = [], {}
    for i in np.where(bad)[0]:
        m = raw.meta.iloc[i]
        mm = int(m['time'].month)
        url = URL.format(var='salinity', v='s', mm=mm)
        urls[mm] = url
        print(f'GET {url}')
        ds = xr.open_dataset(url, decode_times=False)
        field = ds['s_an'].isel(time=0, depth=0).load().to_numpy()
        val, glat, glon, ring = nearest_valid(field, ds['lat'].to_numpy(),
                                              ds['lon'].to_numpy(), m['lat'], m['lon'])
        ds.close()
        rows.append({'Sample_ID': int(m['sample_id'].split(':')[1]), 'month': mm,
                     'lat': m['lat'], 'lon': m['lon'], 'S_insitu': m['sal'],
                     'S_woa': float(val), 'grid_lat': glat, 'grid_lon': glon, 'ring': ring})
    pd.DataFrame(rows).to_csv(dest, index=False)
    side = {'script': 'scripts/sdp/amt_harmonise.py',
            'source': 'NOAA NCEI WOA23 decav 0.25 deg monthly s_an, surface',
            'urls': urls, 'retrieved_utc': datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec='seconds'), 'output_sha256': sdpdata.sha256_of(dest)}
    dest.with_suffix('.json').write_text(json.dumps(side, indent=1))
    return dest


def gsm_oc4(d, Rrs, label):
    fit = gsm.fit_gsm(d.wave, Rrs, d.meta['temp'].to_numpy(), d.meta['sal'].to_numpy())
    oc4 = V.oc4v6(d.wave, Rrs)
    cols = {f'Tchla_GSM_{label}': fit.params[:, 0], f'adg443_GSM_{label}': fit.params[:, 1],
            f'bbp443_GSM_{label}': fit.params[:, 2], f'GSM_converged_{label}': fit.converged,
            f'Tchla_OC4_{label}': oc4}
    h = d.pigments['Tchla'].to_numpy()

    def st(y):
        ok = (y > 0) & np.isfinite(y)
        lx, ly = np.log10(h[ok]), np.log10(y[ok])
        b, a = np.polyfit(lx, ly, 1)
        return {'n': int(ok.sum()), 'slope': float(b), 'intercept': float(a),
                'R2': float(np.corrcoef(lx, ly)[0, 1]**2),
                'bias_dex': float(np.mean(ly - lx)), 'rms_dex': float(np.sqrt(np.mean((ly - lx)**2)))}
    diag = {'converged': int(fit.converged.sum()), 'n_bbp_negative': int((fit.params[:, 2] < 0).sum()),
            'GSM_vs_HPLC': st(fit.params[:, 0]), 'OC4_vs_HPLC': st(oc4)}
    return cols, diag


def folds(d):
    k = sdpdata.load_kramer2022()
    camp = k.meta['campaign'].to_numpy()
    kg = sorted(set(camp))
    ag = sorted(set(d.meta['cruise']))
    fk = {g: i for i, g in enumerate(kg)}
    fa = {g: len(kg) + i for i, g in enumerate(ag)}
    F = pd.concat([
        pd.DataFrame({'sample_id': [f'K22:{i}' for i in range(k.n)], 'set': 'kramer2022',
                      'group': camp, 'fold': [fk[c] for c in camp]}),
        pd.DataFrame({'sample_id': d.meta['sample_id'], 'set': 'amt',
                      'group': d.meta['cruise'], 'fold': d.meta['cruise'].map(fa)})],
        ignore_index=True)
    return F, k.provenance['sha256']


def main():
    fix = salinity_fix()
    d = amt.load_amt()
    print(d.steps.to_string())
    sm = spectral.moving_mean(d.Rrs, SMOOTH, edges='shrink')
    c1, g1 = gsm_oc4(d, sm, 'smooth5')
    c0, g0 = gsm_oc4(d, d.Rrs, 'asis')
    rcols = [f'Rrs_{int(w)}' for w in d.wave]
    meta = d.meta.copy()
    meta['time'] = meta['time'].dt.strftime('%Y-%m-%dT%H:%M:%SZ')
    H = pd.concat([meta, d.pigments, d.flags, pd.DataFrame(c1), pd.DataFrame(c0),
                   pd.DataFrame(d.Rrs, columns=rcols)], axis=1)
    out = D / 'AMT_Rrs_HPLC_harmonised.csv'
    H.to_csv(out, index=False)
    F, ksha = folds(d)
    fout = D / 'AMT_Kramer_folds.csv'
    F.to_csv(fout, index=False)

    git = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True,
                         text=True).stdout.strip()
    side = {'script': 'scripts/sdp/amt_harmonise.py', 'epft_up_version': __version__,
            'git_head': git, 'input': d.provenance,
            'salinity_fix': {'file': str(fix), 'sha256': sdpdata.sha256_of(fix),
                             'values': amt.read_salinity_fix(fix)},
            'kramer_deposit_sha256': ksha,
            'steps_rows_per_cruise': d.steps.to_dict(orient='index'),
            'n_spectra': d.n, 'n_hplc_averaged_gt1': int((d.meta['n_hplc'] > 1).sum()),
            'lod_zeroed_counts': {c.split(':')[1]: int(d.flags[c].sum())
                                  for c in d.flags if c.startswith('lod_zeroed:')},
            'mvchlb': {'proxy': int(d.flags['MVchlb_proxy'].sum()),
                       'scorable': int(d.flags['MVchlb_scorable'].sum())},
            'smoothing': {'width_nm': SMOOTH, 'edges': 'shrink (MATLAB movmean default)',
                          'applied_to': 'Tchla_*_smooth5 columns only; Rrs columns are as '
                                        'delivered'},
            'gsm_oc4': {'smooth5': g1, 'asis': g0},
            'folds': {'n_groups': int(F['group'].nunique()),
                      'kramer2022': F[F.set == 'kramer2022'].groupby('group')['fold'].first()
                      .to_dict(),
                      'amt': F[F.set == 'amt'].groupby('group')['fold'].first().to_dict(),
                      'schemes': {'L8': 'folds 0-7, kramer2022 only',
                                  'H2': 'train kramer2022, test set==amt (per group too)',
                                  'L16': 'all folds 0-15', 'LA': 'test = folds 8-15 together'}},
            'outputs': {p.name: {'file': str(p), 'sha256': sdpdata.sha256_of(p)}
                        for p in (out, fout)},
            'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    (D / 'AMT_Rrs_HPLC_harmonised.json').write_text(json.dumps(side, indent=1, default=str))
    print(json.dumps({k: side[k] for k in ('n_spectra', 'n_hplc_averaged_gt1',
                                           'lod_zeroed_counts', 'mvchlb', 'gsm_oc4')},
                     indent=1)[:3000])
    print(F.groupby(['set', 'group', 'fold']).size().to_string())
    print(f'wrote {out}\nwrote {fout}')


if __name__ == '__main__':
    main()
