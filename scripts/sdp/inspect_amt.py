"""New AMT data #1: inspect Tom Jordan's ``AMT_Rrs_HPLC.csv`` for the SDP work.

Question: can this compilation of Atlantic Meridional Transect (AMT)
above-water Rrs + HPLC matchups **test** our frozen models (as EXPORTS-NA did
in Execution #9) and/or **train** new ones alongside Kramer's 145?

The script reports (no model is applied or trained here):

1. **Inventory:** samples per cruise, dates, HPLC lab, radiometry method,
   positions, T/S, duplicates.
2. **Pigments:** mapping to our 13 (``But`` → ButFuco, ``Hex`` → HexFuco),
   missing columns, zero (below-detection) fractions, internal consistency
   (Tchla vs MVchla + DVchla), ranges against the Kramer deposit.
3. **Radiometry:** band grid, NaNs, non-positive values, the deposit's
   rounding/pre-smoothing signature (§3.4 of the report), and the Lange QC
   screen (§3.5).
4. **Feasibility of the SDP chain:** the GSM-like fit (Execution #2 code and
   tables, in-situ T/S) and OC4v6 against HPLC Tchla, compared with the
   same statistics on the deposit (Fig. 4). This tests whether radiometry and
   pigments are mutually consistent at the level the models need.
5. **Overlap with the deposit:** nearest Kramer sample in space and time.

Outputs: ``reports/figures/sdp/amt_inspect_summary.json``,
``amt_map.png``, ``amt_spectra.png``, ``amt_pigments.png`` and
``amt_chl_consistency.png``.

Run with::

    conda run -n ocean14 python scripts/sdp/inspect_amt.py [path/to/AMT_Rrs_HPLC.csv]
"""
import datetime
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import data as sdpdata, gsm, spectral  # noqa: E402
from epft_up.sdp import validate as V  # noqa: E402
from ingest_kramer2022 import d2_power_vs_quantization  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
DEFAULT = Path.home() / 'Downloads' / 'AMT_Rrs_HPLC.csv'

#: AMT column -> our name (data.PIGMENTS_13 / PIGMENTS_EXTRA where they exist)
PIG_MAP = {'Tchla': 'Tchla', 'Hex': 'HexFuco', 'But': 'ButFuco', 'Allo': 'Allo',
           'Fuco': 'Fuco', 'Perid': 'Perid', 'Zea': 'Zea', 'DVchla': 'DVchla',
           'MVchlb': 'MVchlb', 'Chlc12': 'Chlc12', 'Chlc3': 'Chlc3', 'Neo': 'Neo',
           'Viola': 'Viola', 'Tchlb': 'Tchlb', 'Tchlc': 'Tchlc', 'ABcaro': 'ABcaro',
           'Diadino': 'Diadino', 'Diato': 'Diato', 'MVchla': 'MVchla', 'DVchlb': 'DVchlb',
           'Lut': 'Lut', 'Pras': 'Pras', 'Psc': 'Psc', 'Ppc': 'Ppc'}


def log_fit(x, y):
    """OLS of log10 y on log10 x: slope, intercept, R², bias, RMS (dex)."""
    ok = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    lx, ly = np.log10(x[ok]), np.log10(y[ok])
    b, a = np.polyfit(lx, ly, 1)
    return {'n': int(ok.sum()), 'slope': float(b), 'intercept': float(a),
            'R2': float(np.corrcoef(lx, ly)[0, 1]**2), 'bias_dex': float(np.mean(ly - lx)),
            'rms_dex': float(np.sqrt(np.mean((ly - lx)**2)))}


def km(lat1, lon1, lat2, lon2):
    return np.hypot((lat2 - lat1) * 111.2,
                    (lon2 - lon1) * 111.2 * np.cos(np.radians(0.5 * (lat1 + lat2))))


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(path, na_values=['nan', 'NaN', ''])
    rcols = [c for c in df.columns if c.startswith('Rrs_')]
    wave = np.array([float(c[4:]) for c in rcols])
    meta_cols = [c for c in df.columns if c not in rcols and c not in PIG_MAP]
    Rrs = df[rcols].to_numpy(float)
    # MM/DD/YYYY: AMT23 runs south from the UK in Oct 2013 (10/09 → 10/14 going south)
    df['time'] = pd.to_datetime(df['DateTime'], format='%m/%d/%Y %H:%M', errors='coerce')
    S = {'script': 'scripts/sdp/inspect_amt.py', 'epft_up_version': __version__,
         'input': {'file': str(path), 'sha256': sdpdata.sha256_of(path),
                   'bytes': path.stat().st_size, 'source': 'Tom Jordan (shared 2026-10-24)'},
         'n_rows': int(len(df)), 'columns': {'meta': meta_cols,
                                             'pigments': [c for c in df.columns if c in PIG_MAP],
                                             'n_rrs': len(rcols)}}

    # 1. inventory
    S['unparsed_dates'] = int(df['time'].isna().sum())
    g = df.groupby('Cruise')
    inv = pd.DataFrame({'n': g.size(), 'start': g['time'].min().dt.date.astype(str),
                        'end': g['time'].max().dt.date.astype(str),
                        'lat_min': g['Lat'].min().round(1), 'lat_max': g['Lat'].max().round(1),
                        'HPLC_lab': g['HPLC_lab'].agg(lambda s: '/'.join(sorted(set(s)))),
                        'Rad_method': g['Rad_method'].agg(lambda s: '/'.join(sorted(set(s)))),
                        'Tchla_median': g['Tchla'].median().round(3)})
    S['by_cruise'] = inv.reset_index().to_dict(orient='records')
    S['data_links'] = sorted(df['Data_link'].dropna().unique().tolist())
    S['sample_id_unique'] = bool(df['Sample_ID'].is_unique)
    S['original_id_unique'] = bool(df['Original_ID'].is_unique)
    S['duplicate_spectra'] = int(pd.DataFrame(Rrs).round(10).duplicated().sum())
    S['duplicate_time_position'] = int(df[['time', 'Lat', 'Lon']].duplicated().sum())
    # duplicates: which rows share a spectrum / an Original_ID / a time+position
    key = pd.Series([hash(r.tobytes()) for r in np.round(Rrs, 10)], index=df.index)
    dup_groups = []
    for _, idx in key.groupby(key).groups.items():
        if len(idx) > 1:
            sub = df.loc[idx]
            dup_groups.append({'rows': sub['Sample_ID'].tolist(),
                               'original_id': sub['Original_ID'].tolist(),
                               'cruise': sorted(set(sub['Cruise'])),
                               'time': sub['DateTime'].tolist(),
                               'Tchla': sub['Tchla'].round(4).tolist(),
                               'same_pigments': bool(sub[[c for c in PIG_MAP if c in df]]
                                                     .round(8).nunique().max() == 1)})
    S['duplicate_spectrum_groups'] = dup_groups
    S['n_unique_spectra'] = int(key.nunique())
    S['dup_group_sizes'] = sorted({len(g['rows']) for g in dup_groups})
    pos = (Rrs > 0).all(axis=1)
    S['n_unique_positive_spectra'] = int(key[pos].nunique())
    oid = df[df['Original_ID'].duplicated(keep=False)]
    S['duplicate_original_id'] = oid[['Sample_ID', 'Original_ID', 'Cruise', 'DateTime',
                                      'Tchla']].to_dict(orient='records')
    lows = df[df['Sal'] < 33]
    S['low_salinity_rows'] = lows[['Sample_ID', 'Original_ID', 'Cruise', 'DateTime', 'Lat',
                                   'Lon', 'Temp', 'Sal', 'Tchla']].round(3).to_dict(
                                       orient='records')
    S['T_range'] = [float(df['Temp'].min()), float(df['Temp'].max())]
    S['S_range'] = [float(df['Sal'].min()), float(df['Sal'].max())]
    S['T_S_missing'] = int(df[['Temp', 'Sal']].isna().any(axis=1).sum())

    # 2. pigments
    P = df[[c for c in PIG_MAP if c in df.columns]].rename(columns=PIG_MAP)
    p13 = list(sdpdata.PIGMENTS_13)
    S['pigments13_present'] = [p for p in p13 if p in P.columns]
    S['pigment_missing_fraction'] = {p: float(P[p].isna().mean()) for p in P.columns}
    S['pigment_missing_by_cruise'] = {
        p: df.groupby('Cruise')[k].apply(lambda s: float(s.isna().mean())).to_dict()
        for k, p in PIG_MAP.items() if k in df.columns and P[p].isna().any()}
    S['pigment_zero_fraction'] = {p: float((P[p] == 0).mean()) for p in P.columns}
    S['pigment_negative'] = {p: int((P[p] < 0).sum()) for p in P.columns if (P[p] < 0).any()}
    # MVchlb is missing: Tchlb = MVchlb + DVchlb, and DVchlb is missing too
    S['Tchlb_as_MVchlb_proxy'] = {'Tchlb_missing': float(P['Tchlb'].isna().mean())}
    tsum = P['MVchla'] + P['DVchla']
    rel = (P['Tchla'] - tsum) / P['Tchla']
    S['Tchla_minus_MV_DV_rel'] = {'median': float(np.nanmedian(rel)),
                                  'p05_p95': [float(np.nanpercentile(rel, 5)),
                                              float(np.nanpercentile(rel, 95))]}
    # reporting resolution: smallest non-zero value per pigment
    S['pigment_min_nonzero'] = {p: float(P[p][P[p] > 0].min()) for p in p13
                                if p in P and (P[p] > 0).any()}
    kd = sdpdata.load_kramer2022()
    kp = kd.pigments
    S['Tchla_vs_deposit'] = {
        'amt': {'min': float(P['Tchla'].min()), 'median': float(P['Tchla'].median()),
                'max': float(P['Tchla'].max()),
                'frac_lt_0.1': float((P['Tchla'] < 0.1).mean())},
        'deposit': {'min': float(kp['Tchla'].min()), 'median': float(kp['Tchla'].median()),
                    'max': float(kp['Tchla'].max()),
                    'frac_lt_0.1': float((kp['Tchla'] < 0.1).mean())}}

    # 3. radiometry
    S['wave'] = {'min': float(wave.min()), 'max': float(wave.max()),
                 'uniform_1nm': bool(np.allclose(np.diff(wave), 1.0)),
                 'same_as_deposit': bool(np.array_equal(wave, kd.wave))}
    S['rrs_nan_rows'] = int(np.isnan(Rrs).any(axis=1).sum())
    S['rrs_nonpositive_values'] = int(np.sum(Rrs <= 0))
    S['rrs_rows_with_nonpositive'] = int(np.sum((Rrs <= 0).any(axis=1)))
    S['rrs_nonpositive_by_band'] = {
        f'{int(wave[i])}': int(c) for i, c in enumerate(np.sum(Rrs <= 0, axis=0)) if c}
    good = ~np.isnan(Rrs).any(axis=1)
    S['rrs_range'] = [float(np.nanmin(Rrs)), float(np.nanmax(Rrs))]
    S['on_1e-6_grid'] = bool(np.allclose(np.round(Rrs[good] * 1e6), Rrs[good] * 1e6,
                                         atol=1e-6))
    S['d2_power_over_1e-6_rounding_floor'] = {
        'amt': {k: float(v) for k, v in d2_power_vs_quantization(Rrs[good]).items()},
        'deposit': {k: float(v) for k, v in d2_power_vs_quantization(kd.Rrs).items()}}
    S['d2_power_by_cruise'] = {
        c: {k: float(v) for k, v in
            d2_power_vs_quantization(Rrs[good & (df['Cruise'] == c).to_numpy()]).items()}
        for c in sorted(df['Cruise'].unique())}
    # after the paper's 5 nm moving mean (as Kramer applied to raw spectra)
    sm = spectral.moving_mean(Rrs[good], 5)[:, 2:-2]  # drop NaN edges
    S['d2_power_after_5nm_mean'] = {k: float(v) for k, v in
                                    d2_power_vs_quantization(sm).items()}
    rej, mx = sdpdata.second_derivative_qc(wave, np.nan_to_num(Rrs), step=2.0)
    S['lange_qc_2nm'] = {'n_reject': int(rej[good].sum()),
                         'max_abs_d2_median': float(np.median(mx[good])),
                         'by_cruise': df[good].assign(r=rej[good]).groupby('Cruise')['r']
                         .sum().astype(int).to_dict()}

    # 4. SDP chain feasibility (GSM + OC4 vs HPLC), same code as Execution #2
    ok = good & df[['Temp', 'Sal']].notna().all(axis=1).to_numpy() & (Rrs > 0).all(axis=1)
    S['n_fit_eligible'] = int(ok.sum())
    tables = gsm.load_ref_tables(wave)
    fit = gsm.fit_gsm(wave, Rrs[ok], df.loc[ok, 'Temp'].to_numpy(),
                      df.loc[ok, 'Sal'].to_numpy(), tables=tables)
    chl_gsm = fit.params[:, 0]
    oc4 = V.oc4v6(wave, Rrs[ok])
    hplc = P.loc[ok, 'Tchla'].to_numpy()
    S['gsm'] = {'converged': int(fit.converged.sum()), 'n': int(ok.sum()),
                'n_bbp_negative': int(np.sum(fit.params[:, 2] < 0)),
                'GSM_vs_HPLC': log_fit(hplc, chl_gsm), 'OC4_vs_HPLC': log_fit(hplc, oc4)}
    cr = df.loc[ok, 'Cruise'].to_numpy()
    S['gsm']['by_cruise'] = {c: {'GSM': log_fit(hplc[cr == c], chl_gsm[cr == c]),
                                 'OC4': log_fit(hplc[cr == c], oc4[cr == c])}
                             for c in sorted(set(cr)) if (cr == c).sum() >= 5}
    gz = json.loads((OUT / 'gsm_summary.json').read_text())['fig4']
    S['deposit_fig4'] = {'OC4': gz['OC4v6'], 'GSM_without_worst': gz['GSM_without_worst']}
    dr = fit.dRrs
    S['dRrs_rms'] = {'amt': float(np.sqrt(np.nanmean(dr**2))),
                     'deposit': float(np.sqrt(np.mean(np.load(
                         sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz')
                         ['dRrs']**2)))}

    # 5. overlap with the deposit
    km_min = np.array([np.min(km(la, lo, kd.meta['lat'].to_numpy(), kd.meta['lon'].to_numpy()))
                       for la, lo in zip(df['Lat'], df['Lon'])])
    S['nearest_deposit_km'] = {'min': float(km_min.min()), 'median': float(np.median(km_min)),
                               'n_within_10km': int((km_min < 10).sum())}
    S['run_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')

    # figures
    cruises = sorted(df['Cruise'].unique())
    cmap = plt.get_cmap('tab10')
    col = {c: cmap(i % 10) for i, c in enumerate(cruises)}
    try:
        import cartopy.crs as ccrs
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines(lw=0.5)
        ax.set_extent([-70, 20, -60, 60])
        tr = dict(transform=ccrs.PlateCarree())
    except ImportError:  # pragma: no cover
        fig, ax = plt.subplots(figsize=(8, 6))
        tr = {}
    for c in cruises:
        m = df['Cruise'] == c
        ax.scatter(df.loc[m, 'Lon'], df.loc[m, 'Lat'], s=14, color=col[c], label=c, **tr)
    ax.scatter(kd.meta['lon'], kd.meta['lat'], s=10, marker='x', color='k',
               label='Kramer 2022 (145)', **tr)
    ax.legend(fontsize=7, loc='lower left', ncol=2)
    ax.set_title('AMT_Rrs_HPLC.csv: sample positions')
    fig.savefig(OUT / 'amt_map.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for c in cruises:
        m = (df['Cruise'] == c).to_numpy() & good
        axs[0].plot(wave, Rrs[m].T, color=col[c], lw=0.5, alpha=0.6)
        axs[0].plot([], [], color=col[c], label=c)
    axs[0].set_xlabel('wavelength [nm]')
    axs[0].set_ylabel(r'Rrs [sr$^{-1}$]')
    axs[0].legend(fontsize=7)
    axs[0].set_title('AMT Rrs by cruise')
    fq = sorted(S['d2_power_over_1e-6_rounding_floor']['amt'], key=float)
    axs[1].semilogy([float(f) for f in fq],
                    [S['d2_power_over_1e-6_rounding_floor']['deposit'][f] for f in fq],
                    'k-o', label='Kramer deposit')
    for c in cruises:
        d = S['d2_power_by_cruise'][c]
        axs[1].semilogy([float(f) for f in fq], [d[f] for f in fq], '-o', ms=3, color=col[c],
                        label=c)
    axs[1].axhline(1, color='0.5', ls=':')
    axs[1].set_xlabel('frequency [cycles nm$^{-1}$]')
    axs[1].set_ylabel("power of Rrs'' / 1e-6 rounding floor")
    axs[1].set_title("High-frequency content of Rrs''")
    axs[1].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / 'amt_spectra.png', dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    bins = np.linspace(-2.2, 1, 33)
    axs[0].hist(np.log10(kp['Tchla']), bins, alpha=0.6, label='Kramer deposit (145)')
    axs[0].hist(np.log10(P['Tchla'][P['Tchla'] > 0]), bins, alpha=0.6,
                label=f'AMT ({len(df)})')
    axs[0].set_xlabel('log10 HPLC Tchla [mg m$^{-3}$]')
    axs[0].legend()
    pz = [p for p in p13 if p in P]
    x = np.arange(len(pz))
    axs[1].bar(x - 0.2, [S['pigment_zero_fraction'][p] for p in pz], 0.4, label='AMT zero')
    axs[1].bar(x + 0.2, [float((kp[p] == 0).mean()) for p in pz], 0.4, label='deposit zero')
    axs[1].bar(x - 0.2, [S['pigment_missing_fraction'][p] for p in pz], 0.4,
               bottom=[S['pigment_zero_fraction'][p] for p in pz], color='k',
               label='AMT missing')
    axs[1].set_xticks(x, pz, rotation=45, ha='right')
    axs[1].set_ylabel('fraction of samples')
    axs[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'amt_pigments.png', dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    for ax, y, lab in ((axs[0], oc4, 'OC4v6'), (axs[1], chl_gsm, 'GSM-like')):
        for c in sorted(set(cr)):
            m = cr == c
            ax.loglog(hplc[m], np.clip(y[m], 1e-3, None), 'o', ms=4, color=col[c], label=c)
        st = log_fit(hplc, y)
        ax.plot([1e-3, 10], [1e-3, 10], 'k-', lw=0.6)
        ax.set_title(f"{lab}: slope {st['slope']:.2f}, R² {st['R2']:.2f}, "
                     f"bias {st['bias_dex']:+.2f} dex")
        ax.set_xlabel('HPLC Tchla [mg m$^{-3}$]')
    axs[0].set_ylabel('retrieved Tchla [mg m$^{-3}$]')
    axs[0].legend(fontsize=7)
    axs[0].set_xlim(5e-3, 10)
    axs[0].set_ylim(1e-3, 10)
    fig.tight_layout()
    fig.savefig(OUT / 'amt_chl_consistency.png', dpi=150)
    plt.close(fig)

    with open(OUT / 'amt_inspect_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(S, fh, indent=1, default=str)
    pd.set_option('display.width', 200)
    print(inv.to_string())
    for k in ('n_rows', 'unparsed_dates', 'sample_id_unique', 'original_id_unique',
              'duplicate_spectra', 'duplicate_time_position', 'T_S_missing', 'T_range',
              'S_range', 'pigments13_present', 'Tchla_minus_MV_DV_rel', 'Tchla_vs_deposit',
              'wave', 'rrs_nan_rows', 'rrs_nonpositive_values', 'rrs_rows_with_nonpositive',
              'rrs_nonpositive_by_band', 'rrs_range', 'on_1e-6_grid', 'lange_qc_2nm',
              'n_fit_eligible', 'gsm', 'deposit_fig4', 'dRrs_rms', 'nearest_deposit_km',
              'data_links', 'duplicate_spectrum_groups', 'duplicate_original_id',
              'low_salinity_rows', 'd2_power_after_5nm_mean', 'n_unique_spectra',
              'dup_group_sizes', 'n_unique_positive_spectra'):
        print(k, json.dumps(S[k], default=str)[:2500])
    print('missing:', {k: round(v, 2) for k, v in S['pigment_missing_fraction'].items() if v})
    print('missing by cruise:', json.dumps(S['pigment_missing_by_cruise'])[:900])
    print('zero:', {k: round(v, 2) for k, v in S['pigment_zero_fraction'].items()})
    print('min nonzero:', S['pigment_min_nonzero'])
    print('d2 power amt:', {k: round(v, 2) for k, v in
                            S['d2_power_over_1e-6_rounding_floor']['amt'].items()})
    for c, d in S['d2_power_by_cruise'].items():
        print('  ', c, {k: round(v, 2) for k, v in d.items()})


if __name__ == '__main__':
    main()
