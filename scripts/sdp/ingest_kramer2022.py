"""Execution #1: ingest checks and sanity figures for the Kramer 2022 data.

1. Load PANGAEA 937536 (checksum-verified) and compare the per-campaign
   counts and Tchla statistics with Kramer 2022 Table 1.
2. Apply the Lange et al. (2020) |Rrs''| screen at 610-660 nm (native 1 nm
   grid and Lange's 2 nm grid) and report how many of the 145 it rejects.
3. Test whether the deposited spectra are already smoothed by the paper's
   5 nm moving mean: a 5-point boxcar has transfer-function zeros at 0.2 and
   0.4 cycles/nm, so a pre-smoothed spectrum shows a power deficit there
   relative to 0.3 cycles/nm.
4. Reproduce Kramer 2022 Fig. 1 (map coloured by Tchla) and Fig. 2A (Rrs
   coloured by source).

Writes ``reports/figures/sdp/kramer_fig1_map.png``,
``reports/figures/sdp/kramer_fig2a_rrs.png`` and
``reports/figures/sdp/ingest_summary.json`` (numbers + provenance).

Run with::

    conda run -n ocean14 python scripts/sdp/ingest_kramer2022.py
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

from epft_up.sdp import data as sdpdata  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / 'reports' / 'figures' / 'sdp'

#: Source colours of Kramer 2022 Figs. 2, 4 and 6 (both Tara campaigns are "blue").
SOURCE_COLORS = {
    'ANT': 'red', 'NAAMES': 'orange', 'RemSensPOC': 'gold', 'SABOR': 'green',
    'Tara Oceans': 'blue', 'Tara Med': 'blue', 'BIOSOPE': 'purple',
    'EXPORTS': 'black'}


#: Rounding step of the deposited Rrs [sr^-1] (6 decimals; verified in main()).
QUANT = 1e-6
#: Frequencies [cycles/nm] at which the 2nd-difference power is reported. A
#: 5-point (5 nm) boxcar has transfer-function zeros at 0.2 and 0.4.
FREQS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)


def d2_power_vs_quantization(Rrs, quant=QUANT, freqs=FREQS, half=0.01):
    """Periodogram of the Rrs 2nd difference relative to the rounding floor.

    Rounding to ``quant`` adds (to a very good approximation) white noise of
    variance ``quant**2 / 12``; its 2nd difference has power
    ``(quant**2/12) * (2 sin(pi f))**4``. If a spectrum was smoothed by a
    5-point boxcar *before* rounding, its power at the boxcar zeros (0.2, 0.4
    cycles/nm) falls to this floor while the sidelobes (~0.3) stay above it.

    Parameters
    ----------
    Rrs : ndarray, shape (n_samples, n_wave)
        On a 1 nm grid.
    quant : float, optional
    freqs : sequence of float, optional
    half : float, optional
        Half-width of the frequency bin averaged around each ``f`` [cycles/nm].

    Returns
    -------
    dict
        ``{f: P_data / P_quantization}`` (sample-averaged, Hann-windowed).
    """
    d2 = np.diff(Rrs, n=2, axis=1)
    d2 = d2 - d2.mean(axis=1, keepdims=True)
    n = d2.shape[1]
    win = np.hanning(n)
    P = np.mean(np.abs(np.fft.rfft(d2 * win, axis=1))**2, axis=0)
    f = np.fft.rfftfreq(n, d=1.0)
    Pq = (quant**2 / 12.0) * (2.0 * np.sin(np.pi * f))**4 * np.sum(win**2)
    out = {}
    for f0 in freqs:
        sel = np.abs(f - f0) <= half
        out[f'{f0:.2f}'] = float(P[sel].mean() / Pq[sel].mean())
    return out


def fig1_map(data, path):
    """Kramer 2022 Fig. 1: sample locations coloured by Tchla."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(9, 4.8))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='0.85', edgecolor='none')
    ax.coastlines(lw=0.4)
    sc = ax.scatter(data.meta['lon'], data.meta['lat'], c=data.pigments['Tchla'],
                    s=22, cmap='viridis', norm=LogNorm(0.01, 5), edgecolor='k',
                    lw=0.3, transform=ccrs.PlateCarree(), zorder=3)
    cb = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label(r'HPLC Tchla [mg m$^{-3}$]')
    ax.set_title(f'Kramer et al. (2022) matchups, PANGAEA 937536 (N={data.n})')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig2a_rrs(data, path):
    """Kramer 2022 Fig. 2A: measured Rrs coloured by source."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for camp, color in SOURCE_COLORS.items():
        sel = (data.meta['campaign'] == camp).to_numpy()
        for r in data.Rrs[sel]:
            ax.plot(data.wave, r, color=color, lw=0.6, alpha=0.7)
    handles = [plt.Line2D([], [], color=c, label=k) for k, c in SOURCE_COLORS.items()
               if k != 'Tara Med']
    handles[4].set_label('Tara (Oceans + Med)')
    ax.legend(handles=handles, fontsize=8, frameon=False)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel(r'$R_{rs}$ [sr$^{-1}$]')
    ax.set_xlim(400, 700)
    ax.set_title(r'(A) measured $R_{rs}(\lambda)$, as deposited')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()

    # 1. Table 1
    t1 = sdpdata.table1_comparison(data)
    print(t1.round(3).to_string())
    chl = data.pigments['Tchla']
    overall = {'n': int(data.n), 'chl_min': float(chl.min()),
               'chl_max': float(chl.max()), 'chl_median': float(chl.median())}
    print('overall:', overall, '(paper: N=145, 0.019-4.15, median 0.110)')

    # 2. Lange QC
    qc = {}
    for label, step in (('native_1nm', None), ('lange_2nm', 2.0)):
        rej, mx = sdpdata.second_derivative_qc(data.wave, data.Rrs, step=step)
        qc[label] = {'n_reject': int(rej.sum()),
                     'max_abs_d2_quantiles': {q: float(np.quantile(mx, q))
                                              for q in (0.5, 0.9, 0.99, 1.0)},
                     'rejected_campaigns': data.meta.loc[rej, 'campaign']
                     .value_counts().to_dict()}
        print(f'Lange QC ({label}): rejects {rej.sum()} of {data.n}; '
              f'max|d2| median={np.median(mx):.2e}, max={mx.max():.2e}')

    # 3. pre-smoothing signature: 2nd-difference power vs the rounding floor
    on_grid = bool(np.allclose(np.round(data.Rrs / QUANT), data.Rrs / QUANT,
                               atol=1e-6))
    print(f'all Rrs are multiples of {QUANT:g} sr^-1: {on_grid}')
    smooth = {'rounding_step_sr-1': QUANT, 'all_on_rounding_grid': on_grid,
              'P_data_over_P_quant': {'all': d2_power_vs_quantization(data.Rrs)}}
    for camp in sdpdata.TABLE1.index:
        sel = (data.meta['campaign'] == camp).to_numpy()
        smooth['P_data_over_P_quant'][camp] = d2_power_vs_quantization(data.Rrs[sel])
    print('P_data/P_quant of the 2nd difference, by frequency [cycles/nm]:')
    print('  ' + ' '.join(f'{f:>6.2f}' for f in FREQS))
    for k, v in smooth['P_data_over_P_quant'].items():
        print(f'  {" ".join(f"{x:6.2f}" if x < 1e3 else f"{x:6.0e}" for x in v.values())}  {k}')

    # 4. figures
    fig1_map(data, OUT / 'kramer_fig1_map.png')
    fig2a_rrs(data, OUT / 'kramer_fig2a_rrs.png')

    summary = {
        'script': 'scripts/sdp/ingest_kramer2022.py',
        'provenance': data.provenance,
        'table1': json.loads(t1.to_json(orient='index')),
        'overall': overall,
        'lange_qc': qc,
        'presmoothing_signature': smooth,
        'rrs_min': float(data.Rrs.min()), 'rrs_max': float(data.Rrs.max()),
        'n_rrs_nonpositive': int((data.Rrs <= 0).sum()),
        'pigment_zero_fraction': (data.pigments == 0).mean().round(3).to_dict(),
        'figures': ['reports/figures/sdp/kramer_fig1_map.png',
                    'reports/figures/sdp/kramer_fig2a_rrs.png'],
    }
    with open(OUT / 'ingest_summary.json', 'w') as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f'wrote {OUT / "ingest_summary.json"}')


if __name__ == '__main__':
    main()
