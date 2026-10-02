"""Execution #2: reproduce the Kramer 2022 GSM-like fit and the residual δRrs.

Steps
-----
1. Fit the 145 deposited spectra with ``epft_up.sdp.gsm.fit_gsm`` (Nelder-Mead
   as in ``gsm_invert.m``; in-situ T/S for b_bw). No extra smoothing: the
   deposit already carries the paper's 5 nm moving mean (report §3.4).
2. Fig. 4A/4B: OC4v6 and GSM Tchla vs HPLC (log10 statistics, as the paper's
   OC4 numbers confirm), with and without the one degenerate fit.
3. Fig. 2B/2C: modelled Rrs and δRrs coloured by source.
4. Cross-checks against the references: Zhang (2009) b_sw against the Python
   port's ``betasw124_ZHH2009``; δRrs against the port's
   ``get_rrs_residuals``.
5. Sensitivities of δRrs: (a) WOA23 vs in-situ T/S (Q&A #22; needs
   ``scripts/sdp/fetch_woa_ts.py``); (b) the η band-ratio convention; (c) the
   optimizer (bounded log-parameter least squares); (d) a second, accidental
   5 nm moving mean.
6. Save the δRrs product (with provenance) for Execution #3 to
   ``$OS_COLOR/PANGAEA/Kramer2022/products/gsm_dRrs_insitu.npz``.

Outputs: ``reports/figures/sdp/kramer_fig4_chl.png``,
``reports/figures/sdp/kramer_fig2bc_model_residual.png``,
``reports/figures/sdp/gsm_sensitivity.png`` and
``reports/figures/sdp/gsm_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/reproduce_gsm.py
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

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import data as sdpdata, gsm, spectral  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
SOURCE_COLORS = {
    'ANT': 'red', 'NAAMES': 'orange', 'RemSensPOC': 'gold', 'SABOR': 'green',
    'Tara Oceans': 'blue', 'Tara Med': 'blue', 'BIOSOPE': 'purple',
    'EXPORTS': 'black'}
#: OC4v6 (O'Reilly et al. 1998, SeaWiFS v6) coefficients; bands 443/490/510 over 555.
OC4V6 = (0.3272, -2.9940, 2.7218, -1.2259, -0.5683)


def oc4v6(wave, Rrs):
    """OC4v6 Tchla [mg m^-3] from hyperspectral Rrs sampled at 443/490/510/555 nm."""
    r = lambda w: Rrs[:, int(np.searchsorted(wave, w))]  # noqa: E731
    x = np.log10(np.maximum.reduce([r(443), r(490), r(510)]) / r(555))
    return 10.0**np.polyval(OC4V6[::-1], x)


def log_stats(truth, est):
    """Kramer Fig. 4 statistics in log10 space: R², OLS slope and intercept, n."""
    ok = np.isfinite(est) & (est > 0)
    x, y = np.log10(truth[ok]), np.log10(est[ok])
    slope, icpt = np.polyfit(x, y, 1)
    return {'n': int(ok.sum()), 'R2': float(np.corrcoef(x, y)[0, 1]**2),
            'slope': float(slope), 'intercept': float(icpt),
            'median_ratio': float(np.median(10**(y - x))),
            'rms_log10': float(np.sqrt(np.mean((y - x)**2)))}


def count_nonpositive_trials(fit, i):
    """Re-run sample ``i``'s Nelder-Mead and count trial points with Tchla <= 0."""
    sh = fit.shapes
    args = (sh.aw, sh.A, sh.B, sh.bbw[i], sh.adg_star[i], sh.bbp_star[i])
    seen = {'n': 0, 'neg': 0}

    def cost(p):
        seen['n'] += 1
        if p[0] <= 0:
            seen['neg'] += 1
        return gsm._cost(p, fit.rrs[i], *args)

    from scipy.optimize import minimize
    minimize(cost, np.asarray(gsm.X0), method='Nelder-Mead',
             options={'xatol': 1e-9, 'fatol': 1e-9, 'maxiter': 2000, 'maxfev': 2000})
    return seen


def run_port(data):
    """δRrs from the Python port (``max-danenhower/rrs-SDP-pigments``)."""
    src = sdpdata.kramer2022_dir() / 'ref' / 'rrs-SDP-pigments' / 'src'
    sys.path.insert(0, str(src))
    import warnings
    from sdp import Kramer_hyperRrs as port  # pylint: disable=import-error
    cols = [f'{int(w)}' for w in data.wave]
    Rrs = pd.DataFrame(data.Rrs, columns=cols)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, RrsD = port.get_rrs_residuals(Rrs, data.meta['temp'].to_numpy(),
                                         data.meta['sal'].to_numpy(), data.wave)
        _, bsw_port, _, _ = port.betasw124_ZHH2009(
            data.wave[None, :], data.meta['sal'].to_numpy()[:, None],
            data.meta['temp'].to_numpy()[:, None])
    return np.asarray(RrsD).T, np.asarray(bsw_port)


def rel_rms(diff, ref):
    """RMS of ``diff`` relative to the RMS of ``ref`` (over all finite values)."""
    m = np.isfinite(diff) & np.isfinite(ref)
    return float(np.sqrt(np.mean(diff[m]**2)) / np.sqrt(np.mean(ref[m]**2)))


def fig4(data, chl_oc4, fit, worst, path):
    hplc = data.pigments['Tchla'].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharex=True, sharey=True)
    for ax, est, title, is_gsm in ((axes[0], chl_oc4, '(A) OC4v6', False),
                                   (axes[1], fit.params[:, 0], '(B) GSM-like model', True)):
        for camp, color in SOURCE_COLORS.items():
            sel = (data.meta['campaign'] == camp).to_numpy()
            ax.scatter(hplc[sel], est[sel], s=16, color=color, edgecolor='k', lw=0.3,
                       label=camp if camp != 'Tara Med' else None)
        ax.plot([1e-3, 10], [1e-3, 10], 'k-', lw=0.8)
        s = log_stats(hplc, est)
        txt = f'R² = {s["R2"]:.2f}, slope = {s["slope"]:.2f} (N={s["n"]})'
        if is_gsm:
            keep = np.ones(data.n, bool)
            keep[worst] = False
            s2 = log_stats(hplc[keep], est[keep])
            ax.scatter(hplc[worst], est[worst], s=80,
                       facecolor='none', edgecolor='m', lw=1.2)
            txt += (f'\nwithout the circled fit: R² = {s2["R2"]:.2f}, '
                    f'slope = {s2["slope"]:.2f}')
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=8.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-2, 6)
        ax.set_ylim(2e-4, 6)
        ax.set_xlabel(r'HPLC Tchla [mg m$^{-3}$]')
        ax.set_title(title)
    axes[0].set_ylabel(r'modelled Tchla [mg m$^{-3}$]')
    axes[0].legend(fontsize=7, frameon=False, loc='lower right')
    fig.suptitle('Kramer 2022 Fig. 4 (paper: OC4 R²=0.75, slope 0.87; GSM R²=0.86, slope 0.96)',
                 fontsize=9.5)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig2bc(data, fit, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for camp, color in SOURCE_COLORS.items():
        sel = (data.meta['campaign'] == camp).to_numpy()
        for r_mod, d in zip(fit.Rrs_mod[sel], fit.dRrs[sel]):
            axes[0].plot(fit.wave, r_mod, color=color, lw=0.6, alpha=0.7)
            axes[1].plot(fit.wave, d, color=color, lw=0.6, alpha=0.7)
    axes[0].set_title(r'(B) modelled $R_{rs}(\lambda)$')
    axes[0].set_ylabel(r'$R_{rs}$ [sr$^{-1}$]')
    axes[1].set_title(r'(C) residual $\delta R_{rs}(\lambda)$')
    axes[1].set_ylabel(r'$\delta R_{rs}$ [sr$^{-1}$]')
    axes[1].axhline(0, color='k', lw=0.5)
    for ax in axes:
        ax.set_xlabel('Wavelength [nm]')
        ax.set_xlim(400, 700)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_sensitivity(wave, dRrs, deltas, path):
    """Per-wavelength RMS of each δRrs perturbation vs the RMS of δRrs itself."""
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.semilogy(wave, np.sqrt(np.nanmean(dRrs**2, axis=0)), 'k-', lw=2,
                label=r'RMS $\delta R_{rs}$ (signal)')
    for label, d in deltas.items():
        ax.semilogy(wave, np.sqrt(np.nanmean(d**2, axis=0)) + 1e-12, lw=1, label=label)
    ax.axhline(1e-6 / np.sqrt(12), color='0.5', ls=':', lw=1,
               label='deposit rounding (1e-6/√12)')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel(r'RMS over samples [sr$^{-1}$]')
    ax.set_xlim(400, 700)
    ax.legend(fontsize=7.5, frameon=False, loc='upper left', bbox_to_anchor=(1.01, 1.0))
    ax.set_title(r'Sensitivity of $\delta R_{rs}$ to implementation choices')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    T = data.meta['temp'].to_numpy()
    S = data.meta['sal'].to_numpy()
    hplc = data.pigments['Tchla'].to_numpy()
    tables = gsm.load_ref_tables(data.wave)
    summary = {'script': 'scripts/sdp/reproduce_gsm.py', 'epft_up_version': __version__,
               'input': data.provenance, 'tables': tables.source,
               'run_utc': datetime.datetime.now(datetime.timezone.utc)
               .isoformat(timespec='seconds')}

    # 1. baseline fit
    fit = gsm.fit_gsm(data.wave, data.Rrs, T, S, tables=tables)
    summary['fit'] = {'settings': fit.settings,
                      'n_converged': int(fit.converged.sum()),
                      'n_iter_quantiles': np.percentile(fit.n_iter, [50, 90, 100]).tolist(),
                      'n_negative': dict(zip(gsm.PARAMS, (fit.params < 0).sum(0).tolist())),
                      'negative_bbp_campaigns': data.meta['campaign'][fit.params[:, 2] < 0]
                      .value_counts().to_dict()}
    print(summary['fit'])

    # 2. Fig. 4
    chl_oc4 = oc4v6(data.wave, data.Rrs)
    logres = np.log10(fit.params[:, 0] / hplc)
    worst = int(np.nanargmax(np.abs(logres)))
    keep = np.ones(data.n, bool)
    keep[worst] = False
    summary['fig4'] = {
        'paper': {'OC4': {'R2': 0.75, 'slope': 0.87}, 'GSM': {'R2': 0.86, 'slope': 0.96}},
        'OC4v6': log_stats(hplc, chl_oc4),
        'GSM_all': log_stats(hplc, fit.params[:, 0]),
        'GSM_without_worst': log_stats(hplc[keep], fit.params[keep, 0]),
        'worst': {'index': worst, 'campaign': data.meta['campaign'][worst],
                  'time': str(data.meta['time'][worst]),
                  'hplc': float(hplc[worst]),
                  'params': dict(zip(gsm.PARAMS, fit.params[worst].tolist())),
                  'nm_trials': count_nonpositive_trials(fit, worst)},
        'log10_residual_by_campaign': pd.Series(logres).groupby(
            data.meta['campaign'].to_numpy()).agg(['mean', 'std', 'count'])
        .round(3).to_dict(orient='index')}
    for k in ('OC4v6', 'GSM_all', 'GSM_without_worst'):
        s = summary['fig4'][k]
        print(f'Fig4 {k:18s} R2={s["R2"]:.3f} slope={s["slope"]:.3f} n={s["n"]}')
    print('worst fit:', summary['fig4']['worst'])
    fig4(data, chl_oc4, fit, worst, OUT / 'kramer_fig4_chl.png')
    fig2bc(data, fit, OUT / 'kramer_fig2bc_model_residual.png')

    # 4. cross-checks against the Python port
    dRrs_port, bsw_port = run_port(data)
    bsw_ours = gsm.bsw_zhang2009(data.wave, T, S)
    fit_Rrs_eta = gsm.fit_gsm(data.wave, data.Rrs, T, S, tables=tables, eta_from='Rrs')
    summary['port'] = {
        'bsw_max_rel_diff': float(np.max(np.abs(bsw_ours / bsw_port - 1))),
        'dRrs_vs_ours_rel_rms': rel_rms(dRrs_port - fit.dRrs, fit.dRrs),
        'dRrs_vs_ours_eta_Rrs_rel_rms': rel_rms(dRrs_port - fit_Rrs_eta.dRrs,
                                                fit_Rrs_eta.dRrs),
        'dRrs_vs_ours_eta_Rrs_max_abs': float(np.nanmax(np.abs(dRrs_port
                                                               - fit_Rrs_eta.dRrs))),
        'notes': 'port: eta from above-water Rrs; fmin xtol=ftol=1e-6; no Tchla<=0 guard'}
    print('port cross-check:', summary['port'])

    # 5. sensitivities of δRrs
    deltas = {}
    sens = {}
    woa_csv = sdpdata.kramer2022_dir() / 'woa23_surface_ts.csv'
    if woa_csv.is_file():
        woa = pd.read_csv(woa_csv)
        fit_woa = gsm.fit_gsm(data.wave, data.Rrs, woa['T_woa'].to_numpy(),
                              woa['S_woa'].to_numpy(), tables=tables)
        d = fit_woa.dRrs - fit.dRrs
        deltas['WOA23 − in-situ T/S'] = d
        sens['woa_vs_insitu'] = {
            'dT_insitu_minus_woa': pd.Series(T - woa['T_woa']).describe().round(3).to_dict(),
            'dS_insitu_minus_woa': pd.Series(S - woa['S_woa']).describe().round(3).to_dict(),
            'dRrs_rel_rms': rel_rms(d, fit.dRrs),
            'dRrs_max_abs': float(np.nanmax(np.abs(d))),
            'Tchla_max_rel_change': float(np.nanmax(np.abs(fit_woa.params[:, 0]
                                                           / fit.params[:, 0] - 1))),
            'woa_file_sha256': sdpdata.sha256_of(woa_csv)}
    else:
        sens['woa_vs_insitu'] = 'skipped: run scripts/sdp/fetch_woa_ts.py'
    d = fit_Rrs_eta.dRrs - fit.dRrs
    deltas['η from Rrs (port) − from rrs'] = d
    sens['eta_Rrs_vs_rrs'] = {'dRrs_rel_rms': rel_rms(d, fit.dRrs),
                              'dRrs_max_abs': float(np.nanmax(np.abs(d)))}
    fit_lsq = gsm.fit_gsm(data.wave, data.Rrs, T, S, tables=tables, method='lsq')
    d = fit_lsq.dRrs - fit.dRrs
    deltas['bounded lsq − Nelder–Mead (5 ANT spectra)'] = d
    changed = np.nanmax(np.abs(d), axis=1) > 1e-7
    sens['optimizer_lsq_vs_nm'] = {
        'dRrs_rel_rms': rel_rms(d, fit.dRrs),
        'n_samples_changed_gt_1e-7': int(changed.sum()),
        'changed_campaigns': data.meta['campaign'][changed].value_counts().to_dict(),
        'dRrs_rel_rms_unchanged_samples': rel_rms(d[~changed], fit.dRrs[~changed]),
        'Tchla_lsq_fig4': log_stats(hplc, fit_lsq.params[:, 0])}
    Rrs2 = spectral.moving_mean(data.Rrs, 5)
    inner = np.isfinite(Rrs2[0])
    Rrs2[:, ~inner] = data.Rrs[:, ~inner]     # leave the 2 edge nm untouched
    fit_dbl = gsm.fit_gsm(data.wave, Rrs2, T, S, tables=tables)
    d = fit_dbl.dRrs - fit.dRrs
    deltas['second 5 nm moving mean'] = d
    sens['double_smoothing'] = {'dRrs_rel_rms': rel_rms(d, fit.dRrs),
                                'dRrs_max_abs': float(np.nanmax(np.abs(d)))}
    summary['sensitivity'] = sens
    for k, v in sens.items():
        print('sensitivity', k, v)
    fig_sensitivity(data.wave, fit.dRrs, deltas, OUT / 'gsm_sensitivity.png')

    # 6. save the product for Execution #3
    prod_dir = sdpdata.kramer2022_dir() / 'products'
    prod_dir.mkdir(exist_ok=True)
    prod = prod_dir / 'gsm_dRrs_insitu.npz'
    np.savez_compressed(prod, wave=fit.wave, dRrs=fit.dRrs, Rrs_mod=fit.Rrs_mod,
                        params=fit.params, converged=fit.converged,
                        provenance=json.dumps({k: summary[k] for k in
                                               ('script', 'epft_up_version', 'input',
                                                'tables', 'run_utc')} |
                                              {'settings': fit.settings}, default=str))
    summary['product'] = {'file': str(prod), 'sha256': sdpdata.sha256_of(prod)}
    with open(OUT / 'gsm_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f'wrote {OUT / "gsm_summary.json"} and {prod}')


if __name__ == '__main__':
    main()
