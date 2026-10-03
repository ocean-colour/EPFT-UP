"""Execution #3: reproduce the Kramer 2022 PCR pigment model (§2.5, Table 2, Figs. 3, 6).

Steps
-----
1. Predictors: ``diff(δRrs, 2)`` on the 1 nm grid (299 bands, 401-699 nm)
   from the Execution #2 product. Run on all 145 spectra and on the 144
   without the degenerate SABOR GSM fit (report §4.1).
2. For each of the 13 pigments: ``train_pcr_kramer`` (100 permutations,
   max 30 PCs, 5-fold inner CV on MAE, outputs ≥ 0). Table 2 statistics:
   mean ± SD of R² and of MAE normalized by the mean modelled value.
3. Full reconstruction: median of the 100 models applied to every spectrum,
   clipped at 0, values below the detection limit set to 0. The detection
   limit is a proxy (see ``lod_proxy``). Fig. 6 (Tchla, Fuco, Perid,
   HexFuco, MVchlb, Zea) and Fig. 3 (clustering of the 12 accessory-pigment
   ratios to Tchla, measured vs modelled; Ward linkage on 1 − R).
4. Coefficient check: the port's ``original_a_coefs.xlsx`` /
   ``original_c_coefs.xlsx`` (100 trained models per pigment) applied to our
   δRrs''; correlation of median coefficient spectra.
5. Variants of the paper's supplement for all 13 pigments: measured
   Rrs' + Rrs'' (forward diff + diff 2); δRrs'' at 5 and 10 nm sampling.

Outputs: ``reports/figures/sdp/kramer_table2.csv``,
``reports/figures/sdp/kramer_fig6_pigments.png``,
``reports/figures/sdp/kramer_fig3_dendrograms.png``,
``reports/figures/sdp/pcr_coefficients.png``,
``reports/figures/sdp/pcr_summary.json``, and the coefficient ensembles plus
reconstructions in ``$OS_COLOR/PANGAEA/Kramer2022/products/pcr_rrsD2_1nm.npz``.

Run with::

    conda run -n ocean14 python scripts/sdp/reproduce_pcr.py
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
from scipy.cluster import hierarchy  # noqa: E402
from scipy.spatial.distance import pdist  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import data as sdpdata, models, spectral  # noqa: E402
from epft_up.sdp import validate as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
PIGS = list(sdpdata.PIGMENTS_13)
ACCESSORY = [p for p in PIGS if p != 'Tchla']
SEED = 1

#: Kramer 2022 Table 2: mean R², SD R², mean normalized MAD, SD normalized MAD.
TABLE2 = {
    'Allo': (0.40, 0.19, 1.221, 0.400), 'ButFuco': (0.62, 0.16, 0.588, 0.185),
    'Chlc3': (0.68, 0.13, 0.639, 0.212), 'Chlc12': (0.70, 0.13, 0.703, 0.235),
    'DVchla': (0.55, 0.12, 0.594, 0.103), 'Fuco': (0.65, 0.15, 0.844, 0.274),
    'HexFuco': (0.54, 0.16, 0.692, 0.201), 'MVchlb': (0.42, 0.19, 0.975, 0.295),
    'Neo': (0.42, 0.21, 1.127, 0.354), 'Perid': (0.49, 0.13, 0.783, 0.166),
    'Tchla': (0.72, 0.15, 0.498, 0.127), 'Viola': (0.38, 0.18, 1.101, 0.370),
    'Zea': (0.37, 0.10, 0.472, 0.071)}
#: Kramer 2022 Fig. 3 linkage cut-offs (measured, modelled).
CUT_MEASURED, CUT_MODELLED = 0.65, 0.80
FIG6 = ['Tchla', 'Fuco', 'Perid', 'HexFuco', 'MVchlb', 'Zea']
SOURCE_COLORS = {
    'ANT': 'red', 'NAAMES': 'orange', 'RemSensPOC': 'gold', 'SABOR': 'green',
    'Tara Oceans': 'blue', 'Tara Med': 'blue', 'BIOSOPE': 'purple',
    'EXPORTS': 'black'}


#: Detection-limit proxy (single definition in epft_up.sdp.validate).
lod_proxy = V.lod_proxy


def train_all(X, pig, keep, lods):
    """Train + reconstruct every pigment; returns results and reconstructions."""
    res, recon = {}, {}
    for p in PIGS:
        y = pig[p].to_numpy()[keep]
        r = models.train_pcr_kramer(X[keep], y, seed=SEED)
        res[p] = r
        recon[p], _ = models.predict_ensemble(X[keep], r.coefs, r.intercepts,
                                              lod=lods[p])
    return res, recon


def table2_frame(res):
    rows = []
    for p in PIGS:
        s = res[p].summary()
        t = TABLE2[p]
        rows.append({'pigment': p, 'R2': s['R2'][0], 'R2_sd': s['R2'][1],
                     'R2_paper': t[0], 'R2_sd_paper': t[1],
                     'MADn': s['MAE_norm_pred'][0], 'MADn_sd': s['MAE_norm_pred'][1],
                     'MADn_paper': t[2], 'MADn_sd_paper': t[3],
                     'MADn_obs': s['MAE_norm_obs'][0],
                     'n_pcs_median': float(np.median(res[p].n_pcs)),
                     'dR2_in_sd': (s['R2'][0] - t[0]) / t[1]})
    return pd.DataFrame(rows).set_index('pigment')


def linfit(x, y):
    """OLS y = a + b x; returns slope, intercept, R²."""
    b, a = np.polyfit(x, y, 1)
    return float(b), float(a), float(np.corrcoef(x, y)[0, 1]**2)


def ratio_clusters(conc, cut):
    """Ward clustering of pigment:Tchla ratios on correlation distance (1 − R)."""
    ok = conc['Tchla'] > 0
    ratios = conc.loc[ok, ACCESSORY].div(conc.loc[ok, 'Tchla'], axis=0)
    D = pdist(ratios.to_numpy().T, metric='correlation')
    Z = hierarchy.linkage(D, method='ward')
    labels = hierarchy.fcluster(Z, t=cut, criterion='distance')
    groups = {}
    for name, lab in zip(ACCESSORY, labels):
        groups.setdefault(int(lab), []).append(name)
    return Z, sorted(groups.values(), key=lambda g: ACCESSORY.index(g[0])), int(ok.sum())


#: Kramer 2022 Fig. 6 panel fits (log10 modelled on log10 HPLC): slope, intercept, R².
FIG6_PAPER = {'Tchla': (0.94, -0.004, 0.73), 'Fuco': (0.80, -0.08, 0.60),
              'Perid': (0.78, -0.46, 0.51), 'HexFuco': (0.76, -0.20, 0.64),
              'MVchlb': (0.74, -0.26, 0.51), 'Zea': (0.53, -0.59, 0.55)}


def fig6(meta, pig, recon, keep, path):
    """Kramer Fig. 6: log10 modelled vs log10 HPLC (zeros on either axis dropped)."""
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.4))
    out = {}
    for ax, p in zip(axes.ravel(), FIG6):
        x = pig[p].to_numpy()[keep]
        y = recon[p]
        camp = meta['campaign'].to_numpy()[keep]
        pos = (x > 0) & (y > 0)
        lx, ly = np.log10(x[pos]), np.log10(y[pos])
        for c, col in SOURCE_COLORS.items():
            s = camp[pos] == c
            ax.scatter(lx[s], ly[s], s=12, color=col, edgecolor='none')
        b, a, r2 = linfit(lx, ly)
        bl, al, r2l = linfit(x, y)
        pb, pa, pr2 = FIG6_PAPER[p]
        out[p] = {'log10': {'slope': b, 'intercept': a, 'R2': r2, 'n': int(pos.sum())},
                  'linear_all': {'slope': bl, 'intercept': al, 'R2': r2l, 'n': int(len(x))},
                  'paper_log10': {'slope': pb, 'intercept': pa, 'R2': pr2},
                  'n_zero_measured': int((x == 0).sum()),
                  'n_zero_modelled': int((y == 0).sum())}
        ax.plot([-3, 1], [-3, 1], 'k:', lw=0.8)
        ax.plot([-3, 1], [a - 3 * b, a + b], 'r-', lw=1)
        ax.set_xlim(-3, 1)
        ax.set_ylim(-3, 1)
        ax.text(0.04, 0.96, f'y = {b:.2f}x{a:+.2f}, R² = {r2:.2f} (N={pos.sum()})\n'
                f'paper: y = {pb:.2f}x{pa:+.2f}, R² = {pr2:.2f}',
                transform=ax.transAxes, va='top', fontsize=7.5)
        ax.set_title(p, fontsize=9)
        ax.set_xlabel(r'log$_{10}$ HPLC [mg m$^{-3}$]', fontsize=8)
        ax.set_ylabel(r'log$_{10}$ modelled [mg m$^{-3}$]', fontsize=8)
    fig.suptitle('Kramer 2022 Fig. 6 reproduction (median of 100 PCR models, δRrs″ 1 nm, N=144)',
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return out


def fig3(Zm, Zp, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, Z, cut, title in ((axes[0], Zm, CUT_MEASURED, '(A) measured HPLC ratios'),
                              (axes[1], Zp, CUT_MODELLED, '(B) PCR-modelled ratios')):
        hierarchy.dendrogram(Z, labels=ACCESSORY, ax=ax, color_threshold=cut,
                             leaf_rotation=90)
        ax.axhline(cut, color='r', ls='--', lw=1)
        ax.set_title(f'{title} (cut {cut})', fontsize=10)
        ax.set_ylabel('linkage distance (Ward, 1 − R)')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def load_port_coefs():
    base = (sdpdata.kramer2022_dir() / 'ref' / 'rrs-SDP-pigments' / 'src' / 'sdp'
            / 'resources' / 'sdp_coefs')
    A, C = {}, {}
    for p in PIGS:
        a = pd.read_excel(base / 'original_a_coefs.xlsx', sheet_name=p)
        if not np.array_equal(a['wave'].to_numpy(), np.arange(401, 700)):
            raise ValueError('unexpected wavelength grid in original_a_coefs.xlsx')
        A[p] = a.drop(columns='wave').to_numpy().T          # (100, 299)
        C[p] = pd.read_excel(base / 'original_c_coefs.xlsx', sheet_name=p).to_numpy().ravel()
    return A, C, {'a': str(base / 'original_a_coefs.xlsx'),
                  'c': str(base / 'original_c_coefs.xlsx'),
                  'sha256_a': sdpdata.sha256_of(base / 'original_a_coefs.xlsx'),
                  'sha256_c': sdpdata.sha256_of(base / 'original_c_coefs.xlsx')}


def fig_coefs(wave_d, res, A_port, path):
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    for ax, p in zip(axes, ['Tchla', 'Fuco', 'HexFuco']):
        ours = res[p].coefs
        q = np.percentile(ours, [25, 50, 75], axis=0)
        ax.fill_between(wave_d, q[0], q[2], color='C0', alpha=0.3, lw=0)
        ax.plot(wave_d, q[1], 'C0-', lw=1, label='this work (median, IQR)')
        ax.plot(wave_d, np.median(A_port[p], axis=0), 'C3-', lw=0.8,
                label='port original_a_coefs (median)')
        ax.axhline(0, color='k', lw=0.4)
        r = np.corrcoef(q[1], np.median(A_port[p], axis=0))[0, 1]
        ax.set_title(f'{p}: A(λ) on δRrs″ (r = {r:.2f})', fontsize=9)
        ax.set_ylabel(r'A(λ) [mg m$^{-3}$ sr]')
    axes[0].legend(fontsize=8, frameon=False)
    axes[-1].set_xlabel('Wavelength [nm]')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    prod_dir = sdpdata.kramer2022_dir() / 'products'
    gsm_prod = prod_dir / 'gsm_dRrs_insitu.npz'
    z = np.load(gsm_prod)
    wave_d, X = spectral.difference(z['wave'], z['dRrs'], order=2)
    pig = data.pigments13
    lods = lod_proxy(pig)
    hplc_logres = np.log10(z['params'][:, 0] / pig['Tchla'].to_numpy())
    worst = int(np.argmax(np.abs(hplc_logres)))
    keep_all = np.ones(data.n, dtype=bool)
    keep_144 = keep_all.copy()
    keep_144[worst] = False

    summary = {'script': 'scripts/sdp/reproduce_pcr.py', 'epft_up_version': __version__,
               'input': data.provenance,
               'gsm_product': {'file': str(gsm_prod), 'sha256': sdpdata.sha256_of(gsm_prod)},
               'predictors': {'type': 'diff(dRrs, 2), 1 nm', 'n_features': int(X.shape[1]),
                              'wave_range': [float(wave_d[0]), float(wave_d[-1])]},
               'lod_proxy': lods, 'seed': SEED, 'excluded_for_144': worst,
               'run_utc': datetime.datetime.now(datetime.timezone.utc)
               .isoformat(timespec='seconds')}

    # 2-3. main runs: 145 and 144
    runs = {}
    for label, keep in (('n145', keep_all), ('n144', keep_144)):
        res, recon = train_all(X, pig, keep, lods)
        t2 = table2_frame(res)
        runs[label] = (res, recon, t2, keep)
        print(f'--- Table 2, {label} ---')
        print(t2[['R2', 'R2_sd', 'R2_paper', 'MADn', 'MADn_sd', 'MADn_paper',
                  'n_pcs_median']].round(3).to_string())
    summary['settings'] = runs['n145'][0]['Tchla'].settings
    summary['table2'] = {k: json.loads(v[2].to_json(orient='index')) for k, v in runs.items()}
    runs['n144'][2].round(4).to_csv(OUT / 'kramer_table2.csv')

    res, recon, t2, keep = runs['n144']
    summary['fig6'] = fig6(data.meta, pig, recon, keep, OUT / 'kramer_fig6_pigments.png')
    for p, s6 in summary['fig6'].items():
        print(f"Fig 6 {p:8s} log10 {s6['log10']}  paper {s6['paper_log10']}  "
              f"zeros meas/mod {s6['n_zero_measured']}/{s6['n_zero_modelled']}")

    meas = pig.loc[keep].reset_index(drop=True)
    modl = pd.DataFrame(recon)[PIGS]
    Zm, gm, nm = ratio_clusters(meas, CUT_MEASURED)
    Zp, gp, npred = ratio_clusters(modl, CUT_MODELLED)
    fig3(Zm, Zp, OUT / 'kramer_fig3_dendrograms.png')
    summary['fig3'] = {'measured_groups': gm, 'modelled_groups': gp,
                       'n_measured': nm, 'n_modelled_with_Tchla_gt0': npred}
    print('Fig 3 measured groups:', gm)
    print('Fig 3 modelled groups:', gp)

    # 4. the port's original coefficients on our predictors
    A_port, C_port, port_src = load_port_coefs()
    port = {}
    for p in PIGS:
        y = pig[p].to_numpy()
        med_port, _ = models.predict_ensemble(X, A_port[p], C_port[p], lod=lods[p])
        med_ours, _ = models.predict_ensemble(X, runs['n145'][0][p].coefs,
                                              runs['n145'][0][p].intercepts, lod=lods[p])
        port[p] = {
            'R2_full_recon_port_coefs': float(np.corrcoef(med_port, y)[0, 1]**2),
            'R2_full_recon_ours_n145': float(np.corrcoef(med_ours, y)[0, 1]**2),
            'corr_median_coef_spectra': float(np.corrcoef(
                np.median(A_port[p], axis=0),
                np.median(runs['n145'][0][p].coefs, axis=0))[0, 1]),
            'median_intercept_port': float(np.median(C_port[p])),
            'median_intercept_ours': float(np.median(runs['n145'][0][p].intercepts))}
    summary['port_coefficients'] = {'source': port_src, 'per_pigment': port}
    print(pd.DataFrame(port).T.round(3).to_string())
    fig_coefs(wave_d, runs['n145'][0], A_port, OUT / 'pcr_coefficients.png')

    # 5. supplementary variants
    variants = {}
    _, Xd12 = spectral.derivative_features(data.wave, data.Rrs, orders=(1, 2))
    _, X5 = spectral.derivative_features(z['wave'], z['dRrs'], orders=(2,), step=5)
    _, X10 = spectral.derivative_features(z['wave'], z['dRrs'], orders=(2,), step=10)
    for vname, Xv in (('Rrs_d1_d2_1nm', Xd12), ('dRrs_d2_5nm', X5), ('dRrs_d2_10nm', X10)):
        vres = {p: models.train_pcr_kramer(Xv[keep_144], pig[p].to_numpy()[keep_144],
                                           seed=SEED).summary() for p in PIGS}
        variants[vname] = {'n_features': int(Xv.shape[1]),
                           'R2': {p: vres[p]['R2'] for p in PIGS},
                           'MADn': {p: vres[p]['MAE_norm_pred'] for p in PIGS}}
    base = {p: runs['n144'][0][p].summary()['R2'] for p in PIGS}
    vt = pd.DataFrame({'dRrs_d2_1nm': {p: base[p][0] for p in PIGS},
                       **{v: {p: variants[v]['R2'][p][0] for p in PIGS} for v in variants}})
    print('--- variant mean R² (n=144) ---')
    print(vt.round(3).to_string())
    summary['variants'] = variants

    # products
    prod = prod_dir / 'pcr_rrsD2_1nm.npz'
    np.savez_compressed(
        prod, wave_d=wave_d, pigments=np.array(PIGS),
        **{f'coefs_{lab}_{p}': runs[lab][0][p].coefs for lab in runs for p in PIGS},
        **{f'intercepts_{lab}_{p}': runs[lab][0][p].intercepts for lab in runs for p in PIGS},
        **{f'recon_{lab}_{p}': runs[lab][1][p] for lab in runs for p in PIGS},
        keep_144=keep_144,
        provenance=json.dumps({k: summary[k] for k in ('script', 'epft_up_version',
                                                       'gsm_product', 'predictors',
                                                       'settings', 'seed', 'run_utc')},
                              default=str))
    summary['product'] = {'file': str(prod), 'sha256': sdpdata.sha256_of(prod)}
    with open(OUT / 'pcr_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f'wrote {OUT / "pcr_summary.json"} and {prod}')


if __name__ == '__main__':
    main()
