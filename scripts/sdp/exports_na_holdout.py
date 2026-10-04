"""Execution #9 (partial): EXPORTS North Atlantic 2021 hold-out.

Status: **PROVISIONAL, Tchla only.** The full prompt needs the EXPORTS-NA
HPLC pigments from SeaBASS. Programmatic SeaBASS access is refused
(file_search.cgi: HTTP 403/444, also with Earthdata credentials; the archive
pages list files only via JavaScript), so the pigment matchups could not be
built. See the Execution #9 log entry.

What *is* available, unambiguously: ``Kramer_rrs_testdata.mat`` in the pinned
``sashajane19/Rrs_pigments`` repo. It holds the 17 EXPORTS-NA (May 2021)
hyperspectral Rrs (400-700 nm, 1 nm) with T, S, lat/lon and HPLC
chlorophyll-a; its README says "EXPORTS North Atlantic hyperspectral remote
sensing reflectances and HPLC chlorophyll data". These are the 17 samples
Kramer et al. (2024) added. Caveat: these spectra are **not processed like
the deposit**. They are not rounded to 1e-6 and are far smoother than the
5 nm mean leaves (2nd-difference power at 0.3-0.5 cycles/nm well below the
deposit's floor), and one spectrum is clipped to 0 at 697-700 nm. The
preprocessing difference is itself a test of transfer.

Without retraining (all models fitted on the 145 deposit samples):

* PCR on δRrs'' (Execution #3 ensemble, N = 145) and Kramer's original
  coefficients (port);
* Tchla nulls: log-log calibration of GSM and OC4 Tchla on the 145;
* the shared W with Bayesian intervals (Execution #8, cross-fitted
  calibration), absolute targets;
* the raw GSM Tchla.

If an EXPORTS-NA HPLC table is supplied later (``--hplc path.csv`` with
columns lat, lon and the 13 pigment names), the pigment and ratio targets
are added by matching on position.

Outputs: ``reports/figures/sdp/exports_na_tchla.png`` and
``reports/figures/sdp/exports_na_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/exports_na_holdout.py
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
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import bayes, data as sdpdata, gsm, models, noise, spectral  # noqa: E402
from epft_up.sdp import validate as V  # noqa: E402
from ingest_kramer2022 import d2_power_vs_quantization  # noqa: E402
from reproduce_pcr import load_port_coefs  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
SEED = 1
LAM_NOISE = (0.0, 0.1, 1.0, 10.0, 100.0)
LAM_SMOOTH = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0)
RANKS = tuple(range(1, 11))


def load_exports_na():
    """The 17 EXPORTS-NA spectra from the pinned Rrs_pigments test file."""
    f = sdpdata.kramer2022_dir() / 'ref' / 'Rrs_pigments' / 'Kramer_rrs_testdata.mat'
    d = loadmat(f)
    return {'file': str(f), 'sha256': sdpdata.sha256_of(f),
            'Rrs': d['Rrs'].astype(float), 'T': d['T'].ravel(), 'S': d['S'].ravel(),
            'chl': d['chl'].ravel(), 'lat': d['latlon'][:, 0], 'lon': d['latlon'][:, 1]}


def scores(obs, pred):
    lo, lp = np.log10(obs), np.log10(np.clip(pred, 1e-3, None))
    return {'n': int(len(obs)), 'bias_log10': float(np.mean(lp - lo)),
            'rms_log10': float(np.sqrt(np.mean((lp - lo)**2))),
            'R2_log': float(np.corrcoef(lo, lp)[0, 1]**2),
            'median_ratio': float(np.median(10**(lp - lo))),
            'n_nonpositive_pred': int(np.sum(pred <= 0))}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    z = np.load(sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz')
    pcr = np.load(sdpdata.kramer2022_dir() / 'products' / 'pcr_rrsD2_1nm.npz')
    ex = load_exports_na()
    w = data.wave

    # processing check vs the deposit
    proc = {'on_1e-6_grid': bool(np.allclose(np.round(ex['Rrs'] * 1e6), ex['Rrs'] * 1e6,
                                             atol=1e-6)),
            'd2_power_over_deposit_rounding_floor': d2_power_vs_quantization(ex['Rrs']),
            'n_nonpositive_Rrs': int(np.sum(ex['Rrs'] <= 0)),
            'chl_range': [float(ex['chl'].min()), float(ex['chl'].max())]}

    # GSM on EXPORTS-NA (same code and tables as Execution #2)
    tables = gsm.load_ref_tables(w)
    fit = gsm.fit_gsm(w, ex['Rrs'], ex['T'], ex['S'], tables=tables)
    dR = fit.dRrs
    chl_gsm = fit.params[:, 0]
    oc4_ex = V.oc4v6(w, np.clip(ex['Rrs'], 1e-6, None))
    _, Xd2 = spectral.difference(w, dR, 2)

    preds, sd_w = {}, None
    # PCR (ours, N=145) and Kramer's original coefficients
    lod = V.lod_proxy(data.pigments13)['Tchla']
    preds['PCR (this work, N=145)'], _ = models.predict_ensemble(
        Xd2, pcr['coefs_n145_Tchla'], pcr['intercepts_n145_Tchla'], lod=lod)
    A, C, _ = load_port_coefs()
    preds['PCR (Kramer original coefs)'], _ = models.predict_ensemble(Xd2, A['Tchla'],
                                                                      C['Tchla'], lod=lod)
    # nulls calibrated on the 145
    hplc = data.pigments['Tchla'].to_numpy()
    for name, x145, xex in (('null GSM', z['params'][:, 0], chl_gsm),
                            ('null OC4', V.oc4v6(w, data.Rrs), oc4_ex)):
        b, a = np.polyfit(np.log10(x145), np.log10(hplc), 1)
        preds[name] = 10**(a + b * np.log10(np.clip(xex, 1e-3, None)))
    preds['GSM Tchla (raw)'] = chl_gsm
    # shared W (Bayesian, cross-fitted calibration), absolute targets, trained on 145
    bench = V.Benchmark(data.pigments13, z['params'][:, 0], V.oc4v6(w, data.Rrs),
                        data.meta['campaign'].to_numpy(), n_perm=2, seed=SEED)
    keys = [t for t, v in bench.targets.items() if v['kind'] == 'abs']
    Y = np.column_stack([bench.targets[t]['y'] for t in keys])
    p = z['params']
    aux145 = np.column_stack([np.log10(p[:, 0]), np.log10(p[:, 1]), p[:, 2]])
    Sn = noise.covariance(w, np.median(data.Rrs, axis=0), 'insitu')
    m = bayes.BayesianSharedW(LAM_NOISE, LAM_SMOOTH, RANKS, Sigma_n=Sn)
    m.fit(z['dRrs'], aux145, Y, data.meta['campaign'].to_numpy())
    pe = fit.params
    aux_ex = np.column_stack([np.log10(np.clip(pe[:, 0], 1e-3, None)),
                              np.log10(np.clip(pe[:, 1], 1e-4, None)), pe[:, 2]])
    mu, sd = m.predict(dR, aux_ex)
    j = keys.index('abs:Tchla')
    preds['shared W (Bayesian)'] = np.maximum(mu[:, j], 0.0)
    sd_w = sd[:, j]

    res = {k: scores(ex['chl'], v) for k, v in preds.items()}
    g = bayes.gaussian_scores(ex['chl'], mu[:, j], sd_w)
    res['shared W (Bayesian)'].update({'cov68': g['cov68'], 'cov95': g['cov95'],
                                       'z_sd': g['z_sd'], 'z_mean': float(np.mean(
                                           (ex['chl'] - mu[:, j]) / sd_w)),
                                       'sd_median': float(np.median(sd_w))})
    pd.set_option('display.width', 200)
    print(pd.DataFrame(res).T.round(3).to_string())
    print('processing check:', proc)

    fig, ax = plt.subplots(figsize=(6, 5))
    for k, v in preds.items():
        if k == 'shared W (Bayesian)':
            ax.errorbar(ex['chl'], v, yerr=sd_w, fmt='o', ms=5, lw=0.8, label=k)
        else:
            ax.plot(ex['chl'], np.clip(v, 0.05, None), 'o', ms=4, alpha=0.8, label=k)
    lim = [0.05, 3]
    ax.plot(lim, lim, 'k-', lw=0.6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.4, 1.6)
    ax.set_ylim(0.05, 3)
    ax.set_xlabel(r'HPLC Tchla [mg m$^{-3}$] (EXPORTS-NA 2021)')
    ax.set_ylabel(r'predicted Tchla [mg m$^{-3}$] (≤0.05 shown at 0.05)')
    ax.legend(fontsize=7, frameon=False)
    ax.set_title('EXPORTS-NA hold-out, Tchla only (models trained on the 145)', fontsize=9)
    fig.savefig(OUT / 'exports_na_tchla.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    out = {'script': 'scripts/sdp/exports_na_holdout.py', 'epft_up_version': __version__,
           'status': 'provisional: Tchla only; EXPORTS-NA pigments not obtainable',
           'exports_na_source': {k: ex[k] for k in ('file', 'sha256')},
           'processing_check': proc, 'gsm_converged': int(fit.converged.sum()),
           'scores': res, 'shared_W_hp': m.hp,
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(OUT / 'exports_na_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f'wrote {OUT / "exports_na_summary.json"}')


if __name__ == '__main__':
    main()
