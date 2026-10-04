"""Execution #9: EXPORTS North Atlantic 2021, an independent hold-out campaign.

Test data:

* the 17 EXPORTS-NA spectra of ``Kramer_rrs_testdata.mat`` (pinned
  ``sashajane19/Rrs_pigments``), as used by Kramer et al. (2024);
* their HPLC pigments from SeaBASS (UCSB/CRSEO rosette HPLC, JC214 and DY131),
  matched unambiguously by ``scripts/sdp/exports_na_matchups.py`` (position
  ≤ 0.28 km and exact replicate-mean Tchla).

Caveat: the 17 spectra are not processed like the deposit. They are not
rounded to 1e-6, far smoother at high spectral frequency, and one is clipped
to 0 at 697-700 nm. Their source radiometer is not NASA GSFC's HyperSAS (no
match). Per the user's choice they are used as-is, and the difference is part
of the transfer test.

Without any retraining (everything fitted on the 145 deposit samples), on 13
absolute pigments and 12 log10 pigment:Tchla ratios:

* PCR on δRrs'' (Execution #3 ensemble, N = 145) and Kramer's original
  coefficients (port); for ratios, the ratio of PCR's absolute predictions
  ("derived");
* Tchla nulls: log-log calibrations on GSM and OC4 Tchla (for ratios, the
  log ratio on log Tchla);
* the constant (145-sample mean) log ratio;
* the shared W with Bayesian intervals (Execution #8, cross-fitted
  calibration), absolute and ratio targets.

Outputs: ``reports/figures/sdp/exports_na_scores.csv``,
``exports_na_vs_null.csv`` (paired bootstrap of ΔRMS against the best null),
``exports_na_tchla.png``, ``exports_na_ratios.png`` and
``exports_na_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/exports_na_matchups.py   # first
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
N_BOOT = 2000
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


def r2(a, b):
    return float(np.corrcoef(a, b)[0, 1]**2) if np.std(a) > 0 and np.std(b) > 0 else np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    z = np.load(sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz')
    pcr = np.load(sdpdata.kramer2022_dir() / 'products' / 'pcr_rrsD2_1nm.npz')
    mfile = sdpdata.kramer2022_dir() / 'EXPORTS_NA' / 'exports_na_matchups.csv'
    M = pd.read_csv(mfile).sort_values('test_index').reset_index(drop=True)
    ex = load_exports_na()
    if not np.allclose(M['test_chl'], ex['chl']):
        raise SystemExit('matchup table does not line up with the test spectra')
    w = data.wave
    PIGS = list(sdpdata.PIGMENTS_13)
    ACC = [p for p in PIGS if p != 'Tchla']
    lods = V.lod_proxy(data.pigments13)

    proc = {'on_1e-6_grid': bool(np.allclose(np.round(ex['Rrs'] * 1e6), ex['Rrs'] * 1e6,
                                             atol=1e-6)),
            'd2_power_over_deposit_rounding_floor': d2_power_vs_quantization(ex['Rrs']),
            'n_nonpositive_Rrs': int(np.sum(ex['Rrs'] <= 0))}

    # features on EXPORTS-NA (same code and tables as Execution #2)
    tables = gsm.load_ref_tables(w)
    fit = gsm.fit_gsm(w, ex['Rrs'], ex['T'], ex['S'], tables=tables)
    dR, pe = fit.dRrs, fit.params
    _, Xd2 = spectral.difference(w, dR, 2)
    chl_gsm_ex = np.clip(pe[:, 0], 1e-3, None)
    oc4_ex = V.oc4v6(w, np.clip(ex['Rrs'], 1e-6, None))
    aux_ex = np.column_stack([np.log10(chl_gsm_ex), np.log10(np.clip(pe[:, 1], 1e-4, None)),
                              pe[:, 2]])
    chl_gsm_145 = z['params'][:, 0]
    oc4_145 = V.oc4v6(w, data.Rrs)
    p145 = z['params']
    aux_145 = np.column_stack([np.log10(p145[:, 0]), np.log10(p145[:, 1]), p145[:, 2]])

    truth_abs = {p: M[p].to_numpy(float) for p in PIGS}
    rlog = lambda P, p, T: np.log10(V.replace_zeros(P, lods[p]) / T)  # noqa: E731
    truth_ratio = {p: rlog(truth_abs[p], p, truth_abs['Tchla']) for p in ACC}

    preds_abs, preds_ratio, sd_abs, sd_ratio = {}, {}, {}, {}
    # PCR (ours) and Kramer's original coefficients
    A, C, _ = load_port_coefs()
    for name, getc in (('PCR', lambda p: (pcr[f'coefs_n145_{p}'], pcr[f'intercepts_n145_{p}'])),
                       ('PCR_Kramer_coefs', lambda p: (A[p], C[p]))):
        preds_abs[name] = {p: models.predict_ensemble(Xd2, *getc(p), lod=lods[p])[0]
                           for p in PIGS}
        preds_ratio[name + '_derived'] = {
            p: rlog(preds_abs[name][p], p, V.replace_zeros(preds_abs[name]['Tchla'],
                                                           lods['Tchla'])) for p in ACC}
    # nulls (fitted on the 145)
    pig145 = data.pigments13
    for name, x145, xex in (('null_GSM', chl_gsm_145, chl_gsm_ex), ('null_OC4', oc4_145, oc4_ex)):
        preds_abs[name], preds_ratio[name] = {}, {}
        for p in PIGS:
            y = np.log10(V.replace_zeros(pig145[p].to_numpy(), lods[p]))
            b, a = np.polyfit(np.log10(x145), y, 1)
            preds_abs[name][p] = 10**(a + b * np.log10(xex))
        for p in ACC:
            y = rlog(pig145[p].to_numpy(), p, pig145['Tchla'].to_numpy())
            b, a = np.polyfit(np.log10(x145), y, 1)
            preds_ratio[name][p] = a + b * np.log10(xex)
    preds_abs['GSM_raw'] = {'Tchla': chl_gsm_ex}
    preds_ratio['const'] = {p: np.full(len(M), rlog(pig145[p].to_numpy(), p,
                                                    pig145['Tchla'].to_numpy()).mean())
                            for p in ACC}
    # shared W, Bayesian (cross-fitted), trained on the 145
    bench = V.Benchmark(pig145, chl_gsm_145, oc4_145, data.meta['campaign'].to_numpy(),
                        n_perm=2, seed=SEED)
    Sn = noise.covariance(w, np.median(data.Rrs, axis=0), 'insitu')
    camp = data.meta['campaign'].to_numpy()
    hps = {}
    for kind in ('abs', 'ratio'):
        keys = [t for t, v in bench.targets.items() if v['kind'] == kind]
        Y = np.column_stack([bench.targets[t]['y'] for t in keys])
        m = bayes.BayesianSharedW(LAM_NOISE, LAM_SMOOTH, RANKS, Sigma_n=Sn)
        m.fit(z['dRrs'], aux_145, Y, camp)
        hps[kind] = m.hp
        mu, sd = m.predict(dR, aux_ex)
        target = preds_abs if kind == 'abs' else preds_ratio
        sds = sd_abs if kind == 'abs' else sd_ratio
        target['sharedW'] = {}
        for j, t in enumerate(keys):
            p = t.split(':')[1]
            target['sharedW'][p] = np.maximum(mu[:, j], 0.0) if kind == 'abs' else mu[:, j]
            sds[p] = sd[:, j]

    # scores
    rows = []
    for model, d in preds_abs.items():
        for p, v in d.items():
            o = truth_abs[p]
            lo, lp = np.log10(V.replace_zeros(o, lods[p])), np.log10(V.replace_zeros(v, lods[p]))
            rows.append({'kind': 'abs', 'model': model, 'target': p, 'R2_lin': r2(o, v),
                         'R2_log': r2(lo, lp), 'bias_log': float(np.mean(lp - lo)),
                         'rms_log': float(np.sqrt(np.mean((lp - lo)**2))),
                         'truth_sd_log': float(np.std(lo)),
                         'n_truth_zero': int(np.sum(o == 0))})
    for model, d in preds_ratio.items():
        for p, v in d.items():
            o = truth_ratio[p]
            rows.append({'kind': 'ratio', 'model': model, 'target': p, 'R2_log': r2(o, v),
                         'bias_log': float(np.mean(v - o)),
                         'rms_log': float(np.sqrt(np.mean((v - o)**2))),
                         'truth_sd_log': float(np.std(o)),
                         'n_truth_zero': int(np.sum(truth_abs[p] == 0))})
    S = pd.DataFrame(rows)
    S['crms_log'] = np.sqrt(np.maximum(S['rms_log']**2 - S['bias_log']**2, 0))

    # paired bootstrap (over the 17 samples) of RMS(model) - RMS(best null), log space
    rng = np.random.default_rng(SEED)
    boot = rng.integers(0, len(M), size=(N_BOOT, len(M)))

    def err(kind, model, p):
        if kind == 'abs':
            lo = np.log10(V.replace_zeros(truth_abs[p], lods[p]))
            return np.log10(V.replace_zeros(preds_abs[model][p], lods[p])) - lo
        return preds_ratio[model][p] - truth_ratio[p]
    drows = []
    for kind, cands, nulls in (('abs', ('PCR', 'sharedW'), ('null_GSM', 'null_OC4')),
                               ('ratio', ('PCR_derived', 'sharedW'),
                                ('null_GSM', 'null_OC4', 'const'))):
        for p in (PIGS if kind == 'abs' else ACC):
            ref = min(nulls, key=lambda n: np.mean(err(kind, n, p)**2))
            e0 = err(kind, ref, p)**2
            for c in cands:
                e1 = err(kind, c, p)**2
                d = np.sqrt(e1[boot].mean(1)) - np.sqrt(e0[boot].mean(1))
                drows.append({'kind': kind, 'target': p, 'model': c, 'best_null': ref,
                              'd_rms': float(np.sqrt(e1.mean()) - np.sqrt(e0.mean())),
                              'lo95': float(np.percentile(d, 2.5)),
                              'hi95': float(np.percentile(d, 97.5))})
    D = pd.DataFrame(drows)
    D['verdict'] = np.where(D.hi95 < 0, 'beats null', np.where(D.lo95 > 0, 'worse', 'tie'))
    print(D.round(3).to_string())
    D.round(4).to_csv(OUT / 'exports_na_vs_null.csv', index=False)
    cov = {}
    for kind, sds, tr, pr in (('abs', sd_abs, truth_abs, preds_abs['sharedW']),
                              ('ratio', sd_ratio, truth_ratio, preds_ratio['sharedW'])):
        ys = [tr[p] for p in sds]
        mus = [pr[p] for p in sds]
        g = bayes.gaussian_scores(np.concatenate(ys), np.concatenate(mus),
                                  np.concatenate([sds[p] for p in sds]))
        per = {p: bayes.gaussian_scores(tr[p], pr[p], sds[p]) for p in sds}
        cov[kind] = {'pooled': {k: g[k] for k in ('cov68', 'cov95', 'z_sd', 'crps')},
                     'per_target': {p: {k: v[k] for k in ('cov68', 'cov95', 'z_sd')}
                                    for p, v in per.items()}}
    pd.set_option('display.width', 220)
    piv = S[S.kind == 'ratio'].pivot(index='target', columns='model', values='rms_log')
    print('ratio RMS [dex]:'); print(piv.round(3).to_string())
    piv2 = S[S.kind == 'abs'].pivot(index='target', columns='model', values='rms_log')
    print('abs RMS [dex]:'); print(piv2.round(3).to_string())
    print('coverage:', json.dumps({k: v['pooled'] for k, v in cov.items()}))
    S.round(4).to_csv(OUT / 'exports_na_scores.csv', index=False)

    # figures
    fig, ax = plt.subplots(figsize=(6, 5))
    for k in ('PCR', 'PCR_Kramer_coefs', 'null_GSM', 'null_OC4', 'sharedW'):
        v = preds_abs[k]['Tchla']
        if k == 'sharedW':
            ax.errorbar(truth_abs['Tchla'], v, yerr=sd_abs['Tchla'], fmt='o', ms=5, lw=0.8,
                        label=k)
        else:
            ax.plot(truth_abs['Tchla'], np.clip(v, 0.05, None), 'o', ms=4, alpha=0.8, label=k)
    ax.plot([0.05, 3], [0.05, 3], 'k-', lw=0.6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.4, 1.6)
    ax.set_ylim(0.05, 3)
    ax.set_xlabel(r'HPLC Tchla [mg m$^{-3}$]')
    ax.set_ylabel(r'predicted Tchla [mg m$^{-3}$]')
    ax.legend(fontsize=7, frameon=False)
    ax.set_title('EXPORTS-NA 2021 hold-out: Tchla', fontsize=9)
    fig.savefig(OUT / 'exports_na_tchla.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    show = ['Fuco', 'HexFuco', 'Chlc3', 'Chlc12', 'ButFuco', 'Zea']
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5))
    for ax, p in zip(axes.ravel(), show):
        o = truth_ratio[p]
        for k, mk in (('PCR_derived', 's'), ('null_OC4', '^'), ('const', '_')):
            ax.plot(o, preds_ratio[k][p], mk, ms=5, alpha=0.7, label=k)
        ax.errorbar(o, preds_ratio['sharedW'][p], yerr=sd_ratio[p], fmt='o', ms=5, lw=0.8,
                    label='shared W ±1σ')
        lim = [min(o.min(), -2.5) - 0.1, max(o.max(), 0) + 0.1]
        ax.plot(lim, lim, 'k-', lw=0.6)
        ax.set_title(f'log10 {p}:Tchla', fontsize=9)
        ax.set_xlabel('HPLC', fontsize=8)
        ax.set_ylabel('predicted', fontsize=8)
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle('EXPORTS-NA 2021 hold-out: pigment ratios (no retraining)', fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / 'exports_na_ratios.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    out = {'script': 'scripts/sdp/exports_na_holdout.py', 'epft_up_version': __version__,
           'status': 'complete: 17 samples, 13 pigments, 12 ratios',
           'matchups': {'file': str(mfile), 'sha256': sdpdata.sha256_of(mfile)},
           'exports_na_source': {k: ex[k] for k in ('file', 'sha256')},
           'processing_check': proc, 'gsm_converged': int(fit.converged.sum()),
           'sharedW_hp': hps, 'n_boot': N_BOOT,
           'vs_null': {f'{k}:{m}': g['verdict'].value_counts().to_dict()
                       for (k, m), g in D.groupby(['kind', 'model'])}, 'coverage_sharedW': cov,
           'truth_ranges': {p: [float(truth_abs[p].min()), float(truth_abs[p].max())]
                            for p in PIGS},
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(OUT / 'exports_na_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f'wrote {OUT / "exports_na_summary.json"}')


if __name__ == '__main__':
    main()
