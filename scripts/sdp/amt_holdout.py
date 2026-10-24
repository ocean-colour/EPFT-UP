"""Test on new AMT data #2: frozen-model hold-out on 8 AMT cruises (scheme H2).

Protocol: ``reports/provenance/AMT_benchmark_protocol.md`` (§3 harmonisation
via ``epft_up.sdp.amt.load_amt``; §4 scheme H2; §5 targets; §6 baselines;
§7 metrics). **Nothing is retrained**: every model is the one fitted on the
145 Kramer deposit samples, as in ``scripts/sdp/exports_na_holdout.py``.

Models:
* PCR on δRrs'' (our 100-member ensemble, Execution #3) and Kramer's
  original coefficients; ratios derived from the absolute predictions;
* Tchla nulls (log-log on GSM Tchla, on OC4v6 Tchla) and the constant ratio,
  calibrated on the 145;
* the Bayesian shared W (cross-fitted calibration), absolute and ratio.

Spectra: ``smooth5`` (5 nm moving mean, shrinking window at the grid ends;
primary) and ``asis`` (as delivered; sensitivity).

Scores (pooled and per cruise; MVchlb only on AMT29/30, where it was
measured; DVchla not on AMT27, where it was not separated):
* bias, RMS, centred RMS, R² (log; linear too for absolute targets);
* ΔRMS = RMS(model) − RMS(best null) with two 95% CIs: a paired bootstrap
  over samples, and a **cluster bootstrap over cruises** (primary verdict);
  plus the number of cruises in which the model beats the null;
* shared W: 68/95% coverage, SD(z), mean z, PIT and CRPS.

Campaign-level variance (report §10.3, Rec. 3). τ² (between-campaign) and
σ_w² (within) are split by method of moments from the shared W's
leave-one-campaign-out residuals **on the deposit only**. Two tests:

1. per-sample intervals with τ² added to the predictive variance. The
   cross-fitted variance already contains the between-campaign error on
   average, so this is a conservative bound (it double-counts);
2. **cruise-mean residuals** against ±1.96·√(τ² + σ_w²/n) (campaign-aware)
   versus ±1.96·√(mean σ²/n) (independent errors). This is the
   statement the campaign term is for.

Outputs (``reports/figures/sdp/``): ``amt_holdout_scores.csv``,
``amt_holdout_predictions.csv.gz`` (per sample, model and target),
``amt_holdout_vs_null.csv``, ``amt_holdout_coverage.csv``,
``amt_holdout_campaign.csv``, ``amt_holdout_summary.json`` and figures
``amt_holdout_tchla.png``, ``amt_holdout_ratios.png``,
``amt_holdout_cruise_bias.png``, ``amt_holdout_coverage.png``,
``amt_holdout_smooth_vs_asis.png``.

Run with::

    conda run -n ocean14 python scripts/sdp/amt_harmonise.py   # first
    conda run -n ocean14 python scripts/sdp/amt_holdout.py
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

import matplotlib  # noqa: E402

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import amt as AMT, bayes, data as sdpdata, gsm, models, noise, spectral  # noqa: E402
from epft_up.sdp import validate as V  # noqa: E402
from reproduce_pcr import load_port_coefs  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
CACHE = AMT.amt_dir() / 'products'
SEED = 1
N_BOOT = 2000
LAM_NOISE = (0.0, 0.1, 1.0, 10.0, 100.0)
LAM_SMOOTH = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0)
RANKS = tuple(range(1, 11))
VARIANTS = ('smooth5', 'asis')
KEY_RATIOS = ('Fuco', 'Zea', 'DVchla', 'Chlc12', 'HexFuco')


def r2(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    return float(np.corrcoef(a, b)[0, 1]**2) if len(a) > 2 and a.std() > 0 and b.std() > 0 \
        else np.nan


def gsm_features(d, variant):
    """GSM fit on the AMT spectra (cached by input checksum and variant)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f'amt_gsm_{variant}.npz'
    tag = f"{d.provenance['sha256']}:{variant}:{d.n}"
    if f.exists():
        z = np.load(f, allow_pickle=False)
        if str(z['tag']) == tag:
            return z['dRrs'], z['params'], z['converged'], z['Rrs']
    R = spectral.moving_mean(d.Rrs, 5, edges='shrink') if variant == 'smooth5' else d.Rrs
    fit = gsm.fit_gsm(d.wave, R, d.meta['temp'].to_numpy(), d.meta['sal'].to_numpy())
    np.savez(f, tag=tag, dRrs=fit.dRrs, params=fit.params, converged=fit.converged, Rrs=R)
    return fit.dRrs, fit.params, fit.converged, R


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    data = sdpdata.load_kramer2022()
    z = np.load(sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz')
    pcr = np.load(sdpdata.kramer2022_dir() / 'products' / 'pcr_rrsD2_1nm.npz')
    A, C, _ = load_port_coefs()
    w = data.wave
    PIGS = list(sdpdata.PIGMENTS_13)
    ACC = [p for p in PIGS if p != 'Tchla']
    lods = V.lod_proxy(data.pigments13)
    pig145 = data.pigments13
    camp = data.meta['campaign'].to_numpy()

    d = AMT.load_amt()
    if not np.array_equal(d.wave, w):
        raise SystemExit('AMT grid differs from the deposit')
    cruise = d.meta['cruise'].to_numpy()
    cruises = sorted(set(cruise))
    rlog = lambda P, p, T: np.log10(V.replace_zeros(P, lods[p]) / T)  # noqa: E731
    truth_abs = {p: d.pigments[p].to_numpy(float) for p in PIGS}
    truth_ratio = {p: rlog(truth_abs[p], p, truth_abs['Tchla']) for p in ACC}
    mask = {p: np.ones(d.n, bool) for p in PIGS}
    mask['MVchlb'] = d.flags['MVchlb_scorable'].to_numpy(bool)
    mask['DVchla'] = d.flags['DVchla_scorable'].to_numpy(bool)

    # ---- models fitted once on the 145
    chl_gsm_145 = z['params'][:, 0]
    oc4_145 = V.oc4v6(w, data.Rrs)
    p145 = z['params']
    aux_145 = np.column_stack([np.log10(p145[:, 0]), np.log10(p145[:, 1]), p145[:, 2]])
    nulls = {}
    for name, x145 in (('null_GSM', chl_gsm_145), ('null_OC4', oc4_145)):
        nulls[name] = {'abs': {}, 'ratio': {}}
        for p in PIGS:
            y = np.log10(V.replace_zeros(pig145[p].to_numpy(), lods[p]))
            nulls[name]['abs'][p] = np.polyfit(np.log10(x145), y, 1)
        for p in ACC:
            y = rlog(pig145[p].to_numpy(), p, pig145['Tchla'].to_numpy())
            nulls[name]['ratio'][p] = np.polyfit(np.log10(x145), y, 1)
    const = {p: rlog(pig145[p].to_numpy(), p, pig145['Tchla'].to_numpy()).mean() for p in ACC}
    bench = V.Benchmark(pig145, chl_gsm_145, oc4_145, camp, n_perm=2, seed=SEED)
    Sn = noise.covariance(w, np.median(data.Rrs, axis=0), 'insitu')
    W, hps, keys, tau = {}, {}, {}, {}
    for kind in ('abs', 'ratio'):
        keys[kind] = [t for t, v in bench.targets.items() if v['kind'] == kind]
        Y = np.column_stack([bench.targets[t]['y'] for t in keys[kind]])
        m = bayes.BayesianSharedW(LAM_NOISE, LAM_SMOOTH, RANKS, Sigma_n=Sn)
        m.fit(z['dRrs'], aux_145, Y, camp)
        W[kind], hps[kind] = m, m.hp
        t2, s2 = bayes.campaign_variance(m.oof_residuals(), camp)
        tau[kind] = {t.split(':')[1]: {'tau2': float(a), 'sigw2': float(b)}
                     for t, a, b in zip(keys[kind], t2, s2)}
        print(kind, 'hp', m.hp, 'tau', {k: round(np.sqrt(v['tau2']), 3)
                                        for k, v in tau[kind].items()})

    # ---- predictions per spectral variant
    P = {}
    gsm_diag = {}
    for v in VARIANTS:
        dR, pe, conv, Rv = gsm_features(d, v)
        gsm_diag[v] = {'converged': int(conv.sum()), 'n': d.n}
        _, Xd2 = spectral.difference(w, dR, 2)
        chl = np.clip(pe[:, 0], 1e-3, None)
        oc4 = V.oc4v6(w, Rv)
        aux = np.column_stack([np.log10(chl), np.log10(np.clip(pe[:, 1], 1e-4, None)), pe[:, 2]])
        pa, pr, sa, sr = {}, {}, {}, {}
        for name, getc in (('PCR', lambda p: (pcr[f'coefs_n145_{p}'],
                                               pcr[f'intercepts_n145_{p}'])),
                           ('PCR_Kramer_coefs', lambda p: (A[p], C[p]))):
            pa[name] = {p: models.predict_ensemble(Xd2, *getc(p), lod=lods[p])[0] for p in PIGS}
            pr[name + '_derived'] = {p: rlog(pa[name][p], p, V.replace_zeros(
                pa[name]['Tchla'], lods['Tchla'])) for p in ACC}
        for name, x in (('null_GSM', chl), ('null_OC4', oc4)):
            pa[name] = {p: 10**np.polyval(nulls[name]['abs'][p], np.log10(x)) for p in PIGS}
            pr[name] = {p: np.polyval(nulls[name]['ratio'][p], np.log10(x)) for p in ACC}
        pa['GSM_raw'] = {'Tchla': chl}
        pr['const'] = {p: np.full(d.n, const[p]) for p in ACC}
        for kind in ('abs', 'ratio'):
            mu, sd = W[kind].predict(dR, aux)
            tgt, sds = (pa, sa) if kind == 'abs' else (pr, sr)
            tgt['sharedW'] = {}
            for j, t in enumerate(keys[kind]):
                p = t.split(':')[1]
                tgt['sharedW'][p] = np.maximum(mu[:, j], 0.0) if kind == 'abs' else mu[:, j]
                sds[p] = sd[:, j]
        P[v] = {'abs': pa, 'ratio': pr, 'sd_abs': sa, 'sd_ratio': sr}

    def err(v, kind, model, p):
        if kind == 'abs':
            lo = np.log10(V.replace_zeros(truth_abs[p], lods[p]))
            return np.log10(V.replace_zeros(P[v]['abs'][model][p], lods[p])) - lo
        return P[v]['ratio'][model][p] - truth_ratio[p]

    # ---- scores, pooled and per cruise
    rows = []
    groups = [('pooled', np.ones(d.n, bool))] + [(c, cruise == c) for c in cruises]
    for v in VARIANTS:
        for kind in ('abs', 'ratio'):
            for model, dd in P[v][kind].items():
                for p in dd:
                    e = err(v, kind, model, p)
                    o = truth_abs[p] if kind == 'abs' else truth_ratio[p]
                    lo = np.log10(V.replace_zeros(o, lods[p])) if kind == 'abs' else o
                    for gname, gm in groups:
                        m = gm & mask[p]
                        if m.sum() < 3:
                            continue
                        row = {'variant': v, 'kind': kind, 'model': model, 'target': p,
                               'group': gname, 'n': int(m.sum()),
                               'bias_log': float(e[m].mean()),
                               'rms_log': float(np.sqrt(np.mean(e[m]**2))),
                               'crms_log': float(e[m].std()),
                               'R2_log': r2(lo[m], lo[m] + e[m]),
                               'truth_sd_log': float(lo[m].std())}
                        if kind == 'abs':
                            row['R2_lin'] = r2(o[m], dd[p][m])
                            row['n_truth_zero'] = int(np.sum(o[m] == 0))
                            row['frac_pred_below_lod'] = float(np.mean(dd[p][m] < lods[p]))
                        rows.append(row)
    S = pd.DataFrame(rows)
    S.round(4).to_csv(OUT / 'amt_holdout_scores.csv', index=False)
    # per-sample predictions (for diagnosis and for scoring other groups' models alike)
    prow = []
    for v in VARIANTS:
        for kind in ('abs', 'ratio'):
            sds = P[v]['sd_abs' if kind == 'abs' else 'sd_ratio']
            for model, dd in P[v][kind].items():
                for p, val in dd.items():
                    prow.append(pd.DataFrame({
                        'variant': v, 'sample_id': d.meta['sample_id'], 'cruise': cruise,
                        'kind': kind, 'target': p, 'model': model, 'pred': val,
                        'sd': sds[p] if model == 'sharedW' else np.nan,
                        'obs': truth_abs[p] if kind == 'abs' else truth_ratio[p],
                        'scored': mask[p]}))
    pd.concat(prow).round(6).to_csv(OUT / 'amt_holdout_predictions.csv.gz', index=False)

    # ---- ΔRMS vs best null: sample bootstrap, cluster bootstrap over cruises, per cruise
    drows = []
    # protocol §3 (AMT27 DVchla is missing since rule 7, so the earlier
    # 'DVchla_excl_AMT27' sensitivity scope is now the protocol itself)
    scopes = [('protocol', None)]
    for v, (scope, only) in [(v, sc) for v in VARIANTS for sc in scopes]:
        for kind, cands, nl in (('abs', ('PCR', 'PCR_Kramer_coefs', 'sharedW'),
                                 ('null_GSM', 'null_OC4')),
                                ('ratio', ('PCR_derived', 'PCR_Kramer_coefs_derived', 'sharedW'),
                                 ('null_GSM', 'null_OC4', 'const'))):
            for p in (PIGS if kind == 'abs' else ACC):
                if only is not None and p != only:
                    continue
                m = mask[p] & (cruise != 'AMT27') if only is not None else mask[p]
                idx = np.where(m)[0]
                ref = min(nl, key=lambda n: np.mean(err(v, kind, n, p)[m]**2))
                e0 = err(v, kind, ref, p)**2
                cl = [idx[cruise[idx] == c] for c in cruises if np.any(cruise[idx] == c)]
                bs = rng.integers(0, len(idx), size=(N_BOOT, len(idx)))
                bc = rng.integers(0, len(cl), size=(N_BOOT, len(cl)))
                for c in cands:
                    e1 = err(v, kind, c, p)**2
                    dd = np.sqrt(e1[idx][bs].mean(1)) - np.sqrt(e0[idx][bs].mean(1))
                    dc = np.empty(N_BOOT)
                    for b in range(N_BOOT):
                        ii = np.concatenate([cl[k] for k in bc[b]])
                        dc[b] = np.sqrt(e1[ii].mean()) - np.sqrt(e0[ii].mean())
                    per = [np.sqrt(e1[ci].mean()) - np.sqrt(e0[ci].mean()) for ci in cl]
                    drows.append({'variant': v, 'scope': scope, 'kind': kind, 'target': p,
                                  'model': c,
                                  'best_null': ref, 'n': int(m.sum()), 'n_cruises': len(cl),
                                  'd_rms': float(np.sqrt(e1[idx].mean())
                                                 - np.sqrt(e0[idx].mean())),
                                  'lo95_sample': float(np.percentile(dd, 2.5)),
                                  'hi95_sample': float(np.percentile(dd, 97.5)),
                                  'lo95_cluster': float(np.percentile(dc, 2.5)),
                                  'hi95_cluster': float(np.percentile(dc, 97.5)),
                                  'n_cruises_better': int(np.sum(np.array(per) < 0))})
    D = pd.DataFrame(drows)
    for k in ('sample', 'cluster'):
        D[f'verdict_{k}'] = np.where(D[f'hi95_{k}'] < 0, 'beats null',
                                     np.where(D[f'lo95_{k}'] > 0, 'worse', 'tie'))
    D.round(4).to_csv(OUT / 'amt_holdout_vs_null.csv', index=False)
    D = D[D.scope == 'protocol']

    # ---- coverage of the shared W: base vs + tau^2; cruise-mean test
    crows, mrows, pits = [], [], {}
    for v in VARIANTS:
        for kind in ('abs', 'ratio'):
            sds = P[v]['sd_abs' if kind == 'abs' else 'sd_ratio']
            mus = P[v][kind]['sharedW']
            truth = truth_abs if kind == 'abs' else truth_ratio
            for var in ('base', '+tau2'):
                ys, ms, ss, gs = [], [], [], []
                for p in sds:
                    m = mask[p]
                    s = np.sqrt(sds[p]**2 + (tau[kind][p]['tau2'] if var == '+tau2' else 0.0))
                    ys.append(truth[p][m]); ms.append(mus[p][m]); ss.append(s[m])
                    gs.append(cruise[m])
                    for gname, gm in groups:
                        mm = gm & m
                        if mm.sum() < 3:
                            continue
                        g = bayes.gaussian_scores(truth[p][mm], mus[p][mm], s[mm])
                        crows.append({'variant': v, 'kind': kind, 'variance': var, 'target': p,
                                      'group': gname, 'n': int(mm.sum()),
                                      'cov68': g['cov68'], 'cov95': g['cov95'],
                                      'z_sd': g['z_sd'], 'crps': g['crps'],
                                      'mean_z': float(np.mean((truth[p][mm] - mus[p][mm])
                                                              / s[mm]))})
                y, mu, s, gg = map(np.concatenate, (ys, ms, ss, gs))
                for gname in ['pooled'] + cruises:
                    mm = np.ones(len(y), bool) if gname == 'pooled' else gg == gname
                    g = bayes.gaussian_scores(y[mm], mu[mm], s[mm])
                    crows.append({'variant': v, 'kind': kind, 'variance': var, 'target': 'ALL',
                                  'group': gname, 'n': int(mm.sum()), 'cov68': g['cov68'],
                                  'cov95': g['cov95'], 'z_sd': g['z_sd'], 'crps': g['crps'],
                                  'mean_z': float(np.mean((y[mm] - mu[mm]) / s[mm]))})
                    if gname == 'pooled' and var == 'base':
                        pits[(v, kind)] = g['pit']
            # cruise means: independent vs campaign-aware
            for p in sds:
                t2, w2 = tau[kind][p]['tau2'], tau[kind][p]['sigw2']
                for c in cruises:
                    mm = (cruise == c) & mask[p]
                    if mm.sum() < 3:
                        continue
                    r = truth[p][mm] - mus[p][mm]
                    n = int(mm.sum())
                    s_ind = np.sqrt(np.mean(sds[p][mm]**2) / n)
                    s_cam = np.sqrt(t2 + w2 / n)
                    mrows.append({'variant': v, 'kind': kind, 'target': p, 'cruise': c, 'n': n,
                                  'mean_resid': float(r.mean()), 'sd_indep': float(s_ind),
                                  'sd_campaign': float(s_cam),
                                  'in95_indep': bool(abs(r.mean()) <= 1.96 * s_ind),
                                  'in95_campaign': bool(abs(r.mean()) <= 1.96 * s_cam)})
    Cv = pd.DataFrame(crows)
    Cv.round(4).to_csv(OUT / 'amt_holdout_coverage.csv', index=False)
    Mc = pd.DataFrame(mrows)
    Mc.round(4).to_csv(OUT / 'amt_holdout_campaign.csv', index=False)

    # ---- comparisons: EXPORTS-NA (§9) and L8 (§7)
    exv = pd.read_csv(OUT / 'exports_na_vs_null.csv')
    wt = pd.read_csv(OUT / 'weighting_table.csv')
    comp = {}
    for kind in ('abs', 'ratio'):
        for p in (PIGS if kind == 'abs' else ACC):
            row = {}
            for mname, ex_m, l8_m in (('PCR', 'PCR' if kind == 'abs' else 'PCR_derived',
                                       'PCR_dRrs_d2'),
                                      ('sharedW', 'sharedW', 'RRR_dRrs_aux_group')):
                am = D[(D.variant == 'smooth5') & (D.kind == kind) & (D.target == p)
                       & (D.model == ('PCR' if (mname == 'PCR' and kind == 'abs') else
                                      'PCR_derived' if mname == 'PCR' else 'sharedW'))]
                ex = exv[(exv.kind == kind) & (exv.target == p) & (exv.model == ex_m)]
                l8 = wt[(wt.target == f'{kind}:{p}') & (wt.model == l8_m)]
                row[mname] = {
                    'AMT_cluster': am['verdict_cluster'].iloc[0] if len(am) else None,
                    'AMT_cruises_better': int(am['n_cruises_better'].iloc[0]) if len(am) else None,
                    'AMT_n_cruises': int(am['n_cruises'].iloc[0]) if len(am) else None,
                    'EXPORTS': ex['verdict'].iloc[0] if len(ex) else None,
                    'L8_beats_null': bool(l8['beats_null_loco'].iloc[0]) if len(l8) else None}
            comp[f'{kind}:{p}'] = row

    # ---- summary
    sm = D[D.variant == 'smooth5']
    verdicts = {f'{k}:{m}': g['verdict_cluster'].value_counts().to_dict()
                for (k, m), g in sm.groupby(['kind', 'model'])}
    pooled_cov = Cv[(Cv.target == 'ALL') & (Cv.group == 'pooled')]
    git = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True,
                         text=True).stdout.strip()
    out = {'script': 'scripts/sdp/amt_holdout.py', 'epft_up_version': __version__,
           'git_head': git, 'protocol': 'reports/provenance/AMT_benchmark_protocol.md',
           'amt_input': d.provenance, 'deposit_sha256': data.provenance['sha256'],
           'n_spectra': d.n, 'n_per_cruise': pd.Series(cruise).value_counts().to_dict(),
           'mvchlb_scored_n': int(mask['MVchlb'].sum()), 'gsm': gsm_diag,
           'sharedW_hp': hps, 'tau_deposit': tau, 'n_boot': N_BOOT, 'seed': SEED,
           'verdict_counts_cluster_smooth5': verdicts,
           'coverage_pooled': pooled_cov.drop(columns=['target', 'group']).round(4)
           .to_dict(orient='records'),
           'campaign_mean_test': Mc.groupby(['variant', 'kind'])[['in95_indep', 'in95_campaign']]
           .mean().round(3).reset_index().to_dict(orient='records'),
           'between_cruise_sd_observed_vs_tau': {
               kind: {p: {'tau_deposit': float(np.sqrt(tau[kind][p]['tau2'])),
                          'sd_cruise_means_amt': float(
                              Mc[(Mc.variant == 'smooth5') & (Mc.kind == kind) &
                                 (Mc.target == p)]['mean_resid'].std()),
                          'amt_offset': float(
                              Mc[(Mc.variant == 'smooth5') & (Mc.kind == kind) &
                                 (Mc.target == p)]['mean_resid'].mean())}
                      for p in tau[kind]} for kind in tau},
           'comparison_AMT_EXPORTS_L8': comp,
           'dvchla_scored_n': int(mask['DVchla'].sum()),
           'amt27_dvchla_note': 'AMT27 (DHI) reports DVchla = 0 in all samples while MVchla/'
                                'Tchla is 0.34-0.98: DVchla not separated. Protocol rule 7 '
                                '(2026-10-24) makes it missing; DVchla is scored on 7 cruises.',
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    (OUT / 'amt_holdout_summary.json').write_text(json.dumps(out, indent=1, default=str))

    # ---- figures
    figures(d, P, S, D, Cv, Mc, pits, truth_abs, truth_ratio, cruises, cruise)

    pd.set_option('display.width', 250)
    print(sm.pivot_table(index=['kind', 'target'], columns='model',
                         values=['d_rms', 'n_cruises_better'], aggfunc='first').round(2)
          .to_string())
    print(sm[['kind', 'target', 'model', 'best_null', 'verdict_sample', 'verdict_cluster']]
          .to_string())

    print(pooled_cov.round(3).to_string())
    print(json.dumps(out['campaign_mean_test'], indent=0))
    print(json.dumps(verdicts, indent=0))


def figures(d, P, S, D, Cv, Mc, pits, truth_abs, truth_ratio, cruises, cruise):
    cm = plt.get_cmap('tab10')
    ccol = {c: cm(i) for i, c in enumerate(cruises)}
    v = 'smooth5'
    # Tchla
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    for ax, ks in ((axs[0], ('PCR', 'PCR_Kramer_coefs')), (axs[1], ('null_GSM', 'null_OC4',
                                                                     'sharedW'))):
        for k in ks:
            y = np.clip(P[v]['abs'][k]['Tchla'], 5e-3, None)
            ax.loglog(truth_abs['Tchla'], y, 'o', ms=3.5, alpha=0.7, label=k)
        ax.plot([5e-3, 5], [5e-3, 5], 'k-', lw=0.6)
        ax.set_xlabel(r'HPLC Tchla [mg m$^{-3}$]')
        ax.legend(frameon=False, fontsize=9)
    sel = S[(S.variant == v) & (S.kind == 'abs') & (S.target == 'Tchla') & (S.group == 'pooled')]
    txt = ', '.join(f"{r.model} {r.bias_log:+.2f}/{r.rms_log:.2f}" for r in sel.itertuples()
                    if r.model in ('PCR', 'null_GSM', 'sharedW'))
    axs[0].set_ylabel(r'predicted Tchla [mg m$^{-3}$]')
    fig.suptitle(f'AMT hold-out (182 spectra, 8 cruises), Tchla. bias/RMS dex: {txt}',
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / 'amt_holdout_tchla.png', dpi=200)
    plt.close(fig)

    # key ratios
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5))
    for ax, p in zip(axes.ravel(), KEY_RATIOS):
        o = truth_ratio[p]
        best = D[(D.variant == v) & (D.kind == 'ratio') & (D.target == p)]['best_null'].iloc[0]
        ax.plot(o, P[v]['ratio']['PCR_derived'][p], 's', ms=3.5, alpha=0.6, label='PCR (derived)')
        ax.plot(o, P[v]['ratio'][best][p], '^', ms=3.5, alpha=0.6, label=f'best null ({best})')
        ax.errorbar(o, P[v]['ratio']['sharedW'][p], yerr=P[v]['sd_ratio'][p], fmt='o', ms=3.5,
                    lw=0.5, alpha=0.7, label='shared W ±1σ')
        lim = [min(o.min(), -3) - 0.1, max(o.max(), 0) + 0.1]
        ax.plot(lim, lim, 'k-', lw=0.6)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_title(f'log10 {p}:Tchla', fontsize=11)
        ax.set_xlabel('HPLC')
        ax.set_ylabel('predicted')
        ax.legend(fontsize=7, frameon=False)
    ax = axes.ravel()[-1]
    sub = D[(D.variant == v) & (D.kind == 'ratio')]
    labels = list(sub['target'].unique())
    for k, (mname, mk) in enumerate((('PCR_derived', 's'), ('sharedW', 'o'))):
        q = sub[sub.model == mname].set_index('target').loc[labels]
        y = np.arange(len(labels)) + (k - 0.5) * 0.25
        ax.errorbar(q['d_rms'], y, xerr=[q['d_rms'] - q['lo95_cluster'],
                                          q['hi95_cluster'] - q['d_rms']],
                    fmt=mk, ms=4, lw=0.8, label=mname)
    ax.axvline(0, color='k', lw=0.6)
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('ΔRMS vs best null [dex] (cluster 95% CI)')
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle('AMT hold-out: pigment:Tchla ratios (5 nm-smoothed spectra; no retraining)')
    fig.tight_layout()
    fig.savefig(OUT / 'amt_holdout_ratios.png', dpi=200)
    plt.close(fig)

    # per-cruise bias
    tg = [('abs', 'Tchla')] + [('ratio', p) for p in KEY_RATIOS]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    x = np.arange(len(cruises))
    for ax, (kind, p) in zip(axes.ravel(), tg):
        best = D[(D.variant == v) & (D.kind == kind) & (D.target == p)]['best_null'].iloc[0]
        for k, (mname, mk) in enumerate((('PCR' if kind == 'abs' else 'PCR_derived', 's'),
                                         (best, '^'), ('sharedW', 'o'))):
            q = S[(S.variant == v) & (S.kind == kind) & (S.model == mname) & (S.target == p)]
            q = q.set_index('group').reindex(cruises)
            ax.plot(x + (k - 1) * 0.15, q['bias_log'], mk, ms=6, label=mname)
            pooled = S[(S.variant == v) & (S.kind == kind) & (S.model == mname) &
                       (S.target == p) & (S.group == 'pooled')]['bias_log'].iloc[0]
            ax.axhline(pooled, ls=':', lw=0.8, color=ax.lines[-1].get_color())
        ax.axhline(0, color='k', lw=0.6)
        ax.set_title(f"{'Tchla' if kind == 'abs' else p + ':Tchla'} bias by cruise", fontsize=10)
        ax.set_xticks(x, [c.replace('AMT', '') for c in cruises])
        ax.set_ylabel('bias [dex]')
        ax.legend(fontsize=7, frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel('AMT cruise')
    fig.suptitle('AMT hold-out: per-cruise bias (dotted: pooled)')
    fig.tight_layout()
    fig.savefig(OUT / 'amt_holdout_cruise_bias.png', dpi=200)
    plt.close(fig)

    # coverage
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, kind in zip(axes[:2], ('ratio', 'abs')):
        for k, (var, col) in enumerate((('base', 'C0'), ('+tau2', 'C1'))):
            q = Cv[(Cv.variant == v) & (Cv.kind == kind) & (Cv.variance == var) &
                   (Cv.target == 'ALL')].set_index('group')
            g = ['pooled'] + cruises
            xx = np.arange(len(g))
            ax.bar(xx + (k - 0.5) * 0.38, q.loc[g, 'cov68'], 0.38, color=col, alpha=0.5,
                   label=f'68% ({var})')
            ax.plot(xx + (k - 0.5) * 0.38, q.loc[g, 'cov95'], 'v', color=col, ms=7,
                    label=f'95% ({var})')
        ax.axhline(0.68, color='k', ls=':', lw=0.8)
        ax.axhline(0.95, color='k', ls='--', lw=0.8)
        ax.set_xticks(xx, [s.replace('AMT', '') for s in g], rotation=0)
        ax.set_ylim(0, 1.05)
        ax.set_title(f'shared W interval coverage, {kind}', fontsize=11)
        ax.legend(fontsize=7, frameon=False, ncol=2)
    ax = axes[2]
    ax.hist(pits[(v, 'ratio')], bins=10, range=(0, 1), density=True, alpha=0.7,
            label='ratios (base)')
    ax.hist(pits[(v, 'abs')], bins=10, range=(0, 1), density=True, histtype='step', lw=2,
            label='absolute (base)')
    ax.axhline(1, color='k', lw=0.6)
    q = Mc[(Mc.variant == v) & (Mc.kind == 'ratio')]
    ax.set_title(f"PIT (pooled). Ratio cruise means inside 95%:\n"
                 f"independent errors {q['in95_indep'].mean():.2f}, "
                 f"campaign-aware {q['in95_campaign'].mean():.2f}", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / 'amt_holdout_coverage.png', dpi=200)
    plt.close(fig)

    # smoothed vs as-is
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, kind in zip(axes, ('abs', 'ratio')):
        for mname, mk in ((('PCR' if kind == 'abs' else 'PCR_derived'), 's'), ('sharedW', 'o'),
                          ('null_OC4', '^')):
            q = S[(S.kind == kind) & (S.model == mname) & (S.group == 'pooled')]
            a = q[q.variant == 'smooth5'].set_index('target')['rms_log']
            b = q[q.variant == 'asis'].set_index('target')['rms_log'].reindex(a.index)
            ax.plot(a, b, mk, ms=6, label=mname)
        lim = [0, max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lim, lim, 'k-', lw=0.6)
        ax.set_xlabel('RMS [dex], 5 nm-smoothed spectra')
        ax.set_ylabel('RMS [dex], as-delivered spectra')
        ax.set_title(f'{kind}: each point is one target', fontsize=11)
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / 'amt_holdout_smooth_vs_asis.png', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    main()
