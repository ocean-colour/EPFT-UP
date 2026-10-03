"""Execution #4: skill beyond Tchla, composition (ratios), and leakage (LOCO).

For every pigment, on identical splits, compare Kramer's PCR on δRrs'' with:

* ``null_GSM``: log10 P = a + b log10 Tchla_GSM (Tchla from the Execution #2 fit);
* ``null_OC4``: the same with OC4v6 band-ratio Tchla (no hyperspectral info);
* ``oracle``: the same with HPLC Tchla (Tchla-covariance upper bound; not
  achievable from optics);
* ``const``: the training mean (for ratios).

Targets:

* absolute concentrations, scored like Table 2 (linear) and in log10;
* pigment:Tchla ratios. PCR is trained (i) on log10 ratio (zeros →
  ½·LOD) and (ii) on the linear ratio; (iii) is the ratio of PCR's
  absolute predictions (what a user of SDP products would compute). All are
  scored in log10 and linear space.

Schemes: 100 random 75/25 splits (Kramer) and leave-one-campaign-out
(LOCO, 8 campaigns, pooled predictions with a paired bootstrap).
A zero-replacement sensitivity (frac = 0.1, 0.5, 1.0 of LOD) is run for
the log-ratio targets under LOCO.

Primary sample: all 145 spectra (as the paper). The absolute-concentration
comparison is repeated without the degenerate SABOR GSM fit (N = 144).

Outputs (``reports/figures/sdp/``): ``diag_skill_absolute.png``,
``diag_skill_ratio.png``, ``diag_loco_tchla.png``, ``diag_skill_table.csv``,
``diag_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/diagnostics_null_loco.py
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
from epft_up.sdp import data as sdpdata, spectral, validate as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
PIGS = list(sdpdata.PIGMENTS_13)
ACC = [p for p in PIGS if p != 'Tchla']
SEED = 1
N_PERM = 100
ZERO_FRAC = 0.5
SENS_FRACS = (0.1, 0.5, 1.0)


def run_scheme(fitters, splits, scheme, n):
    """Cross-validate each fitter; return {name: results}."""
    return {name: V.cross_validate(f, splits) for name, f in fitters.items()}


def summarize(res, obs, scorer, scheme, n, null_names, model='PCR'):
    """Per-model scores and paired model-vs-null comparisons for one scheme."""
    out = {}
    if scheme == 'random':
        per = {m: V.per_split_scores(r, obs, scorer) for m, r in res.items()}
        for m, s in per.items():
            out[m] = {k: [float(np.nanmean(v)), float(np.nanstd(v, ddof=1))]
                      for k, v in s.items()}
        for nn in null_names:
            out[f'{model}_vs_{nn}'] = V.paired(per[model], per[nn])
    else:
        pooled = {m: V.pooled_predictions(r, n) for m, r in res.items()}
        for m, pr in pooled.items():
            out[m] = scorer(obs, pr)
        for nn in null_names:
            out[f'{model}_vs_{nn}'] = V.bootstrap_paired(obs, pooled[model],
                                                         pooled[nn], scorer)
    return out


def absolute_block(X, pig, tchla_gsm, tchla_oc4, lods, splits_by_scheme, n):
    """Absolute concentrations: PCR vs Tchla nulls, linear and log scoring."""
    out = {}
    hplc_tchla = pig['Tchla'].to_numpy()
    for p in PIGS:
        y = pig[p].to_numpy()
        ylog = np.log10(V.replace_zeros(y, lods[p], ZERO_FRAC))
        fitters = {'PCR': V.PCRFitter(X, y, 'pigment', seed=SEED),
                   'null_GSM': V.LogLinearFitter(tchla_gsm, ylog),
                   'null_OC4': V.LogLinearFitter(tchla_oc4, ylog)}
        if p != 'Tchla':
            fitters['oracle'] = V.LogLinearFitter(hplc_tchla, ylog)
        nulls = [m for m in fitters if m != 'PCR']
        out[p] = {}
        for scheme, splits in splits_by_scheme.items():
            res = run_scheme(fitters, splits, scheme, n)
            lin = summarize(res, y, V.score_linear, scheme, n, nulls)
            # log scoring: zeros in obs and predictions -> frac * LOD
            res_log = {m: [(te, np.log10(V.replace_zeros(pr, lods[p], ZERO_FRAC)))
                           for te, pr in r] for m, r in res.items()}
            log = summarize(res_log, ylog, V.score_log, scheme, n, nulls)
            out[p][scheme] = {'linear': lin, 'log': log}
            if scheme == 'random' and p in ('Tchla', 'Fuco'):
                print(f'  {p} {scheme} linear R2:',
                      {m: round(lin[m]['R2'][0], 3) for m in fitters})
    return out


def ratio_block(X, pig, tchla_gsm, tchla_oc4, lods, splits_by_scheme, n, abs_preds):
    """Pigment:Tchla ratios: PCR (log, linear, derived) vs nulls."""
    out = {}
    tchla = pig['Tchla'].to_numpy()
    for p in ACC:
        r_lin = pig[p].to_numpy() / tchla
        r_log = np.log10(V.replace_zeros(pig[p].to_numpy(), lods[p], ZERO_FRAC) / tchla)
        fitters_log = {'PCR_log': V.PCRFitter(X, r_log, None, seed=SEED),
                       'null_GSM': V.LogLinearFitter(tchla_gsm, r_log, output='log'),
                       'null_OC4': V.LogLinearFitter(tchla_oc4, r_log, output='log'),
                       'const': V.ConstantFitter(r_log)}
        out[p] = {}
        for scheme, splits in splits_by_scheme.items():
            res = run_scheme(fitters_log, splits, scheme, n)
            # derived ratio from PCR's absolute predictions on the same splits
            res['PCR_derived'] = [
                (te, np.log10(V.replace_zeros(pp, lods[p], ZERO_FRAC)
                              / V.replace_zeros(pt, lods['Tchla'], ZERO_FRAC)))
                for (te, pp), (_, pt) in zip(abs_preds[scheme][p], abs_preds[scheme]['Tchla'])]
            res['PCR_lin'] = V.cross_validate(V.PCRFitter(X, r_lin, 'pigment', seed=SEED),
                                              splits)
            res_log = {m: r for m, r in res.items() if m != 'PCR_lin'}
            res_log['PCR_lin'] = [(te, np.log10(V.replace_zeros(pr * tchla[te], lods[p],
                                                                ZERO_FRAC) / tchla[te]))
                                  for te, pr in res['PCR_lin']]
            nulls = ['null_GSM', 'null_OC4', 'const', 'PCR_derived', 'PCR_lin']
            log = summarize(res_log, r_log, V.score_log, scheme, n, nulls, model='PCR_log')
            res_linspace = {m: [(te, 10.0**pr) for te, pr in r] for m, r in res_log.items()}
            lin = summarize(res_linspace, r_lin, V.score_linear, scheme, n,
                            ['null_GSM', 'null_OC4', 'const'], model='PCR_log')
            out[p][scheme] = {'log': log, 'linear': lin}
    return out


def zero_sensitivity(X, pig, tchla_gsm, lods, loco, n):
    """LOCO pooled log-ratio R² of PCR_log and null_GSM for several zero fracs."""
    out = {}
    tchla = pig['Tchla'].to_numpy()
    for p in ACC:
        out[p] = {}
        for frac in SENS_FRACS:
            r_log = np.log10(V.replace_zeros(pig[p].to_numpy(), lods[p], frac) / tchla)
            pcr = V.pooled_predictions(V.cross_validate(
                V.PCRFitter(X, r_log, None, seed=SEED), loco), n)
            nul = V.pooled_predictions(V.cross_validate(
                V.LogLinearFitter(tchla_gsm, r_log, output='log'), loco), n)
            out[p][str(frac)] = {'PCR_log_R2': V.score_log(r_log, pcr)['R2'],
                                 'null_GSM_R2': V.score_log(r_log, nul)['R2'],
                                 'PCR_log_RMS': V.score_log(r_log, pcr)['RMS'],
                                 'null_GSM_RMS': V.score_log(r_log, nul)['RMS']}
    return out


def collect_abs_preds(X, pig, splits_by_scheme):
    """PCR absolute predictions on each scheme's splits (for derived ratios)."""
    return {scheme: {p: V.cross_validate(V.PCRFitter(X, pig[p].to_numpy(), 'pigment',
                                                     seed=SEED), splits)
                     for p in PIGS}
            for scheme, splits in splits_by_scheme.items()}


def verdicts(abs_res, rat_res):
    """Per-pigment yes/no: skill beyond the best Tchla null, by target and scheme.

    Random: PCR beats the null in ≥ 90% of paired splits on R².
    LOCO: the bootstrap 95% interval of ΔR² excludes 0.
    The best null is the one (GSM or OC4) with the higher R² in that scheme.
    """
    rows = []
    for p in PIGS:
        row = {'pigment': p}
        for scheme in ('random', 'loco'):
            for space in ('linear', 'log'):
                blk = abs_res[p][scheme][space]

                def r2(m):
                    v = blk[m]['R2']
                    return v[0] if isinstance(v, list) else v

                best = max(('null_GSM', 'null_OC4'), key=r2)
                cmp = blk[f'PCR_vs_{best}']
                if scheme == 'random':
                    ok = cmp['R2']['frac_model_better'] >= 0.9
                    d = cmp['R2']['mean_diff']
                else:
                    ok = cmp['dR2_ci'][0] > 0
                    d = cmp['dR2']
                key = f'abs_{space}_{scheme}'
                row[f'{key}_R2_PCR'] = r2('PCR')
                row[f'{key}_R2_null'] = r2(best)
                row[f'{key}_dR2'] = d
                row[f'{key}_beyond'] = bool(ok)
                if 'oracle' in blk:
                    row[f'{key}_R2_oracle'] = r2('oracle')
        if p in rat_res:
            for scheme in ('random', 'loco'):
                blk = rat_res[p][scheme]['log']

                def r2(m):
                    v = blk[m]['R2']
                    return v[0] if isinstance(v, list) else v

                best = max(('null_GSM', 'null_OC4'), key=r2)
                cmp = blk[f'PCR_log_vs_{best}']
                if scheme == 'random':
                    ok = cmp['R2']['frac_model_better'] >= 0.9
                    d = cmp['R2']['mean_diff']
                else:
                    ok = cmp['dR2_ci'][0] > 0
                    d = cmp['dR2']
                key = f'ratio_log_{scheme}'
                row[f'{key}_R2_PCR'] = r2('PCR_log')
                row[f'{key}_R2_null'] = r2(best)
                row[f'{key}_R2_derived'] = r2('PCR_derived')
                row[f'{key}_dR2'] = d
                row[f'{key}_beyond'] = bool(ok)
        rows.append(row)
    return pd.DataFrame(rows).set_index('pigment')


def fig_absolute(abs_res, path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharey=True)
    models_ = [('PCR', 'k', 'o'), ('null_GSM', 'C0', 's'), ('null_OC4', 'C1', '^'),
               ('oracle', 'C2', 'D')]
    x = np.arange(len(PIGS))
    for i, scheme in enumerate(('random', 'loco')):
        for j, space in enumerate(('linear', 'log')):
            ax = axes[j, i]
            for k, (m, c, mk) in enumerate(models_):
                vals = []
                for p in PIGS:
                    blk = abs_res[p][scheme][space]
                    v = blk.get(m, {}).get('R2', np.nan)
                    vals.append(v[0] if isinstance(v, list) else v)
                ax.plot(x + (k - 1.5) * 0.12, vals, mk, color=c, label=m, ms=5,
                        mfc='none' if m == 'oracle' else c)
            ax.set_xticks(x)
            ax.set_xticklabels(PIGS, rotation=60, fontsize=8)
            ax.set_ylim(0, 1)
            ax.grid(axis='y', lw=0.3)
            ax.set_title(f'absolute, {space} scoring, {scheme}', fontsize=9)
    axes[0, 0].set_ylabel('R²')
    axes[1, 0].set_ylabel('R²')
    axes[0, 0].legend(fontsize=7.5, frameon=False, ncol=2)
    fig.suptitle('PCR on δRrs″ vs Tchla-only null models (N = 145)', fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_ratio(rat_res, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    models_ = [('PCR_log', 'k', 'o'), ('PCR_derived', '0.5', 'x'),
               ('null_GSM', 'C0', 's'), ('null_OC4', 'C1', '^')]
    x = np.arange(len(ACC))
    for ax, scheme in zip(axes, ('random', 'loco')):
        for k, (m, c, mk) in enumerate(models_):
            vals = []
            for p in ACC:
                v = rat_res[p][scheme]['log'][m]['R2']
                vals.append(v[0] if isinstance(v, list) else v)
            ax.plot(x + (k - 1.5) * 0.12, vals, mk, color=c, label=m, ms=5)
        ax.set_xticks(x)
        ax.set_xticklabels([f'{p}:Tchla' for p in ACC], rotation=60, fontsize=8)
        ax.set_ylim(0, 1)
        ax.grid(axis='y', lw=0.3)
        ax.set_title(f'log10 pigment:Tchla ratio, {scheme}', fontsize=9)
    axes[0].set_ylabel('R² (log10 space)')
    axes[0].legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_loco_tchla(abs_preds, pig, meta, tchla_gsm, loco_names, path):
    """Held-out-campaign Tchla: PCR vs GSM-null predictions."""
    n = len(tchla_gsm)
    pcr = V.pooled_predictions(abs_preds['loco']['Tchla'], n)
    obs = pig['Tchla'].to_numpy()
    fig, ax = plt.subplots(figsize=(5.5, 5))
    colors = dict(zip(loco_names, plt.cm.tab10(np.arange(len(loco_names)))))
    for c in loco_names:
        s = (meta['campaign'] == c).to_numpy()
        ax.scatter(obs[s], np.maximum(pcr[s], 1e-3), s=14, color=colors[c], label=c)
    ax.plot([1e-2, 6], [1e-2, 6], 'k-', lw=0.7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'HPLC Tchla [mg m$^{-3}$]')
    ax.set_ylabel(r'PCR Tchla, campaign held out [mg m$^{-3}$] (≤0 → 1e-3)')
    ax.legend(fontsize=7, frameon=False)
    ax.set_title('Leave-one-campaign-out PCR Tchla', fontsize=9)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    gsm_prod = sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz'
    z = np.load(gsm_prod)
    _, X = spectral.difference(z['wave'], z['dRrs'], order=2)
    pig = data.pigments13
    lods = V.lod_proxy(pig)
    tchla_gsm = z['params'][:, 0]
    tchla_oc4 = V.oc4v6(data.wave, data.Rrs)
    n = data.n
    loco_names, loco = V.loco_splits(data.meta['campaign'].to_numpy())
    splits = {'random': V.random_splits(n, N_PERM, 0.75, SEED), 'loco': loco}

    summary = {'script': 'scripts/sdp/diagnostics_null_loco.py',
               'epft_up_version': __version__, 'input': data.provenance,
               'gsm_product': {'file': str(gsm_prod), 'sha256': sdpdata.sha256_of(gsm_prod)},
               'lod_proxy': lods, 'zero_frac': ZERO_FRAC, 'seed': SEED, 'n_perm': N_PERM,
               'loco_campaigns': {c: int((data.meta['campaign'] == c).sum())
                                  for c in loco_names},
               'run_utc': datetime.datetime.now(datetime.timezone.utc)
               .isoformat(timespec='seconds')}

    print('absolute concentrations (N=145) ...')
    abs_res = absolute_block(X, pig, tchla_gsm, tchla_oc4, lods, splits, n)
    print('PCR absolute predictions for derived ratios ...')
    abs_preds = collect_abs_preds(X, pig, splits)
    print('ratios ...')
    rat_res = ratio_block(X, pig, tchla_gsm, tchla_oc4, lods, splits, n, abs_preds)
    print('zero-replacement sensitivity (LOCO) ...')
    sens = zero_sensitivity(X, pig, tchla_gsm, lods, loco, n)

    print('absolute concentrations without the degenerate SABOR fit (N=144) ...')
    keep = np.ones(n, dtype=bool)
    keep[int(np.argmax(np.abs(np.log10(tchla_gsm / pig['Tchla'].to_numpy()))))] = False
    pig144 = pig.loc[keep].reset_index(drop=True)
    n144 = int(keep.sum())
    loco_names144, loco144 = V.loco_splits(data.meta['campaign'].to_numpy()[keep])
    splits144 = {'random': V.random_splits(n144, N_PERM, 0.75, SEED), 'loco': loco144}
    abs_res144 = absolute_block(X[keep], pig144, tchla_gsm[keep], tchla_oc4[keep], lods,
                                splits144, n144)

    table = verdicts(abs_res, rat_res)
    table144 = verdicts(abs_res144, rat_res)   # ratio columns from N=145
    cols_show = [c for c in table.columns if c.endswith(('_R2_PCR', '_R2_null', '_beyond'))
                 and ('linear' in c or 'ratio' in c)]
    pd.set_option('display.width', 250)
    print(table[cols_show].round(2).to_string())
    table.round(4).to_csv(OUT / 'diag_skill_table.csv')

    summary.update({'absolute': abs_res, 'absolute_n144': abs_res144, 'ratio': rat_res,
                    'zero_sensitivity_loco': sens,
                    'verdicts_n145': json.loads(table.to_json(orient='index')),
                    'verdicts_n144_absolute': json.loads(
                        table144[[c for c in table144.columns if c.startswith('abs')]]
                        .to_json(orient='index'))})
    fig_absolute(abs_res, OUT / 'diag_skill_absolute.png')
    fig_ratio(rat_res, OUT / 'diag_skill_ratio.png')
    fig_loco_tchla(abs_preds, pig, data.meta, tchla_gsm, loco_names,
                   OUT / 'diag_loco_tchla.png')
    with open(OUT / 'diag_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f'wrote {OUT / "diag_summary.json"}')


if __name__ == '__main__':
    main()
