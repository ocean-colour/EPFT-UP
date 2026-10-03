"""Execution #8: Bayesian shared W: predictive intervals, input-noise propagation, calibration.

Model: :class:`epft_up.sdp.bayes.BayesianSharedW` (Execution #7's shared
low-rank W on δRrs ⊕ GSM parameters, with hyperparameters by inner
leave-campaign-out, plus per-target evidence-maximized Bayesian regression on
its latent scores). It is trained on the deposit's (clean) spectra and
evaluated under leave-one-campaign-out (LOCO).

Scenarios (test spectra; ``K`` Monte Carlo noise draws each):

* ``clean_1nm``: the deposit as is;
* ``insitu_1nm``: + in-situ-like ``pct:0.02`` (+ floor), through the 5 nm mean;
* ``pace_white_1nm``: + PACE OCI per-band σ (ocpy), white, through the 5 nm mean;
* ``pace_corr_1nm``: + smooth PACE-like error (ℓ = 30 nm) + 10% white;
* ``clean_5nm``, ``pace_white_5nm``: everything (GSM residual, model) at 5 nm
  sampling (400, 405, …, 700 nm).

Every noisy draw is propagated through the whole chain: noisy Rrs → GSM
refit → δRrs and GSM parameters → prediction. Each retrieval's interval
combines the model's predictive variance with the input-noise variance,
estimated as the across-draw variance of its predictive mean.

Variants:

* ``unaware``: trained exactly as in #7 (cross-fitted calibration, the default);
* ``unaware_insample``: the same with in-sample calibration (shows its optimism);
* ``aware``: the same, with the scenario's own noise covariance added untuned
  to the penalty (n · wᵀ Σ_test w).

References: Kramer's PCR on δRrs'' (trained clean, applied to the same noisy
draws) and the Tchla-only nulls (trained on clean GSM/OC4 Tchla, applied to
the noisy-draw Tchla).

Scores, pooled over the LOCO folds and draws: R² (log for ratios, linear for
concentrations) and ΔR² vs the best null; 68/95% coverage; the SD of the
standardized errors; CRPS; PIT histograms.

Outputs (``reports/figures/sdp/``): ``uncertainty_scores.csv``,
``uncertainty_summary.csv``, ``uncertainty_coverage.png``,
``uncertainty_pit.png``, ``uncertainty_skill.png``,
``uncertainty_examples.png``, ``uncertainty_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/uncertainty.py
"""
import datetime
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import bayes, data as sdpdata, gsm, models, noise, spectral, theory  # noqa: E402
from epft_up.sdp import validate as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
SEED = 1
K = int(os.environ.get("SDP_K", 30))   # Monte Carlo draws (override for quick tests)
LAM_NOISE = (0.0, 0.1, 1.0, 10.0, 100.0)
LAM_SMOOTH = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0)
RANKS = tuple(range(1, 11))
SCENARIOS = {'clean_1nm': (1, None), 'insitu_1nm': (1, 'insitu'),
             'pace_white_1nm': (1, 'pace_white'), 'pace_corr_1nm': (1, 'pace_corr'),
             'clean_5nm': (5, None), 'pace_white_5nm': (5, 'pace_white')}
KEY_RATIOS = ['ratio:Fuco', 'ratio:Zea', 'ratio:DVchla', 'ratio:Chlc12']


# ------------------------------------------------------------------ features
def gsm_features(wave, Rrs, T, S, tables):
    """δRrs, GSM auxiliaries (log Tchla, log a_dg443, b_bp443) and Tchla_GSM.

    Non-converged Nelder-Mead fits (possible on noisy spectra) are refitted
    with the bounded optimizer; Tchla and a_dg are floored for the logs.
    """
    fit = gsm.fit_gsm(wave, Rrs, T, S, tables=tables)
    bad = ~fit.converged | ~np.isfinite(fit.params).all(axis=1)
    params, dRrs = fit.params.copy(), fit.dRrs.copy()
    if bad.any():
        f2 = gsm.fit_gsm(wave, Rrs[bad], T[bad], S[bad], tables=tables, method='lsq')
        params[bad], dRrs[bad] = f2.params, f2.dRrs
    chl = np.clip(params[:, 0], 1e-3, None)
    adg = np.clip(params[:, 1], 1e-4, None)
    aux = np.column_stack([np.log10(chl), np.log10(adg), params[:, 2]])
    return dRrs, aux, chl, int(bad.sum())


def oc4_any_grid(wave, Rrs):
    """OC4v6 on any grid (linear interpolation to the needed bands).

    For noisy spectra the max-band ratio is clipped to [0.01, 100] so a
    non-positive Rrs(555) degrades the retrieval instead of returning NaN
    (inactive on the deposit).
    """
    w1 = np.arange(400.0, 701.0)
    R1 = np.vstack([np.interp(w1, wave, r) for r in np.atleast_2d(Rrs)])
    at = lambda w: R1[:, int(np.searchsorted(w1, w))]  # noqa: E731
    num = np.maximum.reduce([at(443), at(490), at(510)])
    den = at(555)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(den > 0, num / den, 100.0)
    x = np.log10(np.clip(ratio, 0.01, 100.0))
    return 10.0**np.polyval(V.OC4V6[::-1], x)


def build_scenarios(data, rng):
    """Clean and noisy-draw features for every scenario."""
    T, S = data.meta['temp'].to_numpy(), data.meta['sal'].to_numpy()
    feats, info = {}, {}
    for res in (1, 5):
        wave = np.arange(400.0, 701.0, float(res))
        Rrs = data.Rrs[:, ::res]
        tables = gsm.load_ref_tables(wave)
        d, a, c, nb = gsm_features(wave, Rrs, T, S, tables)
        feats[(res, None)] = {'wave': wave, 'dRrs': d[None], 'aux': a[None], 'chl': c[None],
                              'oc4': oc4_any_grid(wave, Rrs)[None], 'Rrs': Rrs}
        info[f'{res}nm_clean_refits'] = nb
        for name, (r, model) in SCENARIOS.items():
            if r != res or model is None:
                continue
            t0 = time.time()
            draws = noise.draw(wave, Rrs, model, K, rng)
            D, A, C, O, nbad = [], [], [], [], 0
            for k in range(K):
                d, a, c, nb = gsm_features(wave, draws[k], T, S, tables)
                D.append(d), A.append(a), C.append(c), O.append(oc4_any_grid(wave, draws[k]))
                nbad += nb
            feats[(res, model)] = {'wave': wave, 'dRrs': np.array(D), 'aux': np.array(A),
                                   'chl': np.array(C), 'oc4': np.array(O), 'Rrs': Rrs}
            info[f'{name}_refits'] = nbad
            print(f'  features {name}: {time.time() - t0:.0f}s (bounded refits: {nbad})',
                  flush=True)
    return feats, info


# --------------------------------------------------------------------- models
def targets(bench, kind):
    keys = [t for t, v in bench.targets.items() if v['kind'] == kind]
    return keys, np.column_stack([bench.targets[t]['y'] for t in keys])


def run_loco(data, bench, feats):
    """LOCO predictions for every scenario, model and target."""
    camp = data.meta['campaign'].to_numpy()
    names, loco = V.loco_splits(camp)
    kinds = {k: targets(bench, k) for k in ('ratio', 'abs')}
    preds = {}       # (scenario, model, target) -> (mu[K, n], sd[K, n])
    rng = np.random.default_rng(SEED)
    for fold, (tr, te) in zip(names, loco):
        t0 = time.time()
        for res in (1, 5):
            clean = feats[(res, None)]
            wave = clean['wave']
            Sn_shape = noise.covariance(wave, np.median(clean['Rrs'], axis=0), 'insitu')
            D2 = theory.difference_matrix(len(wave), 2)
            for kind, (keys, Y) in kinds.items():
                trained = {}
                unaware = bayes.BayesianSharedW(LAM_NOISE, LAM_SMOOTH, RANKS, Sigma_n=Sn_shape)
                unaware.fit(clean['dRrs'][0][tr], clean['aux'][0][tr], Y[tr], camp[tr])
                trained['unaware'] = unaware
                ins = bayes.BayesianSharedW(LAM_NOISE, LAM_SMOOTH, RANKS, Sigma_n=Sn_shape,
                                            calibration='insample')
                ins.fit(clean['dRrs'][0][tr], clean['aux'][0][tr], Y[tr], camp[tr])
                trained['unaware_insample'] = ins
                for name, (r, model) in SCENARIOS.items():
                    if r != res or model is None:
                        continue
                    Sab = noise.covariance(wave, np.median(clean['Rrs'], axis=0), model)
                    aw = bayes.BayesianSharedW(LAM_NOISE, LAM_SMOOTH, RANKS, Sigma_n=Sn_shape,
                                               Sigma_abs=Sab)
                    aw.fit(clean['dRrs'][0][tr], clean['aux'][0][tr], Y[tr], camp[tr])
                    trained[f'aware:{model}'] = aw
                # PCR reference (paper), trained clean on this resolution's δRrs''
                Xd2_tr = clean['dRrs'][0][tr] @ D2.T
                pcr = {}
                for j, t in enumerate(keys):
                    coef, icpt, _ = models.fit_pcr_one(
                        Xd2_tr, Y[tr, j], rng,
                        constraint='pigment' if kind == 'abs' else None)
                    pcr[t] = (coef, icpt)
                for name, (r, model) in SCENARIOS.items():
                    if r != res:
                        continue
                    F = feats[(res, model)]
                    use = {'unaware': trained['unaware'],
                           'unaware_insample': trained['unaware_insample']}
                    if model is not None:
                        use['aware'] = trained[f'aware:{model}']
                    for mname, mdl in use.items():
                        mus, sds = [], []
                        for k in range(F['dRrs'].shape[0]):
                            mu, sd = mdl.predict(F['dRrs'][k][te], F['aux'][k][te])
                            mus.append(mu), sds.append(sd)
                        mus, sds = np.array(mus), np.array(sds)          # (K, nte, q)
                        for j, t in enumerate(keys):
                            m = mus[:, :, j]
                            if kind == 'abs':
                                m = np.maximum(m, 0.0)
                            _store(preds, (name, f'sharedW_{mname}', t), te, m, sds[:, :, j])
                    for j, t in enumerate(keys):
                        coef, icpt = pcr[t]
                        m = np.array([F['dRrs'][k][te] @ D2.T @ coef + icpt
                                      for k in range(F['dRrs'].shape[0])])
                        if kind == 'abs':
                            m = np.maximum(m, 0.0)
                        _store(preds, (name, 'PCR_dRrs_d2', t), te, m, None)
                        # nulls: trained on clean retrieved Tchla, applied to the draw's Tchla
                        tg = bench.targets[t]
                        for nname, key in (('null_GSM', 'chl'), ('null_OC4', 'oc4')):
                            x_tr = np.log10(clean[key][0][tr])
                            b, a = np.polyfit(x_tr, tg['ylog'][tr], 1)
                            lx = np.log10(np.clip(F[key][:, te], 1e-3, None))
                            m = a + b * lx
                            if kind == 'abs':
                                m = 10.0**m
                            _store(preds, (name, nname, t), te, m, None)
        print(f'  LOCO fold {fold}: {time.time() - t0:.0f}s', flush=True)
    return preds


def _store(preds, key, te, mu, sd):
    n = 145
    if key not in preds:
        Kd = mu.shape[0]
        preds[key] = (np.full((Kd, n), np.nan), None if sd is None else np.full((Kd, n), np.nan))
    preds[key][0][:, te] = mu
    if sd is not None:
        preds[key][1][:, te] = sd


# --------------------------------------------------------------------- scores
def score_all(bench, preds):
    rows, pits = [], {}
    for (scen, mname, t), (mu, sd) in preds.items():
        tg = bench.targets[t]
        y = tg['y']
        Kd = mu.shape[0]
        Y = np.broadcast_to(y, mu.shape)
        if tg['kind'] == 'ratio':
            r2 = np.corrcoef(Y.ravel(), mu.ravel())[0, 1]**2
        else:
            r2 = V.score_linear(Y.ravel(), mu.ravel())['R2']
        row = {'scenario': scen, 'model': mname, 'target': t, 'kind': tg['kind'],
               'R2': float(r2), 'rmse': float(np.sqrt(np.mean((mu - Y)**2)))}
        if sd is not None:
            var_in = mu.var(axis=0, ddof=1) if Kd > 1 else np.zeros(mu.shape[1])
            tot = np.sqrt(sd**2 + var_in[None, :])
            g = bayes.gaussian_scores(Y, mu, tot)
            gm = bayes.gaussian_scores(Y, mu, sd)
            row.update({'cov68': g['cov68'], 'cov95': g['cov95'], 'z_sd': g['z_sd'],
                        'crps': g['crps'], 'cov68_model_only': gm['cov68'],
                        'cov95_model_only': gm['cov95'],
                        'sd_model_median': float(np.median(sd)),
                        'sd_input_median': float(np.median(np.sqrt(var_in))),
                        'sd_total_median': float(np.median(tot))})
            if tg['kind'] == 'ratio':
                pits.setdefault((scen, mname), []).append(g['pit'])
        rows.append(row)
    df = pd.DataFrame(rows)
    # skill above the best null (per scenario/target)
    nul = df[df.model.isin(['null_GSM', 'null_OC4'])].groupby(['scenario', 'target'])['R2'].max()
    df['R2_null_best'] = [nul.get((s, t), np.nan) for s, t in zip(df.scenario, df.target)]
    df['dR2_null'] = df['R2'] - df['R2_null_best']
    return df, {k: np.concatenate(v) for k, v in pits.items()}


def summarize(df):
    rows = []
    for (scen, m, kind), g in df.groupby(['scenario', 'model', 'kind'], sort=False):
        row = {'scenario': scen, 'model': m, 'kind': kind, 'mean_R2': g['R2'].mean(),
               'mean_dR2_null': g['dR2_null'].mean(),
               'n_dR2_gt_0.05': int((g['dR2_null'] > 0.05).sum())}
        for t in KEY_RATIOS:
            if t in set(g.target):
                row[f'R2_{t.split(":")[1]}'] = float(g[g.target == t]['R2'].iloc[0])
        if 'cov68' in g and g['cov68'].notna().any():
            row.update({'cov68': g['cov68'].mean(), 'cov95': g['cov95'].mean(),
                        'cov68_model_only': g['cov68_model_only'].mean(),
                        'cov95_model_only': g['cov95_model_only'].mean(),
                        'z_sd': g['z_sd'].mean(), 'crps': g['crps'].mean()})
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- figures
SCEN_ORDER = list(SCENARIOS)


def fig_coverage(summ, path):
    s = summ[(summ.kind == 'ratio') & summ.model.str.startswith('sharedW')]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, lev, nom in ((axes[0], 'cov68', 0.68), (axes[1], 'cov95', 0.95)):
        for k, m in enumerate(['sharedW_unaware', 'sharedW_aware']):
            g = s[s.model == m].set_index('scenario').reindex(SCEN_ORDER)
            x = np.arange(len(SCEN_ORDER)) + (k - 0.5) * 0.35
            ax.bar(x, g[lev], width=0.33, label=m.replace('sharedW_', 'shared W, '))
            ax.plot(x, g[f'{lev}_model_only'], 'k_', ms=12, mew=1.5,
                    label='model variance only' if k == 0 else None)
        ax.axhline(nom, color='r', ls='--', lw=1)
        ax.set_xticks(range(len(SCEN_ORDER)))
        ax.set_xticklabels(SCEN_ORDER, rotation=30, fontsize=8)
        ax.set_title(f'{int(nom * 100)}% interval coverage (12 log ratios, LOCO)', fontsize=9)
        ax.set_ylim(0, 1)
    axes[0].legend(fontsize=7.5, frameon=False, loc='lower left')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_pit(pits, path):
    fig, axes = plt.subplots(2, len(SCEN_ORDER), figsize=(14, 4.6), sharey=True)
    for i, m in enumerate(['sharedW_unaware', 'sharedW_aware']):
        for j, scen in enumerate(SCEN_ORDER):
            ax = axes[i, j]
            if (scen, m) in pits:
                ax.hist(pits[(scen, m)], bins=10, range=(0, 1), density=True, color=f'C{i}')
            ax.axhline(1, color='k', lw=0.6)
            ax.set_ylim(0, 3)
            if i == 0:
                ax.set_title(scen, fontsize=8)
            if j == 0:
                ax.set_ylabel(m.replace('sharedW_', ''), fontsize=8)
    fig.suptitle('PIT histograms, 12 log ratios pooled, LOCO (flat = calibrated)', fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_skill(df, path):
    fig, axes = plt.subplots(1, len(KEY_RATIOS) + 1, figsize=(15, 3.6), sharey=True)
    mods = [('sharedW_unaware', 'C0'), ('sharedW_aware', 'C1'), ('PCR_dRrs_d2', 'C3')]
    x = np.arange(len(SCEN_ORDER))
    for ax, t in zip(axes, KEY_RATIOS + ['mean']):
        for m, c in mods:
            if t == 'mean':
                g = df[(df.kind == 'ratio') & (df.model == m)].groupby('scenario')['R2'].mean()
                gn = df[df.kind == 'ratio'].groupby('scenario')['R2_null_best'].mean()
            else:
                g = df[(df.target == t) & (df.model == m)].set_index('scenario')['R2']
                gn = df[df.target == t].groupby('scenario')['R2_null_best'].first()
            ax.plot(x, g.reindex(SCEN_ORDER), 'o-', color=c, label=m, ms=4)
        ax.plot(x, gn.reindex(SCEN_ORDER), 'k--', lw=1, label='best Tchla null')
        ax.set_xticks(x)
        ax.set_xticklabels(SCEN_ORDER, rotation=60, fontsize=7)
        ax.set_title(t if t != 'mean' else 'mean over 12 ratios', fontsize=9)
        ax.grid(lw=0.3)
    axes[0].set_ylabel('LOCO log-R²')
    axes[0].legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_examples(bench, preds, data, path):
    t = 'ratio:Fuco'
    y = bench.targets[t]['y']
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    for ax, scen in zip(axes, ['clean_1nm', 'insitu_1nm', 'pace_white_1nm']):
        mu, sd = preds[(scen, 'sharedW_aware' if scen != 'clean_1nm' else 'sharedW_unaware', t)]
        var_in = mu.var(axis=0, ddof=1) if mu.shape[0] > 1 else 0 * mu[0]
        m0, tot = mu[0], np.sqrt(sd[0]**2 + var_in)
        ax.errorbar(y, m0, yerr=tot, fmt='o', ms=3, lw=0.6, alpha=0.7)
        lim = [y.min() - 0.3, y.max() + 0.3]
        ax.plot(lim, lim, 'k-', lw=0.6)
        ax.set_title(f'{scen}: one draw, ±1σ (total)', fontsize=9)
        ax.set_xlabel('log10 Fuco:Tchla, HPLC')
    axes[0].set_ylabel('predicted (LOCO)')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    z = np.load(sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz')
    bench = V.Benchmark(data.pigments13, z['params'][:, 0], V.oc4v6(data.wave, data.Rrs),
                        data.meta['campaign'].to_numpy(), n_perm=2, seed=SEED)
    rng = np.random.default_rng(SEED)
    t0 = time.time()
    feats, info = build_scenarios(data, rng)
    chk = float(np.max(np.abs(feats[(1, None)]['dRrs'][0] - z['dRrs'])))
    print(f'features done ({time.time() - t0:.0f}s); clean 1 nm δRrs vs #2 product: {chk:.2e}')
    preds = run_loco(data, bench, feats)
    df, pits = score_all(bench, preds)
    summ = summarize(df)
    pd.set_option('display.width', 250)
    cols = ['scenario', 'model', 'kind', 'mean_R2', 'mean_dR2_null', 'n_dR2_gt_0.05',
            'R2_Fuco', 'R2_Zea', 'R2_DVchla', 'R2_Chlc12', 'cov68', 'cov95',
            'cov68_model_only', 'cov95_model_only', 'z_sd', 'crps']
    print(summ[[c for c in cols if c in summ]].round(3).to_string())
    df.round(4).to_csv(OUT / 'uncertainty_scores.csv', index=False)
    summ.round(4).to_csv(OUT / 'uncertainty_summary.csv', index=False)
    fig_coverage(summ, OUT / 'uncertainty_coverage.png')
    fig_pit(pits, OUT / 'uncertainty_pit.png')
    fig_skill(df, OUT / 'uncertainty_skill.png')
    fig_examples(bench, preds, data, OUT / 'uncertainty_examples.png')
    out = {'script': 'scripts/sdp/uncertainty.py', 'epft_up_version': __version__,
           'input': data.provenance, 'K_draws': K, 'seed': SEED,
           'grid': {'lam_noise': LAM_NOISE, 'lam_smooth': LAM_SMOOTH, 'ranks': RANKS},
           'scenarios': {k: {'res_nm': v[0], 'noise': v[1]} for k, v in SCENARIOS.items()},
           'noise_settings': {'insitu_pct': noise.INSITU_PCT, 'pace_ell_nm': noise.PACE_ELL,
                              'pace_white_frac': noise.PACE_WHITE_FRAC},
           'feature_info': info, 'clean_1nm_check_max_abs': chk,
           'summary': json.loads(summ.to_json(orient='records')),
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(OUT / 'uncertainty_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f'wrote {OUT / "uncertainty_summary.json"}')


if __name__ == '__main__':
    main()
