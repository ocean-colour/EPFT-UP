"""Execution #5: what does the residual buy? Predictor (source-space) comparison.

Every predictor space is scored with the same Kramer PCR and the standard
:class:`epft_up.sdp.validate.Benchmark`: 13 absolute pigments plus 12 log10
pigment:Tchla ratios, random (100) and leave-one-campaign-out splits, each
against the Tchla-only null models (and a constant ratio).

Predictor spaces (N = 145):

==================  ==========================================================
``Rrs``             measured above-water Rrs (Lange-style PCR on raw spectra)
``Rrs_d2``          ``diff(Rrs, 2)``
``dRrs``            Kramer's residual δRrs = Rrs - Rrs_mod
``dRrs_d2``         ``diff(δRrs, 2)`` (the paper)
``M1``              El Hourany & Kramer (2026) spline residual of rrs
``M3_w{7,11,21}``   El Hourany M3: Savitzky-Golay (order 3) 2nd derivative of
                    rrs, windows 7, 11 and 21 nm (the paper gives no window)
``dRrs+GSM``        δRrs ⊕ {log10 Tchla_GSM, log10 a_dg(443), b_bp(443)}
``dRrs_d2+GSM``     δRrs'' ⊕ the same
``GSM3``            the three GSM parameters alone
``dRrs_d2_flat30``  δRrs'' with A,B Gaussian-smoothed (FWHM 30 nm) in the GSM fit
``dRrs_d2_flat80``  the same with FWHM 80 nm
``dRrs_flat80``     δRrs with FWHM-80 nm tables
``dRrs_awsm``       δRrs with a_w, A, B smoothed like the data (5 nm mean)
``dRrs_d2_awsm``    its ``diff(·, 2)``
==================  ==========================================================

Also: how much of the Rrs'' variance the GSM model's own curvature
(Rrs_mod'') accounts for, i.e. what δRrs'' removes from Rrs''.

Outputs (``reports/figures/sdp/``): ``source_spaces_table.csv`` (all
rows), ``source_spaces_summary.csv``, ``source_spaces_heatmap.png``,
``source_spaces_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/source_spaces.py
    conda run -n ocean14 python scripts/sdp/source_spaces.py --replot   # figure only
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
from epft_up.sdp import data as sdpdata, gsm, spectral, validate as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
SEED = 1


def build_spaces(data, z):
    """Predictor matrices (n_samples, n_features) for every source space."""
    w = data.wave
    T, S = data.meta['temp'].to_numpy(), data.meta['sal'].to_numpy()
    rrs = gsm.Rrs_to_rrs(data.Rrs)
    dRrs = z['dRrs']
    params = z['params']
    gsm3 = np.column_stack([np.log10(params[:, 0]), np.log10(params[:, 1]), params[:, 2]])
    d2 = lambda y: spectral.difference(w, y, order=2)[1]  # noqa: E731

    spaces = {
        'Rrs': data.Rrs,
        'Rrs_d2': d2(data.Rrs),
        'dRrs': dRrs,
        'dRrs_d2': d2(dRrs),
        'M1': spectral.spline_residual(w, rrs),
        'M3_w7': spectral.savgol_second_derivative(w, rrs, 7),
        'M3_w11': spectral.savgol_second_derivative(w, rrs, 11),
        'M3_w21': spectral.savgol_second_derivative(w, rrs, 21),
        'dRrs+GSM': np.hstack([dRrs, gsm3]),
        'dRrs_d2+GSM': np.hstack([d2(dRrs), gsm3]),
        'GSM3': gsm3,
    }
    tables = gsm.load_ref_tables(w)
    flat_info = {}
    for fw in (30, 80):
        fit = gsm.fit_gsm(w, data.Rrs, T, S, tables=gsm.smoothed_tables(tables, fw))
        flat_info[fw] = {'n_converged': int(fit.converged.sum()),
                         'Tchla_log_R2_vs_HPLC': float(np.corrcoef(
                             np.log10(data.pigments['Tchla']),
                             np.log10(np.where(fit.params[:, 0] > 0, fit.params[:, 0],
                                               np.nan)))[0, 1]**2)
                         if np.all(fit.params[:, 0] > 0) else None,
                         'dRrs_corr_with_original': float(np.corrcoef(
                             fit.dRrs.ravel(), dRrs.ravel())[0, 1]),
                         'dRrs_d2_corr_with_original': float(np.corrcoef(
                             d2(fit.dRrs).ravel(), d2(dRrs).ravel())[0, 1])}
        spaces[f'dRrs_d2_flat{fw}'] = d2(fit.dRrs)
        if fw == 80:
            spaces['dRrs_flat80'] = fit.dRrs
    # Consistent smoothing: the deposited Rrs carry a 5 nm moving mean, the a_w
    # table does not. Smooth a_w (and A, B) the same way before the GSM fit.
    sm = lambda x: np.r_[x[:2], spectral.moving_mean(x, 5)[2:-2], x[-2:]]  # noqa: E731
    t_sm = gsm.RefTables(w, sm(tables.A), sm(tables.B), sm(tables.aw),
                         dict(tables.source, smoothed='5 nm moving mean (A, B, aw)'))
    fit = gsm.fit_gsm(w, data.Rrs, T, S, tables=t_sm)
    flat_info['aw_smoothed'] = {'n_converged': int(fit.converged.sum()),
                                'dRrs_d2_corr_with_original': float(np.corrcoef(
                                    d2(fit.dRrs).ravel(), d2(dRrs).ravel())[0, 1])}
    spaces['dRrs_awsm'] = fit.dRrs
    spaces['dRrs_d2_awsm'] = d2(fit.dRrs)
    return spaces, flat_info, fit


def curvature_budget(data, z):
    """How much of Rrs'' (across-sample variance) the GSM curvature Rrs_mod'' carries."""
    w = data.wave
    a = spectral.difference(w, data.Rrs, 2)[1]
    m = spectral.difference(w, z['Rrs_mod'], 2)[1]
    d = a - m

    def var(x):
        return np.sum((x - x.mean(axis=0))**2)

    return {'frac_var_Rrs_d2_removed': float(1.0 - var(d) / var(a)),
            'frac_var_Rrs_d2_in_model': float(var(m) / var(a)),
            'mean_band_corr_Rrs_d2_vs_dRrs_d2': float(np.mean(
                [np.corrcoef(a[:, j], d[:, j])[0, 1] for j in range(a.shape[1])])),
            'note': 'variance across samples, summed over 299 bands'}


def summarize(df):
    """Per predictor: mean R² and counts of null-beating targets."""
    rows = []
    for model, g in df.groupby('model', sort=False):
        a, r = g[g.kind == 'abs'], g[g.kind == 'ratio']
        rows.append({
            'model': model,
            'abs_R2_random': a['R2_random'].mean(), 'abs_R2_loco': a['R2_loco'].mean(),
            'abs_dR2_null_random': a['dR2_random'].mean(),
            'abs_dR2_null_loco': a['dR2_loco'].mean(),
            'abs_beats_null_random': int(a['beats_null_random'].sum()),
            'abs_beats_null_loco': int(a['beats_null_loco'].sum()),
            'ratio_logR2_random': r['R2_random'].mean(), 'ratio_logR2_loco': r['R2_loco'].mean(),
            'ratio_dR2_null_loco': r['dR2_loco'].mean(),
            'ratio_beats_null_random': int(r['beats_null_random'].sum()),
            'ratio_beats_null_loco': int(r['beats_null_loco'].sum()),
            'ratio_beats_const_rms_loco': int((r['RMS_loco'] < r['RMS_const_loco']).sum()),
        })
    return pd.DataFrame(rows).set_index('model')


def heatmap(df, path):
    order = list(dict.fromkeys(df['model']))
    targets = list(dict.fromkeys(df.index))
    targets = ([t for t in targets if t.startswith('abs:')]
               + [t for t in targets if t.startswith('ratio:')])
    n_abs = sum(t.startswith('abs:') for t in targets)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), sharey=True)
    for ax, scheme in zip(axes, ('random', 'loco')):
        M = np.full((len(order), len(targets)), np.nan)
        for i, m in enumerate(order):
            g = df[df.model == m]
            M[i] = [g.loc[t, f'dR2_{scheme}'] for t in targets]
        im = ax.imshow(M, cmap='RdBu', vmin=-0.5, vmax=0.5, aspect='auto')
        ax.set_xticks(range(len(targets)))
        ax.set_xticklabels([t.replace('ratio:', 'r:').replace('abs:', '') for t in targets],
                           rotation=90, fontsize=7)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order, fontsize=8)
        ax.axvline(n_abs - 0.5, color='k', lw=0.8)
        ax.set_title(f'R²(model) − R²(best Tchla null), {scheme}', fontsize=9)
        for i, m in enumerate(order):
            g = df[df.model == m]
            for j, t in enumerate(targets):
                if g.loc[t, f'beats_null_{scheme}']:
                    ax.text(j, i, '•', ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=axes, shrink=0.7, label='ΔR²')
    fig.suptitle('Source spaces under Kramer PCR: absolute pigments (left of line) and '
                 'log ratios; • = beats the null (random ≥ 90% of splits / LOCO CI > 0)',
                 fontsize=9)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if '--replot' in sys.argv:
        heatmap(pd.read_csv(OUT / 'source_spaces_table.csv', index_col=0),
                OUT / 'source_spaces_heatmap.png')
        return
    data = sdpdata.load_kramer2022()
    gsm_prod = sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz'
    z = np.load(gsm_prod)
    bench = V.Benchmark(data.pigments13, z['params'][:, 0], V.oc4v6(data.wave, data.Rrs),
                        data.meta['campaign'].to_numpy(), n_perm=100, seed=SEED)
    spaces, flat_info, fit_awsm = build_spaces(data, z)
    budget = curvature_budget(data, z)
    budget_awsm = curvature_budget(data, {'Rrs_mod': fit_awsm.Rrs_mod})
    print('curvature budget (aw smoothed):', budget_awsm)
    print('curvature budget:', budget)
    print('flat-baseline GSM fits:', flat_info)

    frames = []
    for name, X in spaces.items():
        t0 = time.time()
        frames.append(bench.evaluate(bench.pcr_factory(X), name))
        print(f'{name:16s} p={X.shape[1]:3d}  {time.time() - t0:5.1f}s')
    df = pd.concat(frames)
    summ = summarize(df)
    pd.set_option('display.width', 250)
    print(summ.round(3).to_string())

    df.round(4).to_csv(OUT / 'source_spaces_table.csv')
    summ.round(4).to_csv(OUT / 'source_spaces_summary.csv')
    heatmap(df, OUT / 'source_spaces_heatmap.png')
    out = {'script': 'scripts/sdp/source_spaces.py', 'epft_up_version': __version__,
           'input': data.provenance,
           'gsm_product': {'file': str(gsm_prod), 'sha256': sdpdata.sha256_of(gsm_prod)},
           'benchmark': {'n_perm': 100, 'seed': SEED, 'zero_frac': bench.zero_frac,
                         'lods': bench.lods, 'loco': bench.loco_names},
           'n_features': {k: int(v.shape[1]) for k, v in spaces.items()},
           'm1_alpha': spectral.EH_ALPHA, 'curvature_budget': budget,
           'curvature_budget_aw_smoothed': budget_awsm,
           'flat_baselines': flat_info, 'summary': json.loads(summ.to_json(orient='index')),
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(OUT / 'source_spaces_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f'wrote {OUT / "source_spaces_summary.json"}')


if __name__ == '__main__':
    main()
