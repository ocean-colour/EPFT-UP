"""Execution #7: learned wavelength weighting vs Kramer's PCR, judged by the Benchmark.

All models are scored on the standard Benchmark (13 absolute pigments plus 12
log ratios; random and leave-one-campaign-out splits; each against the
Tchla-only null models). Log ratios are the primary target.

==========================  =====================================================
``PCR_dRrs_d2``             Kramer 2022 (reference)
``ridge_dRrs_d2``           ridge (GCV) on z-scored δRrs'' (estimator change only)
``PLS_dRrs_d2``             PLS (inner CV) on z-scored δRrs''
``ridge_dRrs_aux``          ridge (GCV) on z-scored δRrs ⊕ GSM parameters (basis change)
``PLS_dRrs_aux``            PLS on the same
``smooth_dRrs_aux_group``   penalized (noise + smoothness), per pigment (rank 1),
                            hyperparameters by inner leave-campaign-out
``RRR_dRrs_aux_random``     **shared low-rank W**, penalized, inner random 5-fold
``RRR_dRrs_aux_group``      shared low-rank W, inner leave-campaign-out
``RRR_M1_aux_group``        the same on El Hourany M1 instead of δRrs
==========================  =====================================================

The GSM parameters are log10 Tchla_GSM, log10 a_dg(443) and b_bp(443). The
noise covariance Σ_n is the in-situ-like ``pct:0.02`` + floor model through
the 5 nm mean (report §5.2), at the median spectrum.

A full-data fit of the shared W (hyperparameters by leave-campaign-out CV)
supplies the weight-spectrum figures and the wavelength-importance summary.

Outputs (``reports/figures/sdp/``): ``weighting_table.csv``,
``weighting_summary.csv``, ``weighting_ranks.csv``,
``weighting_latent_spectra.png``, ``weighting_target_spectra.png``,
``weighting_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/learned_weighting.py
"""
import datetime
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))
sys.path.insert(0, os.path.dirname(__file__))

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import data as sdpdata, noise, spectral, theory, gsm  # noqa: E402
from epft_up.sdp import validate as V, weighting as W  # noqa: E402
from source_spaces import summarize  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
SEED = 1
LAM_NOISE = (0.0, 0.1, 1.0, 10.0, 100.0)
LAM_SMOOTH = (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0)
RANKS = tuple(range(1, 11))
BANDS = [(400, 450), (450, 500), (500, 550), (550, 600), (600, 650), (650, 700)]


def noise_cov(data):
    M = theory.boxcar_matrix(len(data.wave), 5)
    s = noise.pct_floor_sigma(np.median(data.Rrs, axis=0), 0.02) * np.sqrt(5.0)
    return M @ np.diag(s**2) @ M.T


def per_pigment_smooth_factory(Xs, aux, Sn, campaigns):
    def factory(y, constraint):
        model = W.PenalizedRRR(Xs, aux, Sn, LAM_NOISE, LAM_SMOOTH, ranks=(1,))
        return W.SharedRRRFitter(model, y[:, None], campaigns, inner='group',
                                 seed=SEED)(y, constraint)
    return factory


class KindRouter:
    """Route Benchmark targets to one SharedRRRFitter per target kind."""

    def __init__(self, bench, make_model, campaigns, inner):
        self.fitters = {}
        for kind in ('abs', 'ratio'):
            keys = [t for t, v in bench.targets.items() if v['kind'] == kind]
            Y = np.column_stack([bench.targets[t]['y'] for t in keys])
            self.fitters[kind] = W.SharedRRRFitter(make_model(), Y, campaigns,
                                                   inner=inner, seed=SEED)

    def __call__(self, y, constraint):
        kind = 'abs' if constraint == 'pigment' else 'ratio'
        return self.fitters[kind](y, constraint)


def full_fit(bench, Xs, aux, Sn, campaigns, kind):
    """Shared W on all data, hyperparameters by leave-campaign-out CV."""
    keys = [t for t, v in bench.targets.items() if v['kind'] == kind]
    Y = np.column_stack([bench.targets[t]['y'] for t in keys])
    model = W.PenalizedRRR(Xs, aux, Sn, LAM_NOISE, LAM_SMOOTH, RANKS)
    idx = np.arange(Y.shape[0])
    folds = [idx[campaigns == c] for c in dict.fromkeys(campaigns)]
    lam_n, lam_s, r = model.select(idx, Y, folds)
    B, st, Ym, Ysd = model.fit_full(idx, Y, lam_n, lam_s)
    X, _ = model._design(idx, st)
    Br, Vr = W.PenalizedRRR.reduce_rank(X, B, r)
    p = Xs.shape[1]
    latent = (B @ Vr)[:p]                     # (p, r) latent weight spectra
    weights = Br[:p] * Ysd                    # raw-unit weight spectrum per target
    return {'keys': keys, 'hp': (lam_n, lam_s, r), 'latent': latent,
            'weights': weights, 'aux_weights': Br[p:] * Ysd, 'Vr': Vr}


def band_importance(wave, weights, Xs):
    """Share of |w(λ)| × SD(Xs(λ)) in each wavelength band, per target."""
    contrib = np.abs(weights) * Xs.std(axis=0, ddof=1)[:, None]
    tot = contrib.sum(axis=0)
    return {f'{a}-{b}': (contrib[(wave >= a) & (wave < b)].sum(axis=0) / tot)
            for a, b in BANDS}


def fig_latent(wave, fit, path):
    r = fit['latent'].shape[1]
    fig, axes = plt.subplots(r, 1, figsize=(8, 2.2 * r + 0.5), sharex=True, squeeze=False)
    for i in range(r):
        ax = axes[i, 0]
        ax.plot(wave, fit['latent'][:, i], 'k-', lw=1)
        ax.axhline(0, color='0.6', lw=0.5)
        load = fit['Vr'][:, i]
        top = np.argsort(-np.abs(load))[:4]
        ax.set_title(f'latent weight spectrum {i + 1}; target loadings: '
                     + ', '.join(f'{fit["keys"][j].split(":")[1]} {load[j]:+.2f}' for j in top),
                     fontsize=8)
    axes[-1, 0].set_xlabel('Wavelength [nm]')
    fig.suptitle(f'Shared W on δRrs ⊕ GSM, log ratios (λ_n={fit["hp"][0]:g}, '
                 f'λ_s={fit["hp"][1]:g}, rank {fit["hp"][2]})', fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_targets(wave, fit_ratio, fit_abs, pcr_w, path):
    show = [('ratio:Fuco', fit_ratio), ('ratio:DVchla', fit_ratio), ('ratio:Zea', fit_ratio),
            ('ratio:HexFuco', fit_ratio), ('abs:Tchla', fit_abs), ('abs:Fuco', fit_abs)]
    fig, axes = plt.subplots(3, 2, figsize=(11, 8), sharex=True)
    for ax, (t, fit) in zip(axes.ravel(), show):
        j = fit['keys'].index(t)
        w = fit['weights'][:, j]
        ax.plot(wave, w / np.max(np.abs(w)), 'C0-', lw=1.2, label='shared W (this work)')
        if t in pcr_w:
            v = pcr_w[t]
            ax.plot(wave, v / np.max(np.abs(v)), 'C3-', lw=0.5, alpha=0.6,
                    label='PCR on δRrs″, effective (Dᵀ A)')
        ax.axhline(0, color='0.6', lw=0.5)
        ax.set_title(t, fontsize=9)
        ax.set_ylabel('weight on δRrs (normalized)', fontsize=7)
    axes[0, 0].legend(fontsize=7, frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel('Wavelength [nm]')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    gsm_prod = sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz'
    z = np.load(gsm_prod)
    camp = data.meta['campaign'].to_numpy()
    bench = V.Benchmark(data.pigments13, z['params'][:, 0], V.oc4v6(data.wave, data.Rrs),
                        camp, n_perm=100, seed=SEED)
    dRrs = z['dRrs']
    d2 = spectral.difference(data.wave, dRrs, 2)[1]
    p = z['params']
    aux = np.column_stack([np.log10(p[:, 0]), np.log10(p[:, 1]), p[:, 2]])
    Xaux = np.hstack([dRrs, aux])
    Sn = noise_cov(data)
    m1 = spectral.spline_residual(data.wave, gsm.Rrs_to_rrs(data.Rrs))

    rrr = lambda Xs: (lambda: W.PenalizedRRR(Xs, aux, Sn, LAM_NOISE, LAM_SMOOTH, RANKS))  # noqa: E731
    routers = {'RRR_dRrs_aux_random': KindRouter(bench, rrr(dRrs), camp, 'random'),
               'RRR_dRrs_aux_group': KindRouter(bench, rrr(dRrs), camp, 'group'),
               'RRR_M1_aux_group': KindRouter(bench, rrr(m1), camp, 'group')}
    models = {
        'PCR_dRrs_d2': bench.pcr_factory(d2),
        'ridge_dRrs_d2': lambda y, c: W.RidgeGCVFitter(d2, y, c),
        'PLS_dRrs_d2': lambda y, c: W.PLSFitter(d2, y, c, seed=SEED),
        'ridge_dRrs_aux': lambda y, c: W.RidgeGCVFitter(Xaux, y, c),
        'PLS_dRrs_aux': lambda y, c: W.PLSFitter(Xaux, y, c, seed=SEED),
        'smooth_dRrs_aux_group': per_pigment_smooth_factory(dRrs, aux, Sn, camp),
        **routers}

    frames, timing = [], {}
    for name, fac in models.items():
        t0 = time.time()
        frames.append(bench.evaluate(fac, name))
        timing[name] = time.time() - t0
        print(f'{name:24s} {timing[name]:6.1f}s', flush=True)
    df = pd.concat(frames)
    summ = summarize(df)
    pd.set_option('display.width', 250)
    print(summ.round(3).to_string())

    # selected hyperparameters (ranks) of the shared models, per scheme mix
    ranks = {}
    for name, rt in routers.items():
        for kind, f in rt.fitters.items():
            c = Counter(h[2] for h in f.chosen)
            ranks[f'{name}:{kind}'] = {'rank_counts': dict(sorted(c.items())),
                                       'median_rank': float(np.median([h[2] for h in f.chosen])),
                                       'lam_noise_counts': dict(Counter(h[0] for h in f.chosen)),
                                       'lam_smooth_counts': dict(Counter(h[1] for h in f.chosen))}
    pd.DataFrame({k: {'median_rank': v['median_rank'], **{f'r{r}': v['rank_counts'].get(r, 0)
                                                           for r in RANKS}}
                  for k, v in ranks.items()}).T.to_csv(OUT / 'weighting_ranks.csv')
    print(json.dumps(ranks, indent=0, default=str)[:2500])

    # full-data shared W for the figures
    fit_ratio = full_fit(bench, dRrs, aux, Sn, camp, 'ratio')
    fit_abs = full_fit(bench, dRrs, aux, Sn, camp, 'abs')
    from epft_up.sdp import models as M
    pcr_w = {}
    D = theory.difference_matrix(len(data.wave), 2)
    for t in ('ratio:Fuco', 'ratio:DVchla', 'ratio:Zea', 'ratio:HexFuco', 'abs:Tchla', 'abs:Fuco'):
        y = bench.targets[t]['y']
        r = M.train_pcr_kramer(d2, y, constraint=None if t.startswith('ratio') else 'pigment',
                               seed=SEED)
        pcr_w[t] = np.median(r.coefs, axis=0) @ D
    fig_latent(data.wave, fit_ratio, OUT / 'weighting_latent_spectra.png')
    fig_targets(data.wave, fit_ratio, fit_abs, pcr_w, OUT / 'weighting_target_spectra.png')
    imp = {kind: {k: dict(zip(fit['keys'], np.round(v, 3).tolist()))
                  for k, v in band_importance(data.wave, fit['weights'], dRrs).items()}
           for kind, fit in (('ratio', fit_ratio), ('abs', fit_abs))}
    pcr_imp = {t: {k: float(v[0]) for k, v in band_importance(data.wave, w[:, None], dRrs).items()}
               for t, w in pcr_w.items()}

    df.round(4).to_csv(OUT / 'weighting_table.csv')
    summ.round(4).to_csv(OUT / 'weighting_summary.csv')
    out = {'script': 'scripts/sdp/learned_weighting.py', 'epft_up_version': __version__,
           'input': data.provenance,
           'gsm_product': {'file': str(gsm_prod), 'sha256': sdpdata.sha256_of(gsm_prod)},
           'grid': {'lam_noise': LAM_NOISE, 'lam_smooth': LAM_SMOOTH, 'ranks': RANKS},
           'noise_model': "pct:0.02 + floor through 5 nm mean, median spectrum",
           'timing_s': timing, 'summary': json.loads(summ.to_json(orient='index')),
           'selected': ranks,
           'full_fit': {'ratio_hp': fit_ratio['hp'], 'abs_hp': fit_abs['hp'],
                        'ratio_target_loadings': dict(zip(fit_ratio['keys'],
                                                          np.round(fit_ratio['Vr'], 3).tolist())),
                        'aux_weights_ratio': dict(zip(fit_ratio['keys'],
                                                      np.round(fit_ratio['aux_weights'].T, 4).tolist()))},
           'band_importance_sharedW': imp, 'band_importance_pcr_effective': pcr_imp,
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(OUT / 'weighting_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f'wrote {OUT / "weighting_summary.json"}')


if __name__ == '__main__':
    main()
