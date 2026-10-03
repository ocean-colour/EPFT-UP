"""Execution #6: numbers and figures for the Maths/Statistics section (report §5).

(i)   Linear-operator argument. PCR on δRrs'' = D δRrs has effective weights
      w = Dᵀ A on δRrs. Verify the identity, compare with PCR trained on
      δRrs itself, and characterize the prior implied by
      derivative + z-score + shrinkage: w ~ N(0, Dᵀ S⁻² D).
(ii)  Noise propagation and information content. Noise models: the
      deposit's rounding, IOPtics-style 'pct:0.02' with floor, and PACE OCI
      (ocpy), the last both white and spectrally correlated (ℓ = 30 nm, plus
      10% white). Instrument noise is passed through the paper's 5 nm mean.
      Per-band SNR of δRrs vs δRrs''; Rodgers' degrees of freedom for signal
      d_s at 1, 2.5, 5 and 10 nm sampling; the noise in PCR predictions.
(iii) PCR vs ridge vs PLS as spectral filters on z-scored δRrs'' (Tchla,
      Fuco).

Outputs (``reports/figures/sdp/``): ``maths_effective_weights.png``,
``maths_prior.png``, ``maths_snr.png``, ``maths_dof.png``,
``maths_filter_factors.png``, ``maths_summary.json``.

Run with::

    conda run -n ocean14 python scripts/sdp/maths_section.py
"""
import datetime
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
from sklearn.cross_decomposition import PLSRegression  # noqa: E402
from sklearn.model_selection import KFold  # noqa: E402

from epft_up import __version__  # noqa: E402
from epft_up.sdp import data as sdpdata, models, noise, theory  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'reports' / 'figures' / 'sdp'
PIGS = ('Tchla', 'Fuco', 'HexFuco')
SEED = 1
ELL = 30.0           # correlation length of the "smooth" PACE error model [nm]
PCT = 0.02           # in-situ-like fractional error with floor
STEPS = (1.0, 2.5, 5.0, 10.0)


# ---------------------------------------------------------------- part (i)
def part_i(wave, dRrs, pig):
    n = len(wave)
    D = theory.difference_matrix(n, 2)
    Xd2 = dRrs @ D.T
    out, curves = {}, {}
    for p in PIGS:
        y = pig[p].to_numpy()
        r2 = models.train_pcr_kramer(Xd2, y, seed=SEED)
        r0 = models.train_pcr_kramer(dRrs, y, seed=SEED)
        A = np.median(r2.coefs, axis=0)
        w_eff = theory.effective_weights(r2.coefs, D)          # (100, 301)
        ident = float(np.max(np.abs(Xd2 @ r2.coefs.T - dRrs @ w_eff.T)))
        w_med, w_dir = np.median(w_eff, axis=0), np.median(r0.coefs, axis=0)
        f, fa_eff = theory.frequency_power(w_med)
        _, fa_dir = theory.frequency_power(w_dir)
        i01 = int(np.searchsorted(f, 0.1))
        out[p] = {'identity_max_abs_diff': ident,
                  'sum_w_eff': float(np.median(w_eff.sum(axis=1))),
                  'first_moment_w_eff': float(np.median(w_eff @ (wave - wave.mean()))),
                  'frac_power_f_ge_0.1_eff': float(fa_eff[i01]),
                  'frac_power_f_ge_0.1_direct': float(fa_dir[i01]),
                  'corr_w_eff_w_direct': float(np.corrcoef(w_med, w_dir)[0, 1]),
                  'R2_random_d2': r2.summary()['R2'][0],
                  'R2_random_direct': r0.summary()['R2'][0]}
        curves[p] = {'A': A, 'w_eff': w_med, 'w_dir': w_dir, 'f': f,
                     'fa_eff': fa_eff, 'fa_dir': fa_dir}

    # implied priors
    s_d2 = Xd2.std(axis=0, ddof=1)
    s_0 = dRrs.std(axis=0, ddof=1)
    C_d2 = theory.implied_prior_covariance(D, s_d2)
    C_0 = theory.implied_prior_covariance(np.eye(n), s_0)
    rng = np.random.default_rng(0)
    draws_d2 = rng.standard_normal((2000, n - 2)) / s_d2 @ D
    draws_0 = rng.standard_normal((2000, n)) / s_0
    win = np.hanning(n)
    ps = lambda W: np.mean(np.abs(np.fft.rfft(W * win, axis=1))**2, axis=0)  # noqa: E731
    fr = np.fft.rfftfreq(n)
    P_d2, P_0 = ps(draws_d2), ps(draws_0)
    prior = {'null_space_check_max': float(np.max(np.abs(np.ones(n) @ C_d2))),
             'diag_var_red_over_blue_d2': float(np.mean(np.diag(C_d2)[260:])
                                                / np.mean(np.diag(C_d2)[10:60])),
             'diag_var_red_over_blue_direct': float(np.mean(np.diag(C_0)[260:])
                                                    / np.mean(np.diag(C_0)[10:60])),
             'prior_power_ratio_hi_lo_d2': float(P_d2[fr >= 0.3].mean() / P_d2[(fr > 0.005) & (fr < 0.05)].mean()),
             'prior_power_ratio_hi_lo_direct': float(P_0[fr >= 0.3].mean() / P_0[(fr > 0.005) & (fr < 0.05)].mean())}
    return out, curves, prior, (fr, P_d2, P_0, np.diag(C_d2), np.diag(C_0))


def fig_effective(wave, curves, path):
    fig, axes = plt.subplots(len(PIGS), 2, figsize=(11, 8.5))
    for i, p in enumerate(PIGS):
        c = curves[p]
        ax = axes[i, 0]
        ax.plot(wave, c['w_eff'], 'C3-', lw=0.8, label=r'$w_{\rm eff}=D^{\top}A$ (PCR on δRrs″)')
        ax.plot(wave, c['w_dir'], 'C0-', lw=1.2, label='w (PCR on δRrs directly)')
        ax.axhline(0, color='k', lw=0.4)
        ax.set_ylabel(r'weight on δRrs [mg m$^{-3}$ sr]', fontsize=8)
        ax.set_title(f'{p}: effective weight spectra on δRrs', fontsize=9)
        if i == 0:
            ax.legend(fontsize=7.5, frameon=False)
        ax = axes[i, 1]
        ax.semilogx(c['f'][1:], c['fa_eff'][1:], 'C3-', label='PCR on δRrs″')
        ax.semilogx(c['f'][1:], c['fa_dir'][1:], 'C0-', label='PCR on δRrs')
        ax.set_ylabel('fraction of weight power ≥ f', fontsize=8)
        ax.set_ylim(0, 1.02)
        ax.grid(lw=0.3)
        if i == 0:
            ax.legend(fontsize=7.5, frameon=False)
    axes[-1, 0].set_xlabel('Wavelength [nm]')
    axes[-1, 1].set_xlabel('frequency f [cycles nm⁻¹]')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_prior(wave, fr, P_d2, P_0, dg_d2, dg_0, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    ax = axes[0]
    ax.semilogy(wave, dg_d2 / dg_d2.mean(), 'C3-', label='derivative + z-score: diag(Dᵀ S⁻² D)')
    ax.semilogy(wave, dg_0 / dg_0.mean(), 'C0-', label='z-score only: diag(S⁻²)')
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('prior variance of w(λ) (normalized)')
    ax.legend(fontsize=7.5, frameon=False)
    ax = axes[1]
    ax.loglog(fr[1:], P_d2[1:] / P_d2[1:].mean(), 'C3-', label='derivative + z-score')
    ax.loglog(fr[1:], P_0[1:] / P_0[1:].mean(), 'C0-', label='z-score only')
    f = fr[1:]
    ax.loglog(f, (2 * np.sin(np.pi * f))**4 / np.mean((2 * np.sin(np.pi * f))**4), 'k:',
              lw=0.8, label=r'$(2\sin\pi f)^4$')
    ax.set_xlabel('frequency f [cycles nm⁻¹]')
    ax.set_ylabel('prior power of w (normalized)')
    ax.legend(fontsize=7.5, frameon=False)
    fig.suptitle('Implied isotropic-shrinkage prior on the effective weights', fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


# --------------------------------------------------------------- part (ii)
def interp_matrix(wave, step):
    """Linear interpolation from the 1 nm grid to a grid of spacing ``step``."""
    g = np.arange(wave[0], wave[-1] + 1e-9, step)
    M = np.zeros((len(g), len(wave)))
    for i, x in enumerate(g):
        j = min(int(np.floor(x - wave[0])), len(wave) - 2)
        t = x - wave[j]
        M[i, j], M[i, j + 1] = 1 - t, t
    return g, M


def noise_models(wave, Rrs_median):
    """Noise covariances [sr^-2] on the 1 nm grid, as seen by the predictors.

    Instrument-level noise (pct, PACE) is defined per 1 nm band and passed
    through the paper's 5 nm moving mean M (applied to every spectrum before
    δRrs, by Kramer and by the port): C = M diag(σ²) Mᵀ. The deposit's
    rounding happens after smoothing and stays white. The PACE σ (per band of
    ocpy's PACE_error.csv, ≈2 nm spacing) is rescaled to 1 nm bands (ocpy
    ``include_sampling``: × √2). The
    "correlated" PACE model is a smooth Gaussian-correlated error (ℓ = ELL,
    like an imperfect atmospheric correction) plus a white part of 10% of σ.
    Any real error has some white part, and a purely smooth covariance is
    numerically singular under differencing.
    """
    n = len(wave)
    M = theory.boxcar_matrix(n, 5)
    pace1 = noise.pace_sigma(wave, include_sampling=True)
    pace_native = noise.pace_sigma(wave)
    pct = noise.pct_floor_sigma(Rrs_median, PCT)
    rnd = noise.rounding_sigma(n)
    smooth_white = lambda s: M @ theory.white_covariance(s) @ M.T  # noqa: E731
    return {
        'rounding (deposit)': theory.white_covariance(rnd),
        f'pct:{PCT} + floor (in-situ-like)': smooth_white(pct * np.sqrt(5.0)),
        'PACE white': smooth_white(pace1),
        f'PACE correlated (ℓ={ELL:g} nm) + 10% white':
            theory.correlated_covariance(wave, pace_native, ELL) + smooth_white(0.1 * pace1),
    }


def part_ii(wave, dRrs, w_eff_by_pig, w_dir_by_pig, pig, Rrs):
    n = len(wave)
    nm = noise_models(wave, np.median(Rrs, axis=0))
    D = theory.difference_matrix(n, 2)
    snr = {}
    sig0 = dRrs.std(axis=0, ddof=1)
    sig2 = (dRrs @ D.T).std(axis=0, ddof=1)
    for name, C in nm.items():
        snr[name] = {'dRrs': sig0 / np.sqrt(np.diag(C)),
                     'dRrs_d2': sig2 / np.sqrt(np.diag(D @ C @ D.T))}
    dof = {}
    Cs_full = np.cov(dRrs, rowvar=False)
    for step in STEPS:
        g, I = (wave, np.eye(n)) if step == 1.0 else interp_matrix(wave, step)
        Cs0 = I @ Cs_full @ I.T
        Dh = theory.difference_matrix(len(g), 2)
        Cs2 = Dh @ Cs0 @ Dh.T
        for name, C in nm.items():
            Cn0 = I @ C @ I.T
            d0, l0 = theory.dof_signal(Cs0, Cn0)
            d2, l2 = theory.dof_signal(Cs2, Dh @ Cn0 @ Dh.T)
            dof.setdefault(name, {})[step] = {
                'n_bands': int(len(g)), 'd_s_dRrs': d0, 'd_s_dRrs_d2': d2,
                'n_lambda_gt1_dRrs': int(np.sum(l0 > 1)),
                'n_lambda_gt1_dRrs_d2': int(np.sum(l2 > 1))}
    pred = {}
    for p in PIGS:
        sd_p = float(pig[p].std())
        pred[p] = {'pigment_sd': sd_p}
        for name, C in nm.items():
            e = float(np.sqrt(w_eff_by_pig[p] @ C @ w_eff_by_pig[p]))
            d = float(np.sqrt(w_dir_by_pig[p] @ C @ w_dir_by_pig[p]))
            pred[p][name] = {'noise_sd_pcr_d2': e, 'noise_sd_pcr_direct': d,
                             'ratio_to_pigment_sd_d2': e / sd_p,
                             'ratio_to_pigment_sd_direct': d / sd_p}
    return snr, dof, pred


def fig_snr(wave, snr, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, key, title in ((axes[0], 'dRrs', 'δRrs'), (axes[1], 'dRrs_d2', 'δRrs″ (1 nm)')):
        for name, v in snr.items():
            x = wave if key == 'dRrs' else wave[1:-1]
            ax.semilogy(x, v[key], lw=1, label=name)
        ax.axhline(1, color='k', lw=0.6, ls='--')
        ax.set_xlabel('Wavelength [nm]')
        ax.set_title(f'per-band SNR of {title}: across-sample SD / noise SD', fontsize=9)
    axes[0].set_ylabel('SNR')
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def fig_dof(dof, path):
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    for k, (name, v) in enumerate(dof.items()):
        st = sorted(v)
        ax.plot(st, [v[s]['d_s_dRrs'] for s in st], 'o-', color=f'C{k}', label=f'{name}: δRrs')
        ax.plot(st, [v[s]['d_s_dRrs_d2'] for s in st], 's--', color=f'C{k}', mfc='none',
                label=f'{name}: δRrs″')
    ax.set_xscale('log')
    ax.set_xticks(STEPS)
    ax.set_xticklabels([f'{s:g}' for s in STEPS])
    ax.set_xlabel('sampling [nm]')
    ax.set_ylabel(r'degrees of freedom for signal $d_s$')
    ax.legend(fontsize=6.5, frameon=False)
    ax.set_title('Information content of δRrs and δRrs″ (Rodgers 2000)', fontsize=9)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


# -------------------------------------------------------------- part (iii)
def part_iii(dRrs, pig, k_pcr):
    D = theory.difference_matrix(dRrs.shape[1], 2)
    X = dRrs @ D.T
    Z = (X - X.mean(0)) / X.std(0, ddof=1)
    out, curves = {}, {}
    kf = KFold(5, shuffle=True, random_state=SEED)
    for p in ('Tchla', 'Fuco'):
        y = pig[p].to_numpy()
        lams = np.logspace(-1, 5, 40)
        cv_r = []
        for lam in lams:
            err = 0.0
            for tr, te in kf.split(Z):
                Zt = Z[tr] - Z[tr].mean(0)
                b = theory.ridge_coefficients(Zt, y[tr], lam)
                pr = (Z[te] - Z[tr].mean(0)) @ b + y[tr].mean()
                err += np.sum((np.maximum(pr, 0) - y[te])**2)
            cv_r.append(err)
        lam_best = float(lams[int(np.argmin(cv_r))])
        cv_p = []
        for nc in range(1, 21):
            err = 0.0
            for tr, te in kf.split(Z):
                m = PLSRegression(n_components=nc, scale=False).fit(Z[tr], y[tr])
                err += np.sum((np.maximum(m.predict(Z[te]).ravel(), 0) - y[te])**2)
            cv_p.append(err)
        nc_best = int(np.argmin(cv_p)) + 1
        Zc = Z - Z.mean(0)
        b_pcr = theory.pcr_coefficients(Zc, y, k_pcr[p])
        b_ridge = theory.ridge_coefficients(Zc, y, lam_best)
        b_pls = PLSRegression(n_components=nc_best, scale=False).fit(Zc, y).coef_.ravel()
        f_pcr, s = theory.filter_factors(Zc, y, b_pcr)
        f_ridge, _ = theory.filter_factors(Zc, y, b_ridge)
        f_pls, _ = theory.filter_factors(Zc, y, b_pls)
        out[p] = {'k_pcr': int(k_pcr[p]), 'ridge_lambda_cv': lam_best, 'pls_ncomp_cv': nc_best,
                  'eff_dof_pcr': float(np.sum(f_pcr)), 'eff_dof_ridge': float(np.sum(f_ridge)),
                  'eff_dof_pls': float(np.sum(f_pls)),
                  'pls_max_filter_factor': float(np.nanmax(np.abs(f_pls))),
                  'corr_beta_pcr_ridge': float(np.corrcoef(b_pcr, b_ridge)[0, 1]),
                  'corr_beta_pcr_pls': float(np.corrcoef(b_pcr, b_pls)[0, 1]),
                  'frac_var_top_k': float(np.sum(s[:k_pcr[p]]**2) / np.sum(s**2))}
        curves[p] = (s, f_pcr, f_ridge, f_pls)
    return out, curves


def fig_filters(curves, info, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, (p, (s, fp, fr, fl)) in zip(axes, curves.items()):
        i = np.arange(1, len(s) + 1)
        ax.plot(i, fp, 'k-', drawstyle='steps-mid', label=f'PCR (k={info[p]["k_pcr"]})')
        ax.plot(i, fr, 'C0-', label=f'ridge (λ={info[p]["ridge_lambda_cv"]:.3g}, CV)')
        ax.plot(i, fl, 'C1.', ms=4, label=f'PLS ({info[p]["pls_ncomp_cv"]} comps, CV)')
        ax.axhline(1, color='0.6', lw=0.5)
        ax.set_xlim(0, 80)
        ax.set_ylim(-1.5, 2.5)
        ax.set_xlabel('singular-value index i (z-scored δRrs″)')
        ax.set_title(f'{p}: filter factors $f_i$', fontsize=9)
        ax2 = ax.twinx()
        ax2.semilogy(i, s**2 / np.sum(s**2), color='0.7', lw=0.8)
        ax2.set_ylabel('variance fraction $s_i^2/\\sum s^2$', fontsize=7, color='0.5')
    axes[0].set_ylabel('$f_i$')
    axes[0].legend(fontsize=7.5, frameon=False, loc='upper right')
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = sdpdata.load_kramer2022()
    gsm_prod = sdpdata.kramer2022_dir() / 'products' / 'gsm_dRrs_insitu.npz'
    z = np.load(gsm_prod)
    wave, dRrs, pig = data.wave, z['dRrs'], data.pigments13

    res_i, curves_i, prior, (fr, P_d2, P_0, dg_d2, dg_0) = part_i(wave, dRrs, pig)
    print('(i)', json.dumps(res_i, indent=0)[:1500])
    print('(i) prior', prior)
    fig_effective(wave, curves_i, OUT / 'maths_effective_weights.png')
    fig_prior(wave, fr, P_d2, P_0, dg_d2, dg_0, OUT / 'maths_prior.png')

    snr, dof, pred = part_ii(wave, dRrs, {p: curves_i[p]['w_eff'] for p in PIGS},
                             {p: curves_i[p]['w_dir'] for p in PIGS}, pig, data.Rrs)
    for name, v in dof.items():
        print('(ii) d_s', name, {s: (round(x['d_s_dRrs'], 1), round(x['d_s_dRrs_d2'], 1))
                                 for s, x in v.items()})
    print('(ii) prediction noise', json.dumps(pred, indent=0)[:1500])
    fig_snr(wave, snr, OUT / 'maths_snr.png')
    fig_dof(dof, OUT / 'maths_dof.png')

    k_pcr = {'Tchla': 22, 'Fuco': 18}   # median PCs selected in Execution #3 (N=145)
    res_iii, curves_iii = part_iii(dRrs, pig, k_pcr)
    print('(iii)', res_iii)
    fig_filters(curves_iii, res_iii, OUT / 'maths_filter_factors.png')

    snr_summary = {name: {k: {'median': float(np.median(v[k])),
                              'frac_bands_gt1': float(np.mean(v[k] > 1))}
                          for k in v} for name, v in snr.items()}
    out = {'script': 'scripts/sdp/maths_section.py', 'epft_up_version': __version__,
           'input': data.provenance,
           'gsm_product': {'file': str(gsm_prod), 'sha256': sdpdata.sha256_of(gsm_prod)},
           'settings': {'seed': SEED, 'ell_nm': ELL, 'pct': PCT, 'steps_nm': STEPS},
           'part_i': res_i, 'implied_prior': prior,
           'part_ii': {'snr': snr_summary, 'dof': dof, 'prediction_noise': pred},
           'part_iii': res_iii,
           'run_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
    with open(OUT / 'maths_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f'wrote {OUT / "maths_summary.json"}')


if __name__ == '__main__':
    main()
