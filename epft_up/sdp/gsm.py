"""The Kramer et al. (2022) GSM-like reflectance model and the residual δRrs.

Independent implementation of Kramer 2022 §2.3, Eqs. 1-6. For each spectrum:

* above- to below-water reflectance (Lee et al. 2002)::

      rrs = Rrs / (0.52 + 1.7 Rrs)

* absorption ``a = a_w + a_ph + a_dg`` with
  ``a_ph = A(λ) Tchla^B(λ)`` (NOMAD regression; reference table) and
  ``a_dg = a_dg(443) exp(-S_dg (λ - 443))``, where (Carder et al. 1999)
  ``S_dg = 0.01447 + 0.00033 Rrs(490)/Rrs(555)``;
* backscattering ``bb = b_bw + b_bp`` with ``b_bw = b_sw / 2`` from Zhang et
  al. (2009) at the sample's T and S, and
  ``b_bp = b_bp(443) (443/λ)^η``, where (Lee et al. 2002)
  ``η = 2 (1 - 1.2 exp(-0.9 rrs(440)/rrs(555)))``;
* Gordon et al. (1988), Eq. 1::

      u = bb / (a + bb),   rrs_mod = 0.0949 u + 0.0794 u**2

* the three free parameters (Tchla, a_dg(443), b_bp(443)) are fitted by
  unweighted least squares in rrs, and the residual is
  ``δRrs = Rrs - Rrs_mod`` with ``Rrs_mod = 0.52 rrs_mod / (1 - 1.7 rrs_mod)``.

Points where the paper's text and the reference code disagree, and what is
done here (each verified in ``scripts/sdp/reproduce_gsm.py``):

1. **Sign of S_dg.** The paper's Eq. 5 reads
   ``S_dg = -0.01447 + 0.00033 Rrs490/Rrs555`` inside ``exp(S_dg (λ-443))``.
   Both reference codes use ``exp(-(0.01447 + 0.00033 r)(λ-443))``, which is
   Carder et al. (1999). We follow Carder and the code; the printed Eq. 5 has
   a sign typo.
2. **η band ratio.** The MATLAB code and Lee et al. (2002) use the
   *below-water* ratio rrs(440)/rrs(555); the Python port uses above-water
   Rrs. We use rrs (``eta_from='rrs'``); the port's choice is an option.
3. **Fit.** The reference is MATLAB ``fminsearch`` (Nelder-Mead) from
   (0.15, 0.01, 0.0029) with TolX = TolFun = 1e-9 and 2000 iterations and
   evaluations, unbounded; a non-converged fit is set to NaN. We reproduce it
   with scipy's Nelder-Mead, which uses the same 5% initial simplex
   (``method='kramer'``). The only change is that a trial with Tchla ≤ 0
   (where ``Tchla**B`` is complex in MATLAB and NaN in numpy) is given an
   infinite cost.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize

from epft_up.sdp import data as sdpdata

#: Gordon et al. (1988) quadratic coefficients (Kramer 2022 Eq. 1).
G0, G1 = 0.0949, 0.0794
#: Initial guess of the reference fit: Tchla [mg m^-3], a_dg(443), b_bp(443) [m^-1].
X0 = (0.15, 0.01, 0.0029)
#: Free-parameter names, in order.
PARAMS = ('Tchla', 'adg443', 'bbp443')


# --------------------------------------------------------------------------
# Reflectance conversions (Lee et al. 2002)
# --------------------------------------------------------------------------
def Rrs_to_rrs(Rrs):
    """Above-water Rrs -> below-water rrs (Lee et al. 2002), Kramer Eq. 2."""
    Rrs = np.asarray(Rrs, dtype=float)
    return Rrs / (0.52 + 1.7 * Rrs)


def rrs_to_Rrs(rrs):
    """Below-water rrs -> above-water Rrs, the exact inverse of :func:`Rrs_to_rrs`."""
    rrs = np.asarray(rrs, dtype=float)
    return 0.52 * rrs / (1.0 - 1.7 * rrs)


# --------------------------------------------------------------------------
# Pure seawater scattering: Zhang, Hu & He (2009, Opt. Express 17, 5698)
# --------------------------------------------------------------------------
_NA = 6.0221417930e23     # Avogadro [mol^-1]
_KB = 1.3806503e-23       # Boltzmann [J K^-1]
_M0 = 18e-3               # molar mass of water [kg mol^-1]


def _refractive_index_sw(wave, T, S):
    """Seawater refractive index and dn/dS (Quan & Fry 1995; Ciddor 1996 air)."""
    lam_um = wave / 1e3
    n_air = 1.0 + (5792105.0 / (238.0185 - lam_um**-2)
                   + 167917.0 / (57.362 - lam_um**-2)) / 1e8
    n0, n1, n2, n3, n4 = 1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6
    n5, n6, n7, n8, n9 = 15.868, 0.01155, -0.00423, -4382.0, 1.1455e6
    nsw = (n0 + (n1 + n2 * T + n3 * T**2) * S + n4 * T**2
           + (n5 + n6 * S + n7 * T) / wave + n8 / wave**2 + n9 / wave**3)
    dnds = (n1 + n2 * T + n3 * T**2 + n6 / wave) * n_air
    return nsw * n_air, dnds


def _isothermal_compressibility(T, S):
    """Seawater isothermal compressibility [Pa^-1] (Millero 1980 secant bulk)."""
    kw = (19652.21 + 148.4206 * T - 2.327105 * T**2 + 1.360477e-2 * T**3
          - 5.155288e-5 * T**4)
    a0 = 54.6746 - 0.603459 * T + 1.09987e-2 * T**2 - 6.167e-5 * T**3
    b0 = 7.944e-2 + 1.6483e-2 * T - 5.3009e-4 * T**2
    ks = kw + a0 * S + b0 * S**1.5
    return 1.0 / ks * 1e-5


def _density_sw(T, S):
    """Seawater density [kg m^-3] (UNESCO 1981)."""
    rho_w = (999.842594 + 6.793952e-2 * T - 9.09529e-3 * T**2 + 1.001685e-4 * T**3
             - 1.120083e-6 * T**4 + 6.536332e-9 * T**5)
    return (rho_w
            + (8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T**2 - 8.2467e-7 * T**3
               + 5.3875e-9 * T**4) * S
            + (-5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2) * S**1.5
            + 4.8314e-4 * S**2)


def _dln_water_activity_dS(T, S):
    """d ln(a_w) / dS of seawater (Millero & Leung 1976 fit)."""
    return ((-5.58651e-4 + 2.40452e-7 * T - 3.12165e-9 * T**2 + 2.40808e-11 * T**3)
            + 1.5 * (1.79613e-5 - 9.9422e-8 * T + 2.08919e-9 * T**2
                     - 1.39872e-11 * T**3) * S**0.5
            + 2.0 * (-2.31065e-6 - 1.37674e-9 * T - 1.93316e-11 * T**2) * S)


def _pmh(n):
    """Density derivative of the refractive index (Proutiere-Megnassah-Hermant)."""
    n2 = n**2
    return (n2 - 1.0) * (1.0 + 2.0 / 3.0 * (n2 + 2.0) * (n / 3.0 - 1.0 / 3.0 / n)**2)


def bsw_zhang2009(wave, T, S, delta=0.039):
    """Total scattering coefficient of pure seawater, Zhang et al. (2009).

    Parameters
    ----------
    wave : ndarray, shape (n_wave,)
        Wavelength [nm].
    T, S : float or ndarray, shape (n_samples,)
        Temperature [degC] and salinity [PSU].
    delta : float, optional
        Depolarization ratio (Farinato & Roswell 1976).

    Returns
    -------
    ndarray, shape (n_wave,) or (n_samples, n_wave)
        ``b_sw`` [m^-1]; backscattering is half of it.
    """
    wave = np.asarray(wave, dtype=float)
    scalar = np.ndim(T) == 0 and np.ndim(S) == 0
    T = np.atleast_1d(np.asarray(T, dtype=float))[:, None]
    S = np.atleast_1d(np.asarray(S, dtype=float))[:, None]
    lam = wave[None, :]

    nsw, dnds = _refractive_index_sw(lam, T, S)
    iso = _isothermal_compressibility(T, S)
    rho = _density_sw(T, S)
    dlnaw = _dln_water_activity_dS(T, S)
    dfri = _pmh(nsw)
    cabannes = (6.0 + 6.0 * delta) / (6.0 - 7.0 * delta)
    beta_df = (np.pi**2 / 2.0 * (lam * 1e-9)**-4 * _KB * (T + 273.15) * iso
               * dfri**2 * cabannes)
    flu_con = S * _M0 * dnds**2 / rho / (-dlnaw) / _NA
    beta_cf = 2.0 * np.pi**2 * (lam * 1e-9)**-4 * nsw**2 * flu_con * cabannes
    bsw = 8.0 * np.pi / 3.0 * (beta_df + beta_cf) * (2.0 + delta) / (1.0 + delta)
    return bsw[0] if scalar else bsw


# --------------------------------------------------------------------------
# Reference tables (A, B for a_ph; a_w)
# --------------------------------------------------------------------------
@dataclass
class RefTables:
    """Spectral constants on the model grid.

    Attributes
    ----------
    wave : ndarray
        Wavelength [nm].
    A, B : ndarray
        ``a_ph(λ) = A(λ) Tchla^B(λ)`` (Kramer 2022 Table S6, NOMAD fit).
    aw : ndarray
        Pure-water absorption [m^-1] (file ``aw_mcf16_350_700_1nm``; Mason et
        al. 2016 to 550 nm).
    source : dict
        Files and the reference-repo commit they came from.
    """
    wave: np.ndarray
    A: np.ndarray
    B: np.ndarray
    aw: np.ndarray
    source: dict


def load_ref_tables(wave=None, ref_dir=None):
    """Load A, B and a_w from the pinned ``Rrs_pigments`` checkout.

    Parameters
    ----------
    wave : ndarray, optional
        Output grid [nm]; default :data:`epft_up.sdp.data.WAVE`. Must be a
        subset of the tables' 1 nm grid (350-700 nm).
    ref_dir : str or Path, optional
        Default ``$OS_COLOR/PANGAEA/Kramer2022/ref/Rrs_pigments``.

    Returns
    -------
    RefTables
    """
    wave = sdpdata.WAVE if wave is None else np.asarray(wave, dtype=float)
    ref_dir = (Path(ref_dir) if ref_dir is not None
               else sdpdata.kramer2022_dir() / 'ref' / 'Rrs_pigments')
    f_ab = ref_dir / 'aph_A_B_Coeffs_Sasha_RSE_paper.txt'
    f_aw = ref_dir / 'aw_mcf16_350_700_1nm.txt'
    ab = np.loadtxt(f_ab, delimiter=',', skiprows=1)
    aw = np.loadtxt(f_aw, skiprows=1)
    idx_ab = np.searchsorted(ab[:, 0], wave)
    idx_aw = np.searchsorted(aw[:, 0], wave)
    if not (np.array_equal(ab[idx_ab, 0], wave) and np.array_equal(aw[idx_aw, 0], wave)):
        raise ValueError('requested wavelengths are not on the reference tables\' grid')
    status = sdpdata.ref_repo_status(ref_dir.parent).get('Rrs_pigments', {})
    return RefTables(
        wave=wave, A=ab[idx_ab, 1], B=ab[idx_ab, 2], aw=aw[idx_aw, 1],
        source={'A_B': str(f_ab), 'aw': str(f_aw),
                'sha256_A_B': sdpdata.sha256_of(f_ab),
                'sha256_aw': sdpdata.sha256_of(f_aw),
                'repo_sha': status.get('sha')})


# --------------------------------------------------------------------------
# Spectral shapes and the forward model
# --------------------------------------------------------------------------
def _at(wave, y, w0):
    """Column of ``y`` at wavelength ``w0`` (must be on the grid)."""
    i = int(np.searchsorted(wave, w0))
    if i >= len(wave) or wave[i] != w0:
        raise ValueError(f'{w0} nm is not on the wavelength grid')
    return y[..., i]


def sdg_slope(Rrs490, Rrs555):
    """CDOM+NAP spectral slope [nm^-1] (Carder et al. 1999; see note 1)."""
    return 0.01447 + 0.00033 * np.asarray(Rrs490) / np.asarray(Rrs555)


def eta_lee2002(r440, r555):
    """Particulate backscattering exponent (Lee et al. 2002)."""
    return 2.0 * (1.0 - 1.2 * np.exp(-0.9 * np.asarray(r440) / np.asarray(r555)))


@dataclass
class Shapes:
    """Per-sample fixed spectral shapes of the GSM-like model.

    All arrays are (n_samples, n_wave) except ``A``, ``B``, ``aw`` (n_wave,).
    """
    wave: np.ndarray
    aw: np.ndarray
    A: np.ndarray
    B: np.ndarray
    bbw: np.ndarray
    adg_star: np.ndarray
    bbp_star: np.ndarray
    Sdg: np.ndarray
    eta: np.ndarray


def build_shapes(wave, Rrs, T, S, tables=None, eta_from='rrs'):
    """Fixed spectral shapes for each spectrum (Kramer 2022 Eqs. 3-6).

    Parameters
    ----------
    wave : ndarray, shape (n_wave,)
    Rrs : ndarray, shape (n_samples, n_wave)
        Above-water reflectance [sr^-1].
    T, S : ndarray, shape (n_samples,)
        Temperature [degC] and salinity [PSU] for b_bw.
    tables : RefTables, optional
        Default :func:`load_ref_tables` on ``wave``.
    eta_from : {'rrs', 'Rrs'}, optional
        Band ratio for η: below-water rrs (MATLAB code, Lee 2002; default) or
        above-water Rrs (the Python port).

    Returns
    -------
    Shapes
    """
    wave = np.asarray(wave, dtype=float)
    Rrs = np.atleast_2d(np.asarray(Rrs, dtype=float))
    tables = load_ref_tables(wave) if tables is None else tables
    rrs = Rrs_to_rrs(Rrs)

    Sdg = sdg_slope(_at(wave, Rrs, 490.0), _at(wave, Rrs, 555.0))
    ref = rrs if eta_from == 'rrs' else Rrs
    if eta_from not in ('rrs', 'Rrs'):
        raise ValueError("eta_from must be 'rrs' or 'Rrs'")
    eta = eta_lee2002(_at(wave, ref, 440.0), _at(wave, ref, 555.0))

    adg_star = np.exp(-Sdg[:, None] * (wave[None, :] - 443.0))
    bbp_star = (443.0 / wave[None, :])**eta[:, None]
    bbw = 0.5 * bsw_zhang2009(wave, T, S)
    bbw = np.broadcast_to(bbw, Rrs.shape).copy()
    return Shapes(wave=wave, aw=tables.aw, A=tables.A, B=tables.B, bbw=bbw,
                  adg_star=adg_star, bbp_star=bbp_star, Sdg=Sdg, eta=eta)


def forward_rrs(params, aw, A, B, bbw, adg_star, bbp_star):
    """Below-water rrs of the GSM-like model for one spectrum.

    Parameters
    ----------
    params : sequence of 3 floats
        (Tchla, a_dg(443), b_bp(443)).
    aw, A, B, bbw, adg_star, bbp_star : ndarray, shape (n_wave,)

    Returns
    -------
    ndarray, shape (n_wave,)
    """
    chl, adg443, bbp443 = params
    a = aw + A * chl**B + adg443 * adg_star
    bb = bbw + bbp443 * bbp_star
    u = bb / (a + bb)
    return (G0 + G1 * u) * u


def _cost(params, rrs_obs, aw, A, B, bbw, adg_star, bbp_star):
    """Unweighted sum of squares in rrs (``gsm_cost.m``); +inf for Tchla <= 0."""
    if params[0] <= 0.0:
        return np.inf
    r = rrs_obs - forward_rrs(params, aw, A, B, bbw, adg_star, bbp_star)
    return float(np.dot(r, r))


@dataclass
class GSMFit:
    """Result of fitting the GSM-like model to a set of spectra.

    Attributes
    ----------
    wave : ndarray, shape (n_wave,)
    params : ndarray, shape (n_samples, 3)
        Fitted (Tchla, a_dg(443), b_bp(443)); NaN where the fit failed.
    converged : ndarray of bool, shape (n_samples,)
    n_iter : ndarray of int, shape (n_samples,)
    cost : ndarray, shape (n_samples,)
    rrs, rrs_mod : ndarray, shape (n_samples, n_wave)
        Observed and modelled below-water reflectance.
    Rrs, Rrs_mod, dRrs : ndarray, shape (n_samples, n_wave)
        Observed and modelled above-water Rrs and the residual δRrs.
    shapes : Shapes
    settings : dict
    """
    wave: np.ndarray
    params: np.ndarray
    converged: np.ndarray
    n_iter: np.ndarray
    cost: np.ndarray
    rrs: np.ndarray
    rrs_mod: np.ndarray
    Rrs: np.ndarray
    Rrs_mod: np.ndarray
    dRrs: np.ndarray
    shapes: Shapes
    settings: dict


def _fit_one_kramer(rrs_obs, args, x0, tol, maxiter):
    """Nelder-Mead as in ``gsm_invert.m``."""
    res = minimize(_cost, np.asarray(x0, dtype=float), args=(rrs_obs, *args),
                   method='Nelder-Mead',
                   options={'xatol': tol, 'fatol': tol, 'maxiter': maxiter,
                            'maxfev': maxiter})
    return res.x, bool(res.success), int(res.nit)


def _fit_one_lsq(rrs_obs, args, x0):
    """Bounded trust-region least squares in log-parameters (robustness check)."""
    def resid(logp):
        return rrs_obs - forward_rrs(np.exp(logp), *args)
    res = least_squares(resid, np.log(np.asarray(x0, dtype=float)), method='trf',
                        x_scale='jac', xtol=1e-12, ftol=1e-12, gtol=1e-12,
                        max_nfev=5000)
    return np.exp(res.x), bool(res.success), int(res.nfev)


def fit_gsm(wave, Rrs, T, S, method='kramer', x0=X0, tol=1e-9, maxiter=2000,
            tables=None, eta_from='rrs'):
    """Fit the GSM-like model to each spectrum and form δRrs.

    Parameters
    ----------
    wave : ndarray, shape (n_wave,)
    Rrs : ndarray, shape (n_samples, n_wave)
        Above-water reflectance [sr^-1].
    T, S : ndarray, shape (n_samples,)
    method : {'kramer', 'lsq'}, optional
        ``'kramer'``: Nelder-Mead exactly as ``gsm_invert.m`` (non-converged
        fits -> NaN). ``'lsq'``: bounded least squares on log-parameters, used
        to test the sensitivity of δRrs to the optimizer.
    x0 : sequence of 3 floats, optional
    tol : float, optional
        Nelder-Mead xatol and fatol (MATLAB TolX, TolFun).
    maxiter : int, optional
        Nelder-Mead maxiter and maxfev (MATLAB MaxIter, MaxFunEvals).
    tables : RefTables, optional
    eta_from : {'rrs', 'Rrs'}, optional

    Returns
    -------
    GSMFit
    """
    wave = np.asarray(wave, dtype=float)
    Rrs = np.atleast_2d(np.asarray(Rrs, dtype=float))
    T = np.broadcast_to(np.asarray(T, dtype=float), (Rrs.shape[0],))
    S = np.broadcast_to(np.asarray(S, dtype=float), (Rrs.shape[0],))
    tables = load_ref_tables(wave) if tables is None else tables
    sh = build_shapes(wave, Rrs, T, S, tables=tables, eta_from=eta_from)
    rrs = Rrs_to_rrs(Rrs)

    n = Rrs.shape[0]
    params = np.full((n, 3), np.nan)
    conv = np.zeros(n, dtype=bool)
    nit = np.zeros(n, dtype=int)
    cost = np.full(n, np.nan)
    rrs_mod = np.full_like(rrs, np.nan)
    for i in range(n):
        args = (sh.aw, sh.A, sh.B, sh.bbw[i], sh.adg_star[i], sh.bbp_star[i])
        if method == 'kramer':
            p, ok, k = _fit_one_kramer(rrs[i], args, x0, tol, maxiter)
        elif method == 'lsq':
            p, ok, k = _fit_one_lsq(rrs[i], args, x0)
        else:
            raise ValueError(f"unknown method {method!r}")
        conv[i], nit[i] = ok, k
        if ok:
            params[i] = p
            rrs_mod[i] = forward_rrs(p, *args)
            cost[i] = _cost(p, rrs[i], *args)

    Rrs_mod = rrs_to_Rrs(rrs_mod)
    return GSMFit(wave=wave, params=params, converged=conv, n_iter=nit, cost=cost,
                  rrs=rrs, rrs_mod=rrs_mod, Rrs=Rrs, Rrs_mod=Rrs_mod,
                  dRrs=Rrs - Rrs_mod, shapes=sh,
                  settings={'method': method, 'x0': list(map(float, x0)),
                            'tol': tol, 'maxiter': maxiter, 'eta_from': eta_from,
                            'tables': tables.source})
