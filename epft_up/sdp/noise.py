"""Rrs noise models for the SDP track (Q&A #21).

The two models mirror IOPtics' ``ioptics/noise.py`` (which is not importable
in the ``ocean14`` environment):

* ``'pace'``: the PACE OCI per-band Rrs 1-sigma from
  ``ocpy.satellites.pace.gen_noise_vector``, the same call IOPtics makes. It
  was estimated by JXP from one early OCI L2 granule (ocpy data README) and
  is used here as white, band-independent noise unless stated otherwise.
* ``'pct:X'`` with the IOPtics two-part floor,
  ``σ = X · max(|Rrs|, median|Rrs|)``, a HyperPro-like stand-in for in-situ
  radiometry (``ioptics.noise.error_floor``, reproduced exactly).

Plus the deposit's own rounding noise (1e-6 sr⁻¹ steps, report §3.4).
"""

from __future__ import annotations

import numpy as np

#: Rounding step of the PANGAEA 937536 Rrs [sr^-1].
ROUNDING_STEP = 1e-6


def pace_sigma(wave, include_sampling=False):
    """PACE OCI per-band Rrs 1-sigma [sr^-1] on ``wave`` (ocpy; IOPtics 'pace').

    ``include_sampling=True`` rescales the per-band σ of ocpy's
    ``PACE_error.csv`` (median spacing 2 nm) to the spacing of ``wave``
    (σ × sqrt(2 nm / spacing)), so that averaging the finer bands back to the
    file's sampling recovers its σ.
    """
    from ocpy.satellites import pace
    return np.asarray(pace.gen_noise_vector(np.asarray(wave, float),
                                            include_sampling=include_sampling),
                      dtype=float)


def pct_floor_sigma(Rrs, frac):
    """IOPtics ``'pct:X'`` with the two-part floor: ``frac · max(|Rrs|, median|Rrs|)``.

    Parameters
    ----------
    Rrs : ndarray, shape (n_wave,) or (n_samples, n_wave)
    frac : float

    Returns
    -------
    ndarray, same shape
    """
    one = np.ndim(Rrs) == 1
    Rrs = np.abs(np.atleast_2d(np.asarray(Rrs, float)))
    med = np.nanmedian(Rrs, axis=1, keepdims=True)
    out = frac * np.maximum(Rrs, med)
    return out[0] if one else out


def rounding_sigma(n_wave):
    """White noise of the deposit's rounding: ``ROUNDING_STEP / sqrt(12)`` per band."""
    return np.full(n_wave, ROUNDING_STEP / np.sqrt(12.0))


# --------------------------------------------------------------------------
# Covariances and draws used by the uncertainty experiments (Execution #8)
# --------------------------------------------------------------------------
#: Correlation length [nm] of the smooth ("atmospheric-correction-like") PACE error.
PACE_ELL = 30.0
#: White fraction added to the smooth PACE error (report §5.2).
PACE_WHITE_FRAC = 0.1
#: In-situ-like fractional error.
INSITU_PCT = 0.02
MODELS = ('insitu', 'pace_white', 'pace_corr')


def _smoother(wave, width_nm=5.0):
    """Moving-mean matrix of ~``width_nm`` on ``wave`` (identity at ≥ 5 nm spacing)."""
    from epft_up.sdp import theory
    h = float(np.median(np.diff(wave)))
    w = max(1, int(round(width_nm / h)))
    w += 1 - w % 2
    return theory.boxcar_matrix(len(wave), w) if w > 1 else np.eye(len(wave))


def band_sigma(wave, Rrs, model):
    """Per-band 1-sigma before the 5 nm mean, on ``wave`` (pre-smoothing units).

    ``insitu``: ``INSITU_PCT · max(|Rrs|, median|Rrs|)`` per spectrum,
    expressed per 5 nm (× sqrt(5 nm / spacing) per band, so the 5 nm mean
    returns it). ``pace_white`` / ``pace_corr`` (white part): ocpy's PACE σ
    rescaled to the band spacing.
    """
    h = float(np.median(np.diff(wave)))
    if model == 'insitu':
        return pct_floor_sigma(Rrs, INSITU_PCT) * np.sqrt(max(5.0 / h, 1.0))
    if model in ('pace_white', 'pace_corr'):
        return pace_sigma(wave, include_sampling=True)
    raise ValueError(f'unknown noise model {model!r}')


def covariance(wave, Rrs_ref, model):
    """Noise covariance [sr⁻²] of Rrs as it reaches δRrs (after the 5 nm mean).

    Parameters
    ----------
    wave : ndarray, shape (n,)
    Rrs_ref : ndarray, shape (n,)
        Reference spectrum (sets the in-situ σ).
    model : {'insitu', 'pace_white', 'pace_corr'}

    Returns
    -------
    ndarray, shape (n, n)
    """
    from epft_up.sdp import theory
    M = _smoother(wave)
    s = band_sigma(wave, Rrs_ref, model)
    if model == 'pace_corr':
        C = theory.correlated_covariance(wave, pace_sigma(wave), PACE_ELL)
        return C + M @ np.diag((PACE_WHITE_FRAC * s)**2) @ M.T
    return M @ np.diag(s**2) @ M.T


def draw(wave, Rrs, model, n_draws, rng):
    """Noisy realizations of each spectrum: shape (n_draws, n_samples, n_wave).

    White parts are drawn per band and passed through the 5 nm mean (as the
    paper's preprocessing would see them); the smooth PACE part is drawn from
    its Gaussian-correlated covariance.
    """
    Rrs = np.atleast_2d(np.asarray(Rrs, float))
    ns, nw = Rrs.shape
    M = _smoother(wave)
    out = np.empty((n_draws, ns, nw))
    if model == 'pace_corr':
        from epft_up.sdp import theory
        C = theory.correlated_covariance(wave, pace_sigma(wave), PACE_ELL)
        Lc = np.linalg.cholesky(C + 1e-30 * np.eye(nw) + 1e-12 * np.mean(np.diag(C)) * np.eye(nw))
        s = band_sigma(wave, Rrs, model) * PACE_WHITE_FRAC
        for k in range(n_draws):
            smooth = rng.standard_normal((ns, nw)) @ Lc.T
            white = (rng.standard_normal((ns, nw)) * s) @ M.T
            out[k] = Rrs + smooth + white
        return out
    s = band_sigma(wave, Rrs, model)
    for k in range(n_draws):
        out[k] = Rrs + (rng.standard_normal((ns, nw)) * s) @ M.T
    return out
