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
