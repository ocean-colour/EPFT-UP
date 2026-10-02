"""Spectral preprocessing for the SDP track: smoothing and trimming.

Kramer et al. (2022, §2.2) interpolated every spectrum to 1 nm, smoothed it
with a 5 nm moving-mean bandpass filter, removed the first and last 4 nm, and
then restricted all spectra to 400-700 nm. The PANGAEA deposit is *already*
in that state (report §3.4: its 2nd-difference power sits at the rounding
floor at the 5-point boxcar's transfer-function zeros), so the reproduction
does **not** call :func:`moving_mean` on it. The functions are here to apply
the paper's recipe to new spectra (e.g. the EXPORTS-NA hold-out) and to
quantify what a second, accidental smoothing would do.

Finite differences and the alternative transforms are added in Execution #3
and #5.
"""

from __future__ import annotations

import numpy as np


def moving_mean(y, width=5):
    """Centred moving mean along the last axis; edges are returned as NaN.

    Equivalent to a ``width``-point boxcar on a uniform grid. Points within
    ``width // 2`` of either end have an incomplete window and are set to NaN
    rather than shrunk or padded, so trimming them (:func:`trim_edges`) is
    explicit.

    Parameters
    ----------
    y : ndarray, shape (..., n)
    width : int, optional
        Odd window length in samples (5 = 5 nm on the 1 nm grid).

    Returns
    -------
    ndarray, shape (..., n)
    """
    if width < 1 or width % 2 == 0:
        raise ValueError(f'moving_mean width must be a positive odd integer, got {width}')
    y = np.asarray(y, dtype=float)
    half = width // 2
    c = np.cumsum(np.concatenate([np.zeros(y.shape[:-1] + (1,)), y], axis=-1),
                  axis=-1)
    out = np.full_like(y, np.nan)
    out[..., half:y.shape[-1] - half] = (c[..., width:] - c[..., :-width]) / width
    return out


def trim_edges(wave, y, n_trim=4):
    """Drop ``n_trim`` samples from each end of the spectral axis.

    Parameters
    ----------
    wave : ndarray, shape (n,)
    y : ndarray, shape (..., n)
    n_trim : int, optional
        Kramer 2022 removed 4 nm after smoothing.

    Returns
    -------
    wave_t : ndarray, shape (n - 2 n_trim,)
    y_t : ndarray, shape (..., n - 2 n_trim)
    """
    sl = slice(n_trim, len(wave) - n_trim)
    return np.asarray(wave)[sl], np.asarray(y)[..., sl]


def kramer_preprocess(wave, Rrs, wave_out=None, width=5, n_trim=4):
    """The Kramer 2022 §2.2 recipe for spectra that are *not* yet processed.

    Linear interpolation to 1 nm, a ``width``-nm moving mean, trimming
    ``n_trim`` nm from each end, then restriction to ``wave_out``.

    Parameters
    ----------
    wave : ndarray, shape (n,)
        Native wavelengths [nm] (need not be uniform).
    Rrs : ndarray, shape (n_samples, n)
    wave_out : ndarray, optional
        Output grid; default 400-700 nm at 1 nm.
    width, n_trim : int, optional

    Returns
    -------
    ndarray, shape (n_samples, len(wave_out))
        NaN where ``wave_out`` falls outside the trimmed native range.
    """
    wave = np.asarray(wave, dtype=float)
    Rrs = np.atleast_2d(np.asarray(Rrs, dtype=float))
    if wave_out is None:
        wave_out = np.arange(400.0, 701.0)
    grid = np.arange(np.ceil(wave[0]), np.floor(wave[-1]) + 1.0)
    on_grid = np.vstack([np.interp(grid, wave, r) for r in Rrs])
    if n_trim < width // 2:
        raise ValueError('n_trim must cover the incomplete moving-mean windows '
                         f'(>= {width // 2} for width={width})')
    sm = moving_mean(on_grid, width)
    # Trim from the 1 nm grid's ends; this also removes the NaN edge points.
    g_t, sm_t = trim_edges(grid, sm, n_trim)
    out = np.full((Rrs.shape[0], len(wave_out)), np.nan)
    ok = (wave_out >= g_t[0]) & (wave_out <= g_t[-1])
    for i, r in enumerate(sm_t):
        out[i, ok] = np.interp(wave_out[ok], g_t, r)
    return out
