"""Spectral preprocessing for the SDP track: smoothing and trimming.

Kramer et al. (2022, §2.2) interpolated every spectrum to 1 nm, smoothed it
with a 5 nm moving-mean bandpass filter, removed the first and last 4 nm, and
then restricted all spectra to 400-700 nm. The PANGAEA deposit is *already*
in that state (report §3.4: its 2nd-difference power sits at the rounding
floor at the 5-point boxcar's transfer-function zeros), so the reproduction
does **not** call :func:`moving_mean` on it. The functions are here to apply
the paper's recipe to new spectra (e.g. the EXPORTS-NA hold-out) and to
quantify what a second, accidental smoothing would do.

Finite differences (:func:`difference`, :func:`derivative_features`) follow
Kramer's ``Kramer_Rrs_pigments.m``: plain MATLAB ``diff`` along wavelength,
i.e. ``diff(y, 2)`` = y[i+1] - 2 y[i] + y[i-1] (Catlett & Siegel 2018 Eq. 2
times Δλ²; Δλ = 1 nm) and a *forward* first difference ``diff(y, 1)``. The
alternative transforms are added in Execution #5.
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


def difference(wave, y, order=2, step=1, scale=False):
    """MATLAB-style ``diff`` along wavelength, optionally on a coarser grid.

    Parameters
    ----------
    wave : ndarray, shape (n,)
        Uniform wavelength grid [nm].
    y : ndarray, shape (..., n)
    order : {1, 2}, optional
        ``1``: forward difference ``y[i+1] - y[i]`` (placed at the midpoint);
        ``2``: ``y[i+1] - 2 y[i] + y[i-1]`` (placed at ``wave[i]``).
    step : int, optional
        Subsample every ``step``-th band first (5 or 10 for the paper's
        5 and 10 nm variants).
    scale : bool, optional
        Divide by Δλ**order (derivative units). Kramer does not; it does not
        matter for z-scored predictors.

    Returns
    -------
    wave_d : ndarray
    d : ndarray, shape (..., n_d)
    """
    wave = np.asarray(wave, dtype=float)[::step]
    y = np.asarray(y, dtype=float)[..., ::step]
    dl = wave[1] - wave[0]
    if order == 1:
        d = np.diff(y, n=1, axis=-1)
        w = 0.5 * (wave[1:] + wave[:-1])
    elif order == 2:
        d = np.diff(y, n=2, axis=-1)
        w = wave[1:-1]
    else:
        raise ValueError('order must be 1 or 2')
    if scale:
        d = d / dl**order
    return w, d


def derivative_features(wave, y, orders=(2,), step=1):
    """Concatenate difference spectra of the given orders as predictors.

    ``orders=(2,)`` is Kramer's δRrs'' input; ``orders=(1, 2)`` is the
    "Rrs' + Rrs''" supplementary variant (Catlett & Siegel 2018 style).

    Returns
    -------
    labels : list of str
        ``'d<order>_<wavelength>'`` for each feature column.
    X : ndarray, shape (n_samples, n_features)
    """
    labels, blocks = [], []
    for o in orders:
        w, d = difference(wave, y, order=o, step=step)
        labels += [f'd{o}_{x:g}' for x in w]
        blocks.append(np.atleast_2d(d))
    return labels, np.hstack(blocks)
