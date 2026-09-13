"""The MOANA retrieval (design doc §3): Rrs (+SST) → three cell abundances.

A fixed linear map, implemented exactly as the operational OCSSW code
(`get_Cpicophyt.c`) applies it, with two deliberate departures:

- **raw floats + QC flags** instead of NASA's silent clamp-to-zero and int32
  truncation (Q&A #11). ``nasa_compat=True`` restores NASA's behaviour for the
  bit-exactness test and nothing else;
- the disputed PC assignment is a **flag** (`pc_mapping`): ``'operational'``
  uses `picophyt.json`'s slots as shipped; ``'atbd'`` moves the last
  Prochlorococcus coefficient U7→U17 and the last Synechococcus coefficient
  U16→U13, as the ATBD v1.2 equations have them (design §4.1, report §7.1).

Free invariants any caller can rely on (asserted in the test suite):
after :func:`standardize`, every spectrum has mean 0 and L2 norm
``sqrt(123) ≈ 11.0905`` exactly, so every score satisfies ``|U_i| ≤ 11.0905``.
"""

from __future__ import annotations

import numpy as np

from epft_up.moana.io import load_luts

# --- QC flag bits (design §3; Q&A #11) ---------------------------------------
FLAG_NEGATIVE_PRO = 1       # Pro regression went negative (kept as-is, flagged)
FLAG_EXTRAPOLATED = 2       # input grid did not bracket 414-660 nm
FLAG_TOO_FEW_BANDS = 4      # fewer valid input bands than min_bands
FLAG_BAD_SST = 8            # SST missing/non-positive where Pro needs it
FLAG_NONFINITE_RRS = 16     # non-finite values inside the fit range

#: Number of MOANA wavelengths; norm of every standardised spectrum is
#: sqrt(N_BANDS - 1).
N_BANDS = 124


def _coef_matrix(lut, pc_mapping='operational'):
    """Expand the zero-padded JSON coefficient vectors onto the 45-PC axis.

    Parameters
    ----------
    lut : dict — from :func:`epft_up.moana.io.load_luts`.
    pc_mapping : {'operational', 'atbd'}
        'operational' consumes the JSON slots as `get_Cpicophyt.c` does
        (slot i multiplies U_i). 'atbd' relocates the two disputed
        coefficients: Pro's U7 value → U17, Syn's U16 value → U13.

    Returns
    -------
    dict with keys ``pro``, ``syn``, ``apeuk``; each value is
    ``(intercept, sst_coef_or_None, w)`` where ``w`` is a (45,) float64
    vector of per-PC weights (zeros where the taxon has no term).
    """
    npc, n_stored = lut['npc'], lut['V'].shape[1]

    def expand(slots, sst_slot):
        # slots: JSON vector; slots[0]=intercept, then optionally SST, then U1..
        w = np.zeros(n_stored)
        off = 2 if sst_slot else 1
        w[:npc] = slots[off:off + npc]
        return (slots[0], slots[1] if sst_slot else None, w)

    pro = expand(lut['pro_coef'], sst_slot=True)
    syn = expand(lut['syn_coef'], sst_slot=False)
    apeuk = expand(lut['apeuk_coef'], sst_slot=False)

    if pc_mapping == 'atbd':
        # Same values, different slots (report §7.1). Indices are 0-based.
        for tax, src, dst in ((pro, 7 - 1, 17 - 1), (syn, 16 - 1, 13 - 1)):
            w = tax[2]
            w[dst], w[src] = w[src], 0.0
    elif pc_mapping != 'operational':
        raise ValueError(f"pc_mapping must be 'operational' or 'atbd', "
                         f"got {pc_mapping!r}")
    return {'pro': pro, 'syn': syn, 'apeuk': apeuk}


def interp_to_moana(wave_in, rrs_in, wave_out, min_bands=10):
    """Linearly interpolate spectra onto the MOANA grid, flagging edge cases.

    Parameters
    ----------
    wave_in : (w,) array — input wavelengths [nm], ascending.
    rrs_in : (n, w) or (w,) array — Rrs [sr⁻¹]; NaNs are treated as invalid
        bands and dropped per spectrum before interpolation.
    wave_out : (124,) array — target grid (the LUT's 414..660 @ 2 nm).
    min_bands : int, optional — minimum valid input bands per spectrum;
        below this the spectrum is all-NaN + FLAG_TOO_FEW_BANDS (the OCSSW
        code has no such guard; we refuse to interpolate noise).

    Returns
    -------
    (rrs_out, flags) : ((n, 124) float64, (n,) int) — constant-value
        extrapolation beyond the input range (as np.interp does) is allowed
        but flagged FLAG_EXTRAPOLATED.
    """
    rrs_in = np.atleast_2d(np.asarray(rrs_in, dtype=np.float64))
    wave_in = np.asarray(wave_in, dtype=np.float64)
    n = rrs_in.shape[0]
    out = np.full((n, wave_out.size), np.nan)
    flags = np.zeros(n, dtype=int)
    for i in range(n):
        ok = np.isfinite(rrs_in[i])
        if ok.sum() < min_bands:
            flags[i] |= FLAG_TOO_FEW_BANDS
            continue
        w, r = wave_in[ok], rrs_in[i, ok]
        if w[0] > wave_out[0] or w[-1] < wave_out[-1]:
            flags[i] |= FLAG_EXTRAPOLATED
        out[i] = np.interp(wave_out, w, r)
    return out, flags


def standardize(rrs):
    """Standardise each spectrum over wavelength (Lange Eq. 3; design §3).

    ``Rrs' = (Rrs − mean(Rrs)) / sd(Rrs)`` with the **sample (N−1)** standard
    deviation (Q&A #8, matching OCSSW). Rows with any non-finite value or
    zero spread return NaN.

    Parameters
    ----------
    rrs : (n, w) or (w,) array — spectra on the MOANA grid.

    Returns
    -------
    (n, w) float64 — standardised spectra; each finite row has mean 0 and
    L2 norm ``sqrt(w − 1)``.
    """
    rrs = np.atleast_2d(np.asarray(rrs, dtype=np.float64))
    mu = rrs.mean(axis=1, keepdims=True)
    sd = rrs.std(axis=1, ddof=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        out = (rrs - mu) / sd
    out[~np.isfinite(out).all(axis=1)] = np.nan
    return out


def project(rrs_std, V):
    """Project standardised spectra onto the PCA loadings: ``U = Rrs' · V``.

    Parameters
    ----------
    rrs_std : (n, w) array — from :func:`standardize`.
    V : (w, p) array — loading matrix (LUT ``V`` or a retrained basis).

    Returns
    -------
    (n, p) float64 — scores; NaN rows propagate.
    """
    return np.atleast_2d(rrs_std) @ np.asarray(V)


def reconstruction_residual(rrs_std, V):
    """L2 norm of the part of ``Rrs'`` outside the span of ``V`` (report §9.3).

    A cheap out-of-domain score: spectra unlike the training set cannot be
    represented by the truncated basis and get a large residual. Computed for
    every retrieval because it is one matrix multiply.

    Parameters
    ----------
    rrs_std : (n, w) array — standardised spectra.
    V : (w, p) array — orthonormal loading matrix.

    Returns
    -------
    (n,) float64 — ``‖Rrs' − V Vᵀ Rrs'‖₂`` per spectrum.
    """
    rrs_std = np.atleast_2d(rrs_std)
    recon = (rrs_std @ V) @ np.asarray(V).T
    return np.linalg.norm(rrs_std - recon, axis=1)


def evaluate(scores, sst=None, lut=None, pc_mapping='operational'):
    """Evaluate the three regressions on PC scores (design §3 step 5).

    Parameters
    ----------
    scores : (n, p) array — from :func:`project`; p ≥ npc.
    sst : (n,) array or scalar, optional — SST [°C]; required for a finite
        Prochlorococcus value (enters as log10(SST)).
    lut : dict, optional — from :func:`load_luts`; loaded if omitted.
    pc_mapping : {'operational', 'atbd'} — see :func:`_coef_matrix`.

    Returns
    -------
    dict with (n,) float64 arrays ``pro``, ``syn``, ``apeuk`` (cells mL⁻¹;
    raw floats — Pro may be negative, flagged not clamped) and (n,) int
    ``flags``.
    """
    lut = lut if lut is not None else load_luts()
    coefs = _coef_matrix(lut, pc_mapping)
    scores = np.atleast_2d(np.asarray(scores, dtype=np.float64))
    n = scores.shape[0]
    flags = np.zeros(n, dtype=int)

    # Prochlorococcus: identity link, needs log10(SST).
    b0, b_sst, w = coefs['pro']
    if sst is None:
        log_sst = np.full(n, np.nan)
        flags |= FLAG_BAD_SST
    else:
        sst = np.broadcast_to(np.asarray(sst, dtype=np.float64), (n,)).copy()
        bad = ~np.isfinite(sst) | (sst <= 0)
        flags[bad] |= FLAG_BAD_SST
        sst[bad] = np.nan
        log_sst = np.log10(sst)
    pro = b0 + b_sst * log_sst + scores @ w
    flags[np.isfinite(pro) & (pro < 0)] |= FLAG_NEGATIVE_PRO

    # Synechococcus / picoeukaryotes: log10 link, exponentiated.
    b0, _, w = coefs['syn']
    syn = 10.0 ** (b0 + scores @ w)
    b0, _, w = coefs['apeuk']
    apeuk = 10.0 ** (b0 + scores @ w)

    return {'pro': pro, 'syn': syn, 'apeuk': apeuk, 'flags': flags}


def run_moana(wave_in, rrs_in, sst=None, lut=None, pc_mapping='operational',
              nasa_compat=False, min_bands=10):
    """Full MOANA retrieval on one or many spectra (design §3, steps 1–5).

    Parameters
    ----------
    wave_in : (w,) array — input wavelengths [nm].
    rrs_in : (n, w) or (w,) array — Rrs [sr⁻¹].
    sst : (n,) array or scalar, optional — SST [°C] for Prochlorococcus.
    lut : dict, optional — from :func:`load_luts`; loaded if omitted.
    pc_mapping : {'operational', 'atbd'} — disputed-coefficient placement.
    nasa_compat : bool, optional
        If True, reproduce OCSSW post-processing exactly: negative Pro
        clamped to 0 and all three outputs truncated to int32. Off by
        default (Q&A #11) — use only for the bit-exactness test.
    min_bands : int, optional — see :func:`interp_to_moana`.

    Returns
    -------
    dict — ``pro``, ``syn``, ``apeuk`` ((n,) arrays, cells mL⁻¹),
    ``flags`` ((n,) int bitmask), ``scores`` ((n, 45)),
    ``recon_residual`` ((n,) out-of-domain score),
    ``rrs_std`` ((n, 124) standardised spectra).
    """
    lut = lut if lut is not None else load_luts()
    rrs124, flags_i = interp_to_moana(wave_in, rrs_in, lut['wave'], min_bands)
    bad = ~np.isfinite(rrs124).all(axis=1)
    flags_i[bad & (flags_i == 0)] |= FLAG_NONFINITE_RRS
    rrs_std = standardize(rrs124)
    scores = project(rrs_std, lut['V'])
    out = evaluate(scores, sst=sst, lut=lut, pc_mapping=pc_mapping)
    out['flags'] |= flags_i
    out['scores'] = scores
    out['recon_residual'] = reconstruction_residual(rrs_std, lut['V'])
    out['rrs_std'] = rrs_std
    if nasa_compat:
        # OCSSW: clamp negative Pro (only defined when the others are
        # positive, which 10^x guarantees), then truncate — not round — to
        # int32 (report §4.4). NaNs stay NaN rather than tripping the cast.
        pro = out['pro'].copy()
        pro[np.isfinite(pro) & (pro < 0)] = 0.0
        for k, v in (('pro', pro), ('syn', out['syn']), ('apeuk', out['apeuk'])):
            t = v.copy()
            fin = np.isfinite(t)
            t[fin] = np.trunc(t[fin])
            out[k] = t
    return out
