"""How the target (iii-a) granule statistics depend on the pixel subsample.

Written during the Move acceptance run (2026-09-13): the report's §7.1 box
quotes Synechococcus median Δlog₁₀ / MAD of −0.005 / 0.007 (operational) and
+0.082 / 0.090 (ATBD) "on 100k ocean pixels", while the prompt-14 log quotes
+0.084 for the ATBD median. ``bitexact_pace`` draws its subsample with
``default_rng(seed).choice``, so the printed third decimal moves with
``max_pixels``. This probe runs the experiment at several sizes so the
re-derivation can state which one the report's numbers came from.

Run with::

    conda run -n ocean14 python scripts/probe_pace_subsample.py
"""
import os
import sys

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, _REPO)
sys.path.insert(0, os.path.join(_REPO, 'reports', 'scripts'))
from epft_up.moana.validation import bitexact_pace          # noqa: E402
import moana_report_figs as figs                             # noqa: E402


def figure_method_stats(n_pixels=60_000, seed=0):
    """Replicate plot_mapping_verdict's statistics without drawing.

    The figure differs from bitexact_pace in three ways: its own subsample
    (``nasa > 0`` mask, 60k pixels), no ``nasa_compat`` truncation, and
    Synechococcus only. Returns {mapping: (median dlog, MAD)}.
    """
    import numpy as np
    import xarray as xr
    from epft_up.moana import run_moana
    aop = xr.open_dataset(figs.AOP_GRANULE)
    moana = xr.open_dataset(figs.GRANULE)
    rrs = aop["Rrs"].sel(lat=moana["lat"].values, lon=moana["lon"].values,
                         method="nearest", tolerance=1e-3)
    wave = aop["wavelength"].values.astype(float)
    nasa = moana["syncoccus_moana"].values
    valid = np.isfinite(nasa) & (nasa != figs.LAND) & (nasa != figs.FILL) & (nasa > 0)
    iy, ix = np.nonzero(valid)
    pick = np.random.default_rng(seed).choice(iy.size, min(n_pixels, iy.size),
                                              replace=False)
    iy, ix = iy[pick], ix[pick]
    spectra = rrs.values[iy, ix, :]
    theirs = np.log10(nasa[iy, ix].astype(float))
    aop.close(), moana.close()
    out = {}
    for mapping in ("operational", "atbd"):
        res = run_moana(wave, spectra, sst=None, pc_mapping=mapping)
        ok = np.isfinite(res["syn"]) & (res["syn"] > 0)
        dlog = np.log10(res["syn"][ok]) - theirs[ok]
        out[mapping] = (float(np.median(dlog)), float(np.median(np.abs(dlog))))
    return out


def main():
    """Print syn median/MAD Δlog₁₀ under both mappings for several subsamples."""
    fm = figure_method_stats()
    print("figure method (plot_mapping_verdict, 60k, nasa>0, no nasa_compat): "
          + "   ".join(f"{m} {v[0]:+.4f} (MAD {v[1]:.4f})" for m, v in fm.items()))
    for n in (60_000, 100_000, 200_000, 500_000):
        res = bitexact_pace(max_pixels=n, seed=0, verbose=False)
        op, at = res['operational']['syn'], res['atbd']['syn']
        print(f"max_pixels={n:>7,d}  n={op['n']:>7,d}  "
              f"operational {op['median_dlog']:+.4f} (MAD {op['mad_dlog']:.4f}, "
              f"exact {op['frac_exact']:.3f})   "
              f"atbd {at['median_dlog']:+.4f} (MAD {at['mad_dlog']:.4f})   "
              f"verdict {res['verdict']}")


if __name__ == "__main__":
    main()
