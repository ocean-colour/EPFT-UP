"""Numbers behind the Dove & Freilich (2026) assessment (Report/Dove+2026).

Two independent checks, both from data already in the repository / cache:

1. **The product inside their box.** Dove & Freilich analyse the PACE MOANA
   product over 63–75°W, 30–40°N (April–November 2024). We hold one daily
   0.1° L4M MOANA granule (2025-07-01). Count, inside that box, the pixels
   that are land-sentinel (254), fill, clipped-to-zero *Prochlorococcus*, and
   report the medians of the three products, to compare with the ranges in
   their Figures 3–4. The date differs from theirs; the product's behaviour
   (coefficients unchanged since OCSSW T2023.31) does not.

2. **How much *Prochlorococcus* the SST term manufactures.** MOANA's Pro model
   is ``Pro = p0 + p_SST·log10(SST) + Σ p_i·U_i`` (report §4.2). For a
   cyclonic cold-core anomaly and for the April→August seasonal warming in
   their region, compute the change in retrieved Pro that comes from SST
   alone, holding the optics fixed — the confound the paper's Pro
   conclusions have to survive.

Run with::

    conda run -n ocean14 python reports/scripts/dove2026_check.py
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from epft_up.moana.io import load_luts                                  # noqa: E402
from epft_up.moana.validation import PACE_FILL, PACE_LAND, PACE_VARS    # noqa: E402

#: Dove & Freilich (2026) study region (their "Regional study setting").
BOX = {'lon': (-75.0, -63.0), 'lat': (30.0, 40.0)}
GRANULE = os.path.join(os.environ.get('OS_COLOR', ''), 'PACE', 'moana_validation',
                       'PACE_OCI.20250701.L4m.DAY.MOANA.V3_2.0p1deg.nc')


def product_in_box(path=GRANULE, box=BOX):
    """Pixel accounting of the three MOANA fields inside the study box.

    Returns
    -------
    dict — per taxon: n_box, n_land254, n_fill, n_real, n_zero, frac_zero,
    median/p10/p90 of the real (>0) retrievals.
    """
    import xarray as xr
    out = {}
    with xr.open_dataset(path) as ds:
        lat, lon = ds['lat'].values, ds['lon'].values
        sel = ((lat[:, None] >= box['lat'][0]) & (lat[:, None] <= box['lat'][1])
               & (lon[None, :] >= box['lon'][0]) & (lon[None, :] <= box['lon'][1]))
        for taxon, var in PACE_VARS.items():
            v = ds[var].values[sel].astype(np.float64)
            land = v == PACE_LAND
            fill = (v == PACE_FILL) | ~np.isfinite(v)
            real = ~land & ~fill
            zero = real & (v == 0)
            pos = v[real & (v > 0)]
            out[taxon] = {
                'n_box': int(v.size), 'n_land254': int(land.sum()),
                'n_fill': int(fill.sum()), 'n_real': int(real.sum()),
                'n_zero': int(zero.sum()),
                'frac_zero': float(zero.sum() / max(real.sum(), 1)),
                'median': float(np.median(pos)) if pos.size else float('nan'),
                'p10': float(np.percentile(pos, 10)) if pos.size else float('nan'),
                'p90': float(np.percentile(pos, 90)) if pos.size else float('nan'),
            }
    return out


def sst_term(lut=None):
    """The SST-only change in retrieved Prochlorococcus for given SST changes.

    Returns
    -------
    dict — ``p_sst`` (the coefficient), ``dpro_per_degC`` at 20/25/28 °C,
    and ΔPro for a −1 °C and −2 °C cold-core anomaly at 25 °C and for the
    seasonal 21 → 28 °C warming of the subtropical western North Atlantic.
    """
    lut = lut or load_luts()
    p_sst = float(lut['pro_coef'][1])          # slot 1 multiplies log10(SST)
    def dpro(t0, t1):
        return p_sst * (np.log10(t1) - np.log10(t0))
    return {
        'p_sst': p_sst,
        'dpro_per_degC': {T: p_sst / (T * np.log(10)) for T in (20.0, 25.0, 28.0)},
        'cold_core_minus1C_at_25C': dpro(25.0, 24.0),
        'cold_core_minus2C_at_25C': dpro(25.0, 23.0),
        'seasonal_21_to_28C': dpro(21.0, 28.0),
    }


def main():
    print('== 1. PACE MOANA (2025-07-01, 0.1° daily) inside 63–75°W, 30–40°N ==')
    if os.path.isfile(GRANULE):
        res = product_in_box()
        for t, r in res.items():
            print(f"{t:5s} box={r['n_box']:5d} land254={r['n_land254']:3d} fill={r['n_fill']:5d} "
                  f"real={r['n_real']:4d} zero={r['n_zero']:3d} ({100 * r['frac_zero']:.1f} %)  "
                  f"median={r['median']:9.0f}  p10–p90={r['p10']:8.0f}–{r['p90']:8.0f} cells/mL")
    else:
        print(f'granule not found: {GRANULE}')

    print('\n== 2. Prochlorococcus manufactured by the SST term alone ==')
    s = sst_term()
    print(f"p_SST = {s['p_sst']:.0f}  (Pro += p_SST · log10(SST))")
    for T, d in s['dpro_per_degC'].items():
        print(f"  dPro/dSST at {T:.0f} °C: {d:8.0f} cells/mL per °C")
    print(f"  cold-core −1 °C at 25 °C : {s['cold_core_minus1C_at_25C']:+9.0f} cells/mL")
    print(f"  cold-core −2 °C at 25 °C : {s['cold_core_minus2C_at_25C']:+9.0f} cells/mL")
    print(f"  seasonal 21 → 28 °C      : {s['seasonal_21_to_28C']:+9.0f} cells/mL")


if __name__ == '__main__':
    main()
