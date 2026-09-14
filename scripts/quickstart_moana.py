"""README quickstart: one PACE OCI pixel through the MOANA retrieval.

Proposed for the README (prompt "Website and RTD"/2, 2026-09-14) and kept
here so the example is executable and tested, not prose. Fetches the daily
0.1° PACE Rrs granule for one day (cached under ``$OS_COLOR/PACE/`` after the
first call; needs Earthdata credentials in ``~/.netrc``), takes the cloud-free
pixel nearest a Sargasso Sea point, and retrieves the three picophytoplankton
abundances with QC flags.

Run with::

    conda run -n ocean14 python scripts/quickstart_moana.py
"""
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import numpy as np                                      # noqa: E402
import xarray as xr                                     # noqa: E402

from epft_up.moana import run_moana                     # noqa: E402
from epft_up.moana.validation import fetch_pace_pair    # noqa: E402


def main():
    """Retrieve Pro / Syn / picoeukaryote abundances for one PACE pixel."""
    aop_path, _ = fetch_pace_pair(date='2025-07-01')    # cached after first run
    with xr.open_dataset(aop_path) as aop:
        # A daily composite is ~89 % cloud/fill, so take the cloud-free pixel
        # nearest a Sargasso Sea point (30°N, 60°W) rather than the point itself.
        rrs_all = aop['Rrs'].values                      # (lat, lon, wavelength)
        wave = aop['wavelength'].values
        lat, lon = aop['lat'].values, aop['lon'].values
        iy, ix = np.nonzero(np.isfinite(rrs_all).all(axis=-1))
        j = np.argmin((lat[iy] - 30.0) ** 2 + (lon[ix] + 60.0) ** 2)
        rrs = rrs_all[iy[j], ix[j]]
        print(f"pixel  {lat[iy[j]]:.2f}°N {-lon[ix[j]]:.2f}°W")
    out = run_moana(wave, rrs[None, :], sst=np.array([27.0]))  # SST [°C] for Pro
    for taxon in ('pro', 'syn', 'apeuk'):
        print(f"{taxon:6s} {out[taxon][0]:10.0f} cells/mL")
    print(f"flags  {int(out['flags'][0])}   (0 = clean retrieval)")


if __name__ == "__main__":
    main()
