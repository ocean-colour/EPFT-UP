# Installation

EPFT-UP targets **Python ≥ 3.12** and is developed in the `ocean14` conda
environment. All runtime dependencies are on PyPI.

## From source

```bash
git clone https://github.com/ocean-colour/EPFT-UP.git
cd EPFT-UP
pip install -r requirements.txt
pip install -e .
```

## Earthdata login

The PACE OCI granules used by the MOANA validation and by the
[quickstart](#quickstart) are fetched through
[`earthaccess`](https://earthaccess.readthedocs.io/) and need a free
[NASA Earthdata](https://urs.earthdata.nasa.gov/) account. Put the
credentials in `~/.netrc` (mode `600`):

```
machine urs.earthdata.nasa.gov login <username> password <password>
```

Code paths that need Earthdata skip cleanly when the file is absent.

## The data root

Data-dependent code reads one directory tree, rooted at the environment
variable `$OS_COLOR`, and writes derived products under `$OS_COLOR/EPFT-UP/`.
Set it to any directory; the [Data](../data.md) page lists every input, where
it comes from, and the layout `scripts/fetch_data.py` creates under the root.

```bash
export OS_COLOR=/path/to/your/data/root
```

## Quickstart

One PACE OCI pixel through the MOANA retrieval
([`scripts/quickstart_moana.py`](https://github.com/ocean-colour/EPFT-UP/blob/main/scripts/quickstart_moana.py)
is the executable copy and is exercised by the test suite):

```python
import numpy as np, xarray as xr
from epft_up.moana import run_moana
from epft_up.moana.validation import fetch_pace_pair

aop_path, _ = fetch_pace_pair(date='2025-07-01')       # cached after the first call
with xr.open_dataset(aop_path) as aop:                  # daily 0.1° PACE Rrs, 172 bands
    rrs_all, wave = aop['Rrs'].values, aop['wavelength'].values
    lat, lon = aop['lat'].values, aop['lon'].values
    iy, ix = np.nonzero(np.isfinite(rrs_all).all(axis=-1))   # cloud-free pixels
    j = np.argmin((lat[iy] - 30.0)**2 + (lon[ix] + 60.0)**2)   # nearest to 30°N 60°W
    rrs = rrs_all[iy[j], ix[j]]

out = run_moana(wave, rrs[None, :], sst=np.array([27.0]))   # SST [°C] — Pro needs it
for taxon in ('pro', 'syn', 'apeuk'):
    print(f"{taxon:6s} {out[taxon][0]:10.0f} cells/mL")
print("flags", int(out['flags'][0]))                        # 0 = clean retrieval
```

```
pro        150219 cells/mL
syn          3337 cells/mL
apeuk         893 cells/mL
flags 0
```

*Prochlorococcus* needs SST; here a nominal 27 °C — the validation pipeline
takes it from the cruise's underway record, and the operational product from
GHRSST. `flags` is a QC bit-word (`FLAG_NEGATIVE_PRO`, `FLAG_EXTRAPOLATED`,
`FLAG_TOO_FEW_BANDS`, …; see {mod}`epft_up.moana.algorithm`): EPFT-UP returns
raw floats plus flags rather than NASA's silent clamp-to-zero and integer
truncation.

## Running the tests

Data-independent (Tier-1) tests run anywhere; data-dependent (Tier-2) tests
skip automatically when `$OS_COLOR` or a specific dataset is absent:

```bash
pytest -q
```
