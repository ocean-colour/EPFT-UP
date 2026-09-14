# Data

Everything EPFT-UP's MOANA baseline consumes, where it comes from, under what
licence, and how to obtain it. The repository itself is data-free apart from
the two small NASA reference tables it vendors (see [Provenance](provenance.md)).

## Inputs

| input | used for | size | source | licence | how to obtain |
|---|---|---|---|---|---|
| Flow-cytometry cell abundances, **AMT23 / 24 / 25 / 28** (Tarran et al., BODC) | training truth (AMT24) and held-out truth (23/25/28) | 0.1–0.3 MB each | [AMT23](https://doi.org/10.5285/a2104adc-e990-6789-e053-6c86abc0d557) · [AMT24](https://doi.org/10.5285/a2104adc-e98f-6789-e053-6c86abc0d557) · [AMT25](https://doi.org/10.5285/a2104adc-e98e-6789-e053-6c86abc0d557) · [AMT28](https://doi.org/10.5285/a147c314-688b-55e9-e053-6c86abc0dc81) | NERC Open Data Licence | browser download from BODC (not scriptable); `fetch_data.py` verifies the files |
| In-situ hyperspectral Rrs, AMT23/25/28 (Brewin et al. 2023, BODC) | held-out validation (target ii) | 0.8 MB | [10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6](https://doi.org/10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6) | NERC Open Data Licence | browser download from BODC; verified by `fetch_data.py` |
| Underway SST (+ IOPs, pigments), AMT24 netCDF (Jordan et al. 2025) | the SST term in training (target i) | 545 MB | [10.5281/zenodo.12527954](https://doi.org/10.5281/zenodo.12527954) | CC-BY-4.0 | `fetch_data.py` downloads and checks it |
| PACE OCI daily 0.1° granules: L3M AOP Rrs + L4M MOANA product (2025-07-01) | the granule experiment (target iii-a) and five report figures | 175 MB | NASA OB.DAAC, collections `PACE_OCI_L3M_AOP` / `PACE_OCI_L4M_MOANA` v3.2 | public (Earthdata login) | `fetch_data.py` / `epft_up.moana.validation.fetch_pace_pair` |
| **AMT24 HyperSAS Level-2 radiometry** (PML: LT, Lsky, Ed, provider Rrs; 37 days, ~1 M spectra) | the training radiometry (target i) | 15 GB (4.8 GB of `.sav` actually read) | Plymouth Marine Laboratory, private communication (2026-08) | none — unpublished | **available on request**; permission to deposit a derived, screened subset is being sought from PML |
| MOANA reference tables (`pca_picophyto.h5`, `picophyt.json`) | every retrieval | 25 KB | NASA OCSSW, tag T2023.31 | public | vendored in `epft_up/data/moana/` |

Targets (ii) and (iii-a) and five of the six report figures therefore
reproduce from public sources today. Target (i) — retraining on AMT24 —
additionally needs the PML radiometry.

## Attribution

Using these data carries the following obligations, which EPFT-UP's own
documents honour and which any derived work must carry forward:

- BODC deposits (flow cytometry, Brewin et al. 2023) — the NERC Open Data
  Licence requires the statement **"Contains data supplied by Natural
  Environment Research Council."** and citation of each dataset's DOI.
- Jordan et al. (2025) — CC-BY-4.0: cite Jordan, T. M., et al. (2025),
  *Earth System Science Data* **17**, 493–516,
  [10.5194/essd-17-493-2025](https://doi.org/10.5194/essd-17-493-2025), and
  the Zenodo record.
- PACE OCI — acknowledge the NASA Ocean Biology Processing Group / OB.DAAC
  and cite the collection versions used (v3.2).
- The PML AMT24 radiometry — Brewin et al., private communication; any use
  requires PML's agreement.

## Layout under `$OS_COLOR`

The loaders resolve every path from one root. `scripts/fetch_data.py` creates
this layout; the AMT24 sub-tree accepts a legacy `$OS_COLOR/AMT24/` location
as well.

```
$OS_COLOR/
├── AMT/
│   ├── AMT23/AMT23_JR20131005_AFC_Dataset.csv        BODC flow cytometry
│   ├── AMT24/AMT24_JR20140922_AFC_Dataset.csv
│   │       BODC_data_AMT_modern_and_historical_optical_observations.csv   Brewin 2023 Rrs
│   │       amt24_final_with_debiased_chl.nc            Jordan 2025 (Zenodo)
│   │       Radiometry/level2/<yyyyddd>/AMT24_HSAS_<yyyy-ddd>.sav   PML Level-2 (on request)
│   ├── AMT25/AMT25_JR15001_AFC_Dataset.csv
│   └── AMT28/AMT28_JR18001_AFC_Dataset.csv
├── PACE/moana_validation/PACE_OCI.<yyyymmdd>.L3m.DAY.AOP.V3_2.0p1deg.nc
│                         PACE_OCI.<yyyymmdd>.L4m.DAY.MOANA.V3_2.0p1deg.nc
└── EPFT-UP/                                             derived products (written by EPFT-UP)
    └── moana/rederivation_<date>.json, amt24_matchups_<date>.parquet
```

## `fetch_data.py`

```bash
export OS_COLOR=/path/to/data/root
python scripts/fetch_data.py            # Jordan netCDF + PACE pair; verifies BODC files
python scripts/fetch_data.py --verify   # checksums only, no downloads
```

The script downloads what can be downloaded (Zenodo, NASA), verifies every
public file it knows against pinned SHA-256 checksums, and prints the BODC
DOIs with the expected filenames for the files that must be fetched by hand.
