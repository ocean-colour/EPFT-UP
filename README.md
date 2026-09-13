# EPFT-UP
Empirical Phytoplankton Functional Types with Uncertainty and Provenance

**EPFT-UP** is a Python package for deriving empirical phytoplankton functional
types (PFTs) from ocean optics and related observations, with rigorous treatment
of uncertainty and full provenance for every product. The goal is to generate the
algorithms, metrics, and diagnostics to share with the community.

## Installation

EPFT-UP targets Python ≥ 3.12 and is developed against the `ocean14` conda
environment.

```bash
git clone https://github.com/ocean-colour/EPFT-UP.git
cd EPFT-UP
pip install -r requirements.txt
pip install -e .
```

## MOANA baseline

`epft_up/moana/` is a Python reimplementation of NASA's PACE MOANA
picophytoplankton algorithm (Lange et al. 2020), retrained and validated on the
AMT cruises. It was developed in [IOPtics](https://github.com/ocean-colour/IOPtics)
and imported here from its `moana` branch @ `3aa3b6e` on 2026-09-13; it is kept
as the frozen reproduction baseline that later EPFT-UP algorithms are measured
against. Start with [`reports/MOANA_Claude_Report.md`](reports/MOANA_Claude_Report.md)
(what the algorithm is and how it performs) and
[`docs/design/moana_design.md`](docs/design/moana_design.md) (how this
implementation works). Data-dependent code reads the `$OS_COLOR` data tree and
writes derived products under `$OS_COLOR/EPFT-UP/`.

## Layout

- `epft_up/` — the Python package source (`epft_up/tests/` holds the test suite);
  `epft_up/moana/` is the MOANA baseline, `epft_up/data/moana/` its vendored
  NASA reference tables.
- `docs/design/` — design documents.
- `reports/` — reports, their figures and the scripts that generate them.
- `requests/` — open data requests to collaborators.
- `papers/` — reference PDFs (gitignored; URLs are recorded in the docs).
- `scripts/` — small utilities.
- `claude_prompts/` — prompts and task definitions that drive this work.

## Related work

- [IOPtics](https://github.com/ocean-colour/IOPtics) — testing and evaluating IOP
  (inherent optical property) algorithms.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
