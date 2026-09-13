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

## Layout

- `epft_up/` — the Python package source (`epft_up/tests/` holds the test suite).
- `claude_prompts/` — prompts and task definitions that drive this work.

## Related work

- [IOPtics](https://github.com/ocean-colour/IOPtics) — testing and evaluating IOP
  (inherent optical property) algorithms.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
