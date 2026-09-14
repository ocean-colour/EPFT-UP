# Provenance and citing

## Where the code came from

`epft_up/moana/` was developed in [IOPtics](https://github.com/ocean-colour/IOPtics)
between 2026-08-01 and 2026-09-10 and imported into this repository from the
IOPtics `moana` branch at commit `3aa3b6e` on 2026-09-13, as a fresh copy
(no history rewrite). Each migrated document carries a provenance header
naming that commit and its original path; the import changed nothing but the
`ioptics` → `epft_up` package prefix, and the full re-derivation on
[the reproduction page](models/moana/reproduction.md) is the evidence that
nothing else changed. The complete working record of the IOPtics phase —
20 prompts, 36 Q&A decisions, and their logs — is preserved verbatim in the
repository at
[`claude_prompts/moana_IOPtics_history.md`](https://github.com/ocean-colour/EPFT-UP/blob/main/claude_prompts/moana_IOPtics_history.md).

## Vendored NASA reference tables

```{include} ../../epft_up/data/moana/README.md
:heading-offset: 2
```

## Citing EPFT-UP

The repository carries a
[`CITATION.cff`](https://github.com/ocean-colour/EPFT-UP/blob/main/CITATION.cff)
(GitHub's *Cite this repository* button reads it). Until the first release is
archived on Zenodo, cite the repository:

> Prochaska, J. X. and Claude (2026). *EPFT-UP: Empirical Phytoplankton
> Functional Types with Uncertainty and Provenance* (v0.1.0) [software].
> https://github.com/ocean-colour/EPFT-UP

A Zenodo DOI will be minted at the `v0.1.0` release and added here and to the
README badge row.

## Data attribution

See the [Data](data.md) page for the licence obligations that travel with the
inputs (NERC Open Data Licence statement; CC-BY-4.0 for Jordan et al. 2025;
NASA OB.DAAC acknowledgement).
