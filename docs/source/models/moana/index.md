# MOANA

A from-scratch Python reimplementation of NASA's PACE **MOANA** algorithm
(Multiple Ordination ANAlysis; Lange et al. 2020) — a fixed linear map from
hyperspectral remote-sensing reflectance (+ SST) to near-surface cell
abundances of *Prochlorococcus*, *Synechococcus* and autotrophic
picoeukaryotes — together with its training pipeline, retraining on the AMT24
cruise, and validation on held-out cruises and on the operational PACE
product. It is EPFT-UP's **reproduction baseline**: kept frozen so that later
models are measured against something that still runs and still passes its
tests.

- **Package:** {mod}`epft_up.moana` — `io`, `pipeline`, `algorithm`,
  `train`, `validation`.
- **Entry point:** {func}`epft_up.moana.algorithm.run_moana`.
- **Provenance:** developed in [IOPtics](https://github.com/ocean-colour/IOPtics)
  and imported from its `moana` branch @ `3aa3b6e` on 2026-09-13
  ([details](../../provenance.md)).

```{toctree}
:maxdepth: 1

report
design
reproduction
open_items
```
