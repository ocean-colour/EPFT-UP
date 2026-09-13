# MOANA reference tables (vendored from NASA OCSSW)

These two files are the **complete** set of tunable constants behind the NASA
MOANA algorithm (Multiple Ordination ANAlysis), which returns near-surface cell
abundances of *Prochlorococcus*, *Synechococcus* and autotrophic picoeukaryotes
from hyperspectral `Rrs` plus SST.

They are vendored here deliberately: neither file is published in Lange et al.
(2020) nor in the MOANA ATBD, so without them the algorithm cannot be run at all.
They are small (24 KB + 390 B) and NASA pins them to a 2023 OCSSW tag, so keeping
a copy in-repo makes our results reproducible without an 11 GB OCSSW install.

## Provenance

| | |
|---|---|
| Upstream | NASA Ocean Color Science Software (OCSSW), bundle `common` |
| Canonical path | `$OCDATAROOT/common/` |
| Retrieved | 2026-08-01, OCSSW tag **V2026.4** |
| Files last changed upstream | OCSSW tag **T2023.31** (unchanged since 2023) |
| Direct source | `https://oceandata.sci.gsfc.nasa.gov/manifest/tags/T2023.31/common/<file>` |
| Bundle manifest | `https://oceandata.sci.gsfc.nasa.gov/manifest/tags/V2026.4/common/manifest.json` |

SHA-256, verified against NASA's bundle manifest at download time:

```
478ace8798ce525bc6a09567b3773cadd532056c13d20ce6daf0a529f46deb57  pca_picophyto.h5
cd6b5c1cca879bc97acf6c5a3d4f47bf658d5a839b99ff4b1f534eac881756d4  picophyt.json
```

Re-verify with `shasum -a 256 *` in this directory.

## `pca_picophyto.h5` — the PCA loadings

HDF5, 24,488 bytes, two datasets:

| Dataset | Shape | dtype | Meaning |
|---|---|---|---|
| `component` | (124, 45) | float32 | loading matrix `V[λ, i]`; 45 retained principal components |
| `wavelength` | (124,) | int32 | 414, 416, …, 660 nm (2 nm steps) |

`V`'s columns are orthonormal (`max |offdiag(VᵀV)| = 1.2e-7`), i.e. a properly
truncated eigenvector basis. PC1 is a smooth blue→red ramp (max +0.171 at 414 nm,
min −0.097 at 658 nm), consistent with Lange et al.'s description of PC1 as
backscatter slope plus pure-water absorption.

The 414–660 nm / 2 nm grid here is authoritative and resolves a contradiction in
the ATBD, whose abstract instead claims a 395–705 nm requirement.

## `picophyt.json` — the regression coefficients

```json
{ "npc": "17",
  "pro_coef":   "…19 values…",
  "syn_coef":   "…18 values…",
  "apeuk_coef": "…18 values…" }
```

Vectors are **zero-padded**: their lengths (19, 18, 18) match the OCSSW C
allocations `npc+2`, `npc+1`, `npc+1` with `npc = 17`, so coefficient slot *k*
maps to a fixed principal component. Only the first `npc = 17` scores are used,
even though 45 are stored and computed.

Slot → PC mapping, as consumed by `get_Cpicophyt.c`:

- `pro_coef[0]` = intercept, `pro_coef[1]` = coefficient on `log10(SST)`,
  `pro_coef[i+2]` multiplies score `U_(i+1)`
- `syn_coef[0]` = intercept, `syn_coef[i+1]` multiplies `U_(i+1)`
- `apeuk_coef[0]` = intercept, `apeuk_coef[i+1]` multiplies `U_(i+1)`

## Known inconsistency (read before using)

Decoding `picophyt.json` through the mapping above and comparing with the ATBD
v1.2 equations:

| Taxon | PCs per this file (operational) | PCs per ATBD v1.2 |
|---|---|---|
| *Prochlorococcus* | 1, 2, 6, **7** | 1, 2, 6, **17** |
| *Synechococcus* | 1, 2, 3, 5, 8, 9, 10, 11, 12, **16** | 1, 2, 3, 5, 8, 9, 10, 11, 12, **13** |
| picoeukaryotes | 2, 3, 4, 5, 10, 13, 15 | identical ✓ |

All coefficient *values* agree (to ATBD rounding). Only the PC index of the last
*Prochlorococcus* and last *Synechococcus* term differs. The ATBD's assignment
spans **14** distinct PCs across the three taxa — matching Lange et al.'s
statement that 14 PCs survived the significance cut — whereas this file spans
**15**. That is evidence the operational table may have two coefficients in the
wrong slots, i.e. the shipping PACE products may not implement the published
algorithm.

**Settled empirically 2026-08-16** (`epft_up.moana.validation.bitexact_pace`,
one PACE granule pair): the shipping products use **this file's** assignment —
our Synechococcus matches NASA's at median |Δlog₁₀| ≈ 0.006 under this mapping
versus a +0.08 offset under the ATBD's. Which assignment was *scientifically
intended* remains a question for the authors (report §11.1.1).

Our implementation therefore supports **both** mappings behind a flag, with the
operational (this file's) mapping as the default so we can reproduce NASA's
products bit-for-bit.

## Algorithm summary

Per spectrum: linearly interpolate `Rrs` onto the 124-point grid → standardize
across wavelength, `(Rrs − mean)/sd` with the **sample (N−1)** denominator → take
scores `U_i = Σ_λ V[λ,i]·Rrs'[λ]` (a plain dot product; there is **no** training
mean to subtract) → evaluate the three regressions. *Prochlorococcus* is linear in
cells mL⁻¹ and needs `log10(SST)`; the other two are `log10`-transformed and are
exponentiated afterwards.

## References

- Lange, P. K., et al. (2020). "Radiometric approach for the detection of
  picophytoplankton assemblages across oceanic fronts." *Optics Express* 28(18),
  25682–25705. https://doi.org/10.1364/OE.398127
- MOANA ATBD v1.2 (2024-12-16). https://doi.org/10.5067/0AV5267G5C77
  (local copy: `papers/moana_atbd.pdf`)
- OCSSW implementation: `ocssw_src/src/l2gen/get_Cpicophyt.c`
