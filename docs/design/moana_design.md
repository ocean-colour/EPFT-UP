# MOANA — design document

> **Provenance.** Imported from IOPtics `moana` @ `3aa3b6e` on 2026-09-13; original path `docs/design/moana_design.md`.

**Author:** Claude (Fable 5), for J. Xavier Prochaska
**Date:** 2026-08-16 (rev. 2 — Q&A #33–35 answered and folded in; no open
decisions remain on the matchup definition)
**Status:** design; no `epft_up/moana/` code exists yet. Prompt 12 implements this.
**Companion:** `reports/MOANA_Claude_Report.md` (the *what and why* of the
algorithm); this document is the *how* of our implementation. Decisions cited as
"Q&A #n" live in `claude_prompts/moana_IOPtics_history.md`.

---

## 1. Purpose and scope

We are reimplementing NASA's MOANA algorithm (Lange et al. 2020; operational on
PACE OCI) in Python, in order to:

1. **(i)** retrain it from scratch on the AMT24 cruise data and reproduce
   Lange et al. Tables 1–2;
2. **(ii)** apply the published model to genuinely held-out cruises
   (AMT23/25/28, in-situ hyperspectral Rrs);
3. **(iii)** validate the operational PACE product against in-situ cell counts —
   which has never been done — including a bit-exactness reproduction of one
   PACE granule that settles the PC-mapping discrepancy (§4.1).

MOANA is a **separate track** from the IOP work (Q&A #3): no `AlgorithmSpec`, no
BING wrapper, no shared metrics table. Code lives in `epft_up/moana/`, reference
data in `epft_up/data/moana/`, tests in `epft_up/tests/test_moana.py`.

Two components make up the work, and this document treats them separately
because their inputs and failure modes are different:

- **The retrieval** (§3) — a fixed linear map from `Rrs(λ)` (+ SST) to three
  cell abundances. Fully specified by two vendored NASA files. Deterministic,
  testable to machine precision.
- **The training pipeline** (§5–§6) — from the raw-ish AMT24 HyperSAS Level-2
  stream to a training matrix, and from there to our own PCA + regressions.
  Full of judgment calls; every one of them is either pinned to Lange et al.'s
  published thresholds or recorded as an explicit, configurable decision.

## 2. Sources of truth

| Source | Role | Location |
|---|---|---|
| `pca_picophyto.h5` | PCA loadings `V[124, 45]` + wavelength grid — **the** algorithm constants | `epft_up/data/moana/` (vendored, sha256-pinned) |
| `picophyt.json` | 25 regression coefficients, zero-padded, `npc=17` | `epft_up/data/moana/` |
| OCSSW `get_Cpicophyt.c` | authoritative operational behaviour (interp, N−1, clamp, cast) | read 2026-08-01; not vendored |
| Lange et al. 2020 | science, training design, processing thresholds, skill tables | `papers/lange2020.pdf` (NTRS accepted manuscript¹) |
| MOANA ATBD v1.2 | operational recipe, coefficient table, PC assignment | `papers/moana_atbd.pdf` |
| AMT24 HyperSAS Level-2 | training radiometry (LT, LI/Lsky, ES/Ed + ancillaries) | `$OS_COLOR/AMT24/Radiometry/level2/` — **use the `.sav` files only** (§5.1) |
| AMT24 flow cytometry | training truth (cells mL⁻¹) | `$OS_COLOR/AMT24/AMT24_JR20140922_AFC_Dataset.csv` (BODC; Tarran & Zubkov 2020) |
| Jordan et al. 2025 netCDF | `uway_sst` (1-min, QC'd) + IOPs/pigments | `$OS_COLOR/AMT24/amt24_final_with_debiased_chl.nc` |

¹ *The `papers/` PDFs were lost in the machine migration (`papers/*.pdf` is
gitignored) and restored 2026-08-16. `lange2020.pdf` is now the NTRS accepted
manuscript (24 pp, different pagination from the Optics Express version but
identical content and equation numbering); the ATBD was re-fetched from
oceancolor.gsfc.nasa.gov.*

## 3. The retrieval, as we will implement it

Per spectrum (see the report §4 for derivations and invariants):

1. **Band screen.** Drop non-finite bands. Require a configurable minimum number
   of valid bands (default: all 124 target bands bracketed) — the OCSSW code has
   no such guard and silently extrapolates; we will not repeat that.
2. **Interpolate** linearly onto `wavelength` from the LUT: 414, 416, …, 660 nm
   (124 bands). Out-of-range behaviour is a flag, not silent clamping.
3. **Standardise:** `Rrs' = (Rrs − mean(Rrs)) / sd(Rrs)` across the 124 values,
   **sample (N−1) sd** (Q&A #8). No training mean exists to subtract.
4. **Project:** `U_i = Σ_λ V[λ,i]·Rrs'(λ)`, all 45 stored components (only the
   first 17 are consumed, but the extras are free and §9.3 of the report wants
   the reconstruction residual eventually).
5. **Evaluate:**
   - `Pro = p₀ + p_SST·log₁₀(SST) + Σ pᵢ·Uᵢ` — **linear** in cells mL⁻¹;
   - `log₁₀(Syn) = s₀ + Σ sᵢ·Uᵢ`, then `10^x`;
   - `log₁₀(Apeuk) = a₀ + Σ aᵢ·Uᵢ`, then `10^x`.

**Output policy (Q&A #11):** return raw floats plus a QC bitmask
(`NEGATIVE_PRO`, `EXTRAPOLATED_BANDS`, `TOO_FEW_BANDS`, …). A separate
`nasa_compat=True` mode reproduces OCSSW exactly — clamp negative Pro to 0,
truncate all three to `int32` — for the bit-exactness test and nothing else.

**Vectorisation:** the natural unit is a 2-D array `(n_spectra, n_bands)`; all
steps above are `numpy` matrix operations. No per-pixel Python loops.

### 3.1 Invariants (free tests)

- `V`'s columns are orthonormal to float32 precision.
- After step 3: `mean(Rrs') = 0`, `‖Rrs'‖₂ = √123 ≈ 11.0905` — identically,
  every spectrum.
- Therefore `|U_i| ≤ 11.0905` and `Σ U_i² ≤ 123`. Violation of any of these is
  an implementation bug, and `epft_up/tests/test_moana.py` asserts all three.

## 4. Discrepancies the code must carry

### 4.1 PC mapping: operational vs ATBD (report §7.1)

`pc_mapping="operational"` (default) uses `picophyt.json` slots as-is
(Pro on PC 1,2,6,**7**; Syn on …,**16**). `pc_mapping="atbd"` moves the two
disputed coefficients to U**17** and U**13**. Default is operational so the
PACE bit-exactness test can pass; the flag settles which one NASA actually
runs (Q&A #9). Both mappings share the same coefficient *values*.

### 4.2 Spectral range

414–660 nm @ 2 nm is authoritative (the LUT itself); the ATBD's 395–705 nm
figure is treated as a data-availability screen only (Q&A #7).

### 4.3 NEW — the paper gives picoeukaryotes an SST term; the ATBD and the LUT do not

Found 2026-08-16 while pinning the training equations. Lange et al.'s **Eq. 7**
is `log₁₀(y_Apeuk) = a + b₀·log₁₀(SST) + Σ bᵢuᵢ`, and the accompanying prose
says *"In the final regressions, SST was used as an additional predictor for
Prochlorococcus **and picoeukaryotes**"*. But the ATBD's coefficient table and
`picophyt.json` (C allocation `npc+1`, no SST slot) have **no SST term for
picoeukaryotes**. Since the paper tabulates no coefficient values, we cannot
tell whether the operational model dropped the term or the paper's Eq. 7/prose
is wrong. Consequences:

- `train.py` takes a per-taxon `use_sst` switch, so we can fit both variants
  and see which reproduces the published skill.
- This joins the §11 NASA-contact list in the report (to be added there when
  the report is next revised).

### 4.4 Training-side subtlety: `prcomp` centring

Lange used R's `prcomp`, whose **default centres each column** (per-wavelength
mean over samples) before the SVD. But the operational retrieval projects
`Rrs'` straight onto `V` with no training mean subtracted — mathematically
consistent only if training used `center=FALSE` (or the column means were
negligible). We cannot know which from the documents. `train.py` therefore
supports both (`center=True/False`), and the reproduction experiment compares
each basis against the vendored LUT (§6.3); whichever matches is what Lange
did. This is exactly the kind of silent convention that makes "reproductions"
fail, so it is a first-class experiment, not a footnote.

## 5. The AMT24 Level-2 → training-matrix pipeline

### 5.1 Stage 0 — ingest (`io.py`)

Read the per-day IDL saves
`$OS_COLOR/AMT24/Radiometry/level2/<yyyyddd>/AMT24_HSAS_<yyyy-ddd>.sav` via
`scipy.io.readsav`: arrays `matrix_lt`, `matrix_li`, `matrix_es`, `matrix_rrs`,
each `(n_spectra, 158)`. **Never read the CSVs** — `*_ES.dat` is a
byte-identical copy of `*_LT.dat` on all 37 days (Q&A #28).

Column layout (all four matrices identical): columns 0–13 ancillary
(ship start time, time [decimal hr UTC], lat, lon, pressure, PAR, wind dir,
wind speed [m s⁻¹], ship orientation, salinity, SST, **Δ-azimuth [deg]**,
azimuth, solar zenith), columns 14–154 the spectrum on 306.0–796.0 nm at
exactly 3.5 nm (141 bands), columns 155–157 THS compass, pitch [deg],
roll [deg]. The reader returns a tidy structure (per-day dict or `xarray`)
with the wavelength axis explicit, so no downstream code ever indexes raw
columns.

Known data facts the reader must respect: DOY 267, 273, 298 are absent
(273 and 298 are inside Lange's 30 Sep–1 Nov window — PML follow-up, Q&A #30);
DOY 266's ancillary SST/salinity are garbage (pre-departure transit) but its
radiometry is well-formed; everything is daylight-only at ~0.86 s cadence,
993,417 spectra total.

### 5.2 Stage 1 — geometry screening (`pipeline.py`)

Exact Lange et al. thresholds (their §2.3, our verbatim extraction 2026-08-16):

| screen | keep | source |
|---|---|---|
| tilt (from pitch & roll) | `tilt < 5°` | Lange step 3 |
| solar zenith | `10° < θ₀ < 80°` | Lange step 3 |
| relative azimuth Δϕ | `50° < Δϕ < 170°` | Lange step 4 |

Tilt is computed from the pitch/roll columns
(`tilt = acos(cos(pitch)·cos(roll))`, the standard small-angle-exact form).
Δϕ: the files carry a signed `DELTA-AZIMITH-ANGLE` spanning −165°…+64°; Lange's
convention is `Δϕ = ϕ_sensor − ϕ_sun`. We screen on `|Δϕ|` after verifying the
file's convention against the azimuth and heading columns on a few clear days —
recorded as a pipeline unit check, not assumed.

### 5.3 Stage 2 — glint minimisation and Rrs (`pipeline.py`)

Lange's two-part scheme, applied to the screened stream:

1. **1-minute selection:** partition each day into 1-minute intervals; within
   each interval keep only the spectrum with the minimum mean `Lt(750–800 nm)`.
   This collapses ~70 samples/min to 1 and is itself the first glint filter.
2. **Per-spectrum correction:** solve
   `min_{ρ_sky, L_NIR} Σ_{λ=750}^{800} | Lt(λ) − ρ_sky·Lsky(λ) − L_NIR |`
   (an L1 fit — two parameters, ~14 bands at 3.5 nm; small fixed iteration or
   `scipy.optimize`), then
   `Lw(λ) = Lt(λ) − ρ_sky·Lsky(λ) − L_NIR` and `Rrs(λ) = Lw(λ)/Ed(λ)`.

Bounds: `ρ_sky ∈ [0, 0.1]`, `L_NIR ≥ 0`, both recorded per spectrum as
diagnostics (their distributions are themselves a QC product — a ρ_sky pile-up
at a bound flags bad intervals).

The provider's own `matrix_rrs` (fixed ρ = 0.0280, no NIR offset, no screening —
established empirically, prompt-8 report) is retained as a **cross-check
channel** only (Q&A #29): the pipeline can emit both and difference them.

### 5.4 Stage 3 — spectral QC (`pipeline.py`)

Lange's three post-processing screens, applied to the 1-min Rrs:

1. **Local time in [09:00, 17:00].** "Local" is not defined in the paper; we
   use **solar time** (`UTC + lon/15°` hours), which is the only definition
   that makes sense on a meridional transect, and record it as an assumption.
2. **No negative Rrs in 400–700 nm** (evaluated on the native 3.5 nm grid).
3. **Second-derivative noise filter:** drop spectra with
   `|d²Rrs/dλ²| > 2×10⁻⁴ sr⁻¹ nm⁻²` anywhere in 610–660 nm, central
   differences on the grid actually in use. (The paper prints the threshold's
   unit as sr⁻¹ nm⁻¹; a second derivative is per-nm², and the code comments say
   so. The filter is applied after resampling to 2 nm — same as Lange, whose
   order was interpolate-then-QC.)

### 5.5 Stage 4 — resample to the MOANA grid

Linear interpolation from the 141-band 3.5 nm grid to 414–660 @ 2 nm
(124 bands). We accept the double-interpolation smoothing for now (Q&A #31).
The resampler is a shared utility with the retrieval's step 2 — one
implementation, one set of tests.

### 5.6 Stage 5 — matchup and binning

Truth table: `AMT24_JR20140922_AFC_Dataset.csv` — 814 CTD-bottle samples,
68 stations, 1–31 Oct 2014; taxa arrive as BODC parameter codes
(`P700A90Z`, `P701A90Z`, `PYEUA00A`, …) whose mapping to
Pro/Syn/picoeukaryotes is resolved from the deposit's own metadata documents
(in `$OS_COLOR/AMT24/`) inside `io.py`, with the mapping asserted in tests —
not hard-coded silently.

- **Depth cut (Q&A #33, confirmed):** `depth ≤ 10 m`, shallowest bottle per
  station; both window and rule remain parameters, but these are the defaults
  and the values used for every reported result.
- **Binning (Q&A #32):** for each FCM sample, take the **median** of screened
  Rrs spectra within **±15 min**, per band; keep the 16–84 % spread as a
  per-band uncertainty and the n-in-bin as a weight/QC field. Window and
  estimator are parameters.
- **Lange-strict mode:** Lange matched by *exact date-hour-minute* to the
  1-min-selected spectrum. For the Tables 1–2 reproduction we support
  `window=0` (nearest 1-min spectrum) alongside our robust default, so the
  reproduction is apples-to-apples while our own training uses the sturdier
  estimator.
- **SST (Q&A #34, confirmed):** primary = `uway_sst` from the Jordan netCDF,
  linearly interpolated to the bin's mid-time; the HSAS ancillary SST column is
  a cross-check only (it is corrupt early in the cruise). The pipeline computes
  and stores both, plus their difference as a QC diagnostic.

**Known coverage gap:** the BODC deposit is CTD-bottle-only, while Lange's
training set (n = 73–78) also used the ~30-minute **underway** FCM samples —
the very samples the paper shows are what make *Synechococcus* retrievable
(CTD-only retraining degraded Syn MAE 1.27→1.37). With the confirmed ≤10 m /
shallowest-bottle rule we have 68 training samples — the *sparse-CTD*
configuration. Per Q&A #35 this is now item 1 on the PML follow-up list
(`requests/PML_follow_up.md`); until it is resolved, our Tables 1–2
reproduction has a known, quantified handicap — and Lange's own CTD-only rows
(their Table 2 sensitivity) are the fair comparison line.

### 5.7 Pipeline products

The pipeline emits one netCDF/parquet per run under
`$OS_COLOR/AMT24/derived/` (the repo stays data-free): the screened 1-min Rrs
stream (with per-stage QC flags, ρ_sky, L_NIR), and the matched training matrix
(spectra × 124, plus FCM counts, SST, uncertainties, station metadata,
provenance). Every stage records counts in/out, so the attrition table
(raw → tilt → geometry → glint-selected → QC → matched) is a standard output —
that table is the first thing to compare against Lange's n.

## 6. Retraining (`train.py`)

### 6.1 PCA

On the matched, standardised training matrix `X (n × 124)`:
SVD-based PCA with `center` configurable (§4.4), no scaling. Rank is
`min(n−1, 124)` — with n ≈ 70–130 we cannot get 45 meaningful components unless
n permits; we keep what exists. Then Lange's two cuts, both parameterised:
discard PCs with `sd < 0.1 %` of PC1's sd (paper: left 20); the regression
stage discards more by significance (paper: left 14).

### 6.2 Regressions

Per taxon: multilinear regression of the target (`Pro` linear;
`log₁₀(Syn)`, `log₁₀(Apeuk)`) on the retained PC scores, with per-taxon
`use_sst` (§4.3). Predictor selection follows Lange: start with all candidates,
**backward stepwise elimination of the highest-p predictor, accepting each
removal only if AIC decreases**, iterated to the paper's stopping rule.
Stepwise selection is statistically shaky at this n (report §9.5) — but the
goal here is *reproduction*, so we implement their procedure exactly and keep
ridge/CV variants out of scope for this track.

### 6.3 Comparing our basis to NASA's

Eigenvector sign and order are arbitrary, so the comparison is a
**cosine-similarity matrix** between our loadings and the LUT's 45 columns
(after aligning grids): `|cos| > 0.99` on the leading components under one of
the two centring conventions is the success criterion, with a
Procrustes-style summary for the rest. Skill reproduction targets Lange
Table 1 (hyperspectral, in-situ): bias/MAE/R² within the tolerance the n
mismatch (§5.6) allows, alongside the CTD-only comparison line.

## 7. Code layout and conventions

```
epft_up/moana/
    __init__.py      # public surface: run_moana(), load_luts(), …
    io.py            # .sav reader, LUT loader, FCM/BODC reader, Jordan-SST reader
    pipeline.py      # stages 1–5: screening, glint, QC, resampling, matchup
    algorithm.py     # the retrieval (§3): standardise, project, evaluate, QC flags
    train.py         # PCA + stepwise regressions + basis comparison (§6)
    validation.py    # prompt 13: targets (i)–(iii)  [stub until prompt 13]
epft_up/data/moana/  # vendored LUTs (already in place, sha256-pinned)
epft_up/tests/test_moana.py
```

Per the project conventions and Q&A decisions: **methods, not classes** — the
natural state objects are plain dicts/`xarray.Dataset`s passed explicitly;
imports at top; every method's docstring lists inputs/outputs; inline comments
explain the *why* (especially every Lange threshold, each of which cites the
paper section). All thresholds live in one module-level config dict
(`DEFAULT_PIPELINE`) so a run's provenance is a serialisable dict, matching the
IOPtics provenance habit. Plotting (matplotlib, `Agg`) goes in the report
scripts, not the package.

Dependencies: `numpy`, `scipy` (`io.readsav`, `optimize`), `pandas`, `h5py`,
`xarray` — all already in `requirements.txt`; `earthaccess` only inside
validation target (iii) code paths (skip-guarded, Q&A #16).

## 8. Testing strategy (prompt 12 preview)

- **Pure-math tests, no data needed:** the §3.1 invariants; a synthetic
  spectrum through standardise→project→evaluate against hand-computed values;
  both PC mappings produce identical results when the disputed coefficients are
  zeroed; `nasa_compat` clamping/truncation.
- **LUT tests:** sha256 of the vendored files; orthonormality; grid = 414:2:660.
- **Pipeline tests on synthetic streams:** a fabricated day with known glint
  (`ρ_sky`, `L_NIR` planted) is recovered; screening thresholds cut exactly at
  5°/80°/10°/50°/170°; the second-derivative filter fires on a planted spike.
- **`$OS_COLOR`-gated tests:** one real `.sav` day loads; ES ≠ LT (guards the
  known CSV bug never leaking in); column-layout assertions; FCM code mapping.
- **`~/.netrc`-gated tests:** none until validation (iii).

## 9. Open items

| # | Item | State |
|---|---|---|
| 1 | Q&A #33 (FCM depth cut) and #34 (SST source) | **closed** — confirmed as proposed (≤10 m shallowest; Jordan `uway_sst` primary); folded into §5.6 |
| 2 | Underway FCM samples missing from the BODC deposit | Q&A #35 — on the PML follow-up list (`requests/PML_follow_up.md`, item 1); until then Tables 1–2 reproduction runs in the CTD-only configuration |
| 3 | Paper-vs-ATBD SST term for picoeukaryotes (§4.3) | handled by per-taxon `use_sst`; add to the report's §11 NASA list at next revision |
| 4 | `prcomp` centring convention (§4.4) | settled empirically by §6.3 |
| 5 | DOY 273/298, ES-bug report, processing note, native-grid spectra | consolidated with everything else PML-bound in `requests/PML_follow_up.md` |
| 6 | "Local time" definition in Lange QC | assumed solar time; documented in code |
