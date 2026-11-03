Written by Opus 5.5 and Billy Li.

# Claude Fun — design document

## 1. Purpose and scope

This design document describes Claude Fun, a sub-project of EPFT-UP where we reproduce MOANA's (Lange et al. 2020) capabilities and improve upon it purely by incorporating more data and using a bigger model. We do not attempt to incorporate more physics/biology for this model. The objective is to probe the amount of skill a machine learning model can pick up from Rrs spectra (plus SST, as in MOANA) paired with cell counts alone, given as much currently available data as we can, and without additional ocean modelling, nutrient constraints, or radiative transfer modelling. This document describes the project to Claude Fable 5.1, who will write the code.

To allow for sufficient skill and good scientific usage, we plan to use a Bayesian fully connected neural network for the model. With MOANA, the learning is heavily limited to correlations between large-scale spectral features and water mass/location/time, and thus cell counts. For Claude Fun, the network should be large enough to learn patterns beyond correlations between broad spectral shape and cell counts, including possible spectral features that represent the cells themselves.

The model output is a predicted mean and uncertainty for the cell abundances of *Prochlorococcus*, *Synechococcus*, and picoeukaryotes, jointly Gaussian in log10(cells mL⁻¹) with a 3×3 covariance, so correlations between taxa are kept (§6.2). Every marginal of a multivariate Gaussian is itself Gaussian, so the reported per-taxon mean and σ are exactly what a diagonal model would report, and the comparison with MOANA is unchanged. Unlike MOANA, *Prochlorococcus* is also represented in log space. Because the output is a normal distribution, the model is trained with a Gaussian negative log-likelihood. The Bayesian treatment adds a KL-divergence term between the posterior over network weights and a weight prior (§6.3). No prior distribution on the cell counts themselves is used: the aim is to measure what the model learns from data alone.

The model input is designed to accept as many types of data as possible without losing input information. MOANA uses only the AMT cruises, which are spatially and temporally limited and represent only a few populations. We intend to at least incorporate SeaFlow+PACE, with SeaFlow data from the Ribalet lab (UW) matched to PACE observations. MOANA degrades the input spectra to a fixed grid, whereas we intend the model to read spectra of any resolution. The input is therefore a set of per-band tokens (λ, Rrs, σ_Rrs), where σ_Rrs is the Rrs uncertainty, read by an attention encoder (§5–6). Instrument bandwidth (FWHM) is deliberately not an input: with only a few instruments it would act as an instrument label and could add bias and spurious correlation rather than information.

For the current stage, we train this model on the AMT24 cruise data and compare it to MOANA. We then apply the trained model to AMT23/25/28, and finally test it on SeaFlow+PACE pairs. Spectral attribution maps show which wavelengths drive each prediction (§8.7). In the following stage (next design document, outline in §14), we will train on both AMT and SeaFlow+PACE, and develop a data processor that ingests the different input data forms.

### 1.1 What this stage can and cannot show

AMT24 gives few independent labels. The confirmed depth rule yields one sample per CTD station (68 stations in the deposit), but the HyperSAS stream is daylight-only and Stage 3 keeps 09:00–17:00 solar time. From the AFC file itself (Stage 5): all 68 casts have a ≤10 m bottle, but only 35 surface samples fall in daylight (solar 09:00–17:00), and 2 of those are on days with no radiometry (DOY 273, 298). **At most 33 labels can be matched to radiometry** (28 within Lange's 30 Sep–1 Nov window); the Stage 1–3 screens can only lower this, and the attrition table gives the final number. That is enough to validate the pipeline, the architecture, uncertainty calibration and resolution robustness, and to compare against MOANA on equal footing. It is not enough to support a large network or to demonstrate learning of fine spectral features; that test belongs to the next stage. This stage therefore uses deliberately small defaults (§6) and treats "bigger model" as a next-stage question.

## 2. Sources of truth

All AMT paths are relative to one config key, `AMT_ROOT = $OS_COLOR/AMT` (confirmed: every cruise sits under the `AMT/` parent folder).

| source | path / origin | used for |
|---|---|---|
| AMT24 HyperSAS Level-2 (PML) | `<AMT_ROOT>/AMT24/Radiometry/level2/<yyyyddd>/AMT24_HSAS_<yyyy-ddd>.sav` — `.sav` files only | training inputs |
| AMT24 flow cytometry (BODC; Tarran & Zubkov 2020) | `<AMT_ROOT>/AMT24/AMT24_JR20140922_AFC_Dataset.csv` | training labels |
| AMT23 / 25 / 28 flow cytometry (BODC) | `<AMT_ROOT>/AMT23/AMT23_JR20131005_AFC_Dataset.csv`, `<AMT_ROOT>/AMT25/AMT25_JR15001_AFC_Dataset.csv`, `<AMT_ROOT>/AMT28/AMT28_JR18001_AFC_Dataset.csv` | transfer-test labels (§8.5) |
| Brewin et al. 2023 AMT optical compilation (BODC) | `<AMT_ROOT>/AMT24/BODC_data_AMT_modern_and_historical_optical_observations.csv` | transfer-test Rrs, σ_Rrs and SST for AMT23/25/28 (§8.5); calibration of the σ floor (§5.2). Covers AMT23/25/26/28 only — no AMT24 rows |
| Jordan et al. 2025 netCDF (Zenodo) | `<AMT_ROOT>/AMT24/amt24_final_with_debiased_chl.nc` | primary SST, `uway_sst` (Stage 5) |
| BODC series documentation for the AFC deposits | BODC metadata reports (e.g. AMT24: series 2009640) | parameter-code → taxon mapping, qualifier flags (Stage 5) |
| Lin et al. 2022, *Opt. Express* 30, 45648 | paper | context for the σ floor (§5.2) |
| Lange et al. 2020, Table 1 (and Table 2 CTD-only rows) | transcribed to `epft_up/claude_fun/reference/lange2020_tables.csv` | published reference numbers (§8.3) |
| `epft_up.moana` | in-repo | MOANA retrained on identical folds; published MOANA on held-out cruises (§8.3) |
| SeaFlow v2 cell counts | Ribalet lab, UW server | test only this stage (§8.6) |
| PACE OCI L2 Rrs | NASA, via `earthaccess` | test only this stage (§8.6) |

The vendored MOANA LUTs (`epft_up/data/moana/`) are used only through `epft_up.moana` for the baseline, not by Claude Fun's model.

## 3. The AMT24 Level-2 → training-matrix pipeline

Copied from `EPFT-UP/docs/design/moana_design.md`, with these **Claude Fun changes**:

1. **Stage 3 noise filter** is evaluated on a QC-only 2 nm resample (identical to Lange); only the keep/drop flag is carried forward.
2. **Stage 4 resampling is QC/baseline-only.** The model consumes spectra on the native 3.5 nm grid.
3. **Stage 5** stores a per-band σ_Rrs, keeps the member 1-min spectra of each bin (for replicate augmentation, §7.2), reads BODC qualifier flags, and handles zero / below-detection counts.
4. The training matrix is spectra × 141 native bands (the model then selects its wavelength range, §5.1), not × 124.

### Stage 0 — ingest (`io.py`)

Read the per-day IDL saves
`<AMT_ROOT>/AMT24/Radiometry/level2/<yyyyddd>/AMT24_HSAS_<yyyy-ddd>.sav` via
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

### Stage 1 — geometry screening (`pipeline.py`)

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

### Stage 2 — glint minimisation and Rrs (`pipeline.py`)

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

### Stage 3 — spectral QC (`pipeline.py`)

Lange's three post-processing screens, applied to the 1-min Rrs:

1. **Local time in [09:00, 17:00].** "Local" is not defined in the paper; we
   use **solar time** (`UTC + lon/15°` hours), which is the only definition
   that makes sense on a meridional transect, and record it as an assumption.
2. **No negative Rrs in 400–700 nm** (evaluated on the native 3.5 nm grid).
3. **Second-derivative noise filter:** drop spectra with
   `|d²Rrs/dλ²| > 2×10⁻⁴ sr⁻¹ nm⁻²` anywhere in 610–660 nm, central
   differences on the grid actually in use. (The paper prints the threshold's
   unit as sr⁻¹ nm⁻¹; a second derivative is per-nm², and the code comments say
   so.) *Claude Fun change:* the threshold was set by Lange on a 2 nm grid, and a
   central-difference second derivative is not grid-invariant (its noise scales
   as 1/h²), so the filter is evaluated on a **temporary 2 nm linear resample**
   (the Stage 4 utility) — interpolate-then-QC, exactly as Lange — and only the
   resulting keep/drop flag is carried forward. The model never sees the
   resampled spectrum.

### Stage 4 — resample to the MOANA grid (QC and baseline only)

Linear interpolation from the 141-band 3.5 nm grid to 414–660 @ 2 nm
(124 bands). In Claude Fun this is used only (a) inside the Stage 3 noise
filter, (b) by the MOANA baseline (§8.3), and (c) by the resolution-robustness
test (§8.4). Claude Fun's model input stays on the native grid. The resampler is
a shared utility with `epft_up.moana` — one implementation, one set of tests.

### Stage 5 — matchup and binning

Truth table: `AMT24_JR20140922_AFC_Dataset.csv` — 814 CTD-bottle samples,
68 casts, 25 Sep–1 Nov 2014; taxa arrive as BODC parameter codes whose mapping to
Pro/Syn/picoeukaryotes is resolved from the deposit's own metadata inside
`io.py`, with the mapping asserted in tests — not hard-coded silently. Expected
mapping (BODC series 2009640, AMT24 AFC): **`P700A90Z` = *Synechococcus*,
`P701A90Z` = *Prochlorococcus***, `PYEUA00A` = picoeukaryotes (< 3 µm). Note the
non-obvious order: 700 is Syn, 701 is Pro.

File facts the AFC reader must respect (checked on all four cruises):

- Column order differs between cruises (AMT23 has no `Bot_Ref`; AMT28 adds
  bacterial parameters). Map columns by header code, never by position. Every
  parameter column is followed by its own `QV:BODC` flag column, so the flag
  columns share one name and are paired by position relative to their
  parameter.
- Time is `dd/mm/yyyy hh:mm` GMT and is the **cast start**, not the bottle time
  (see "Sample time" below).
- AMT23 abundances are fractional; AMT24/25/28 are rounded to whole cells mL⁻¹.
- No `<` (below-detection) qualifier appears in any file. Below-detection is
  represented by **zero**; `N` marks nulls.
- Bottle flags other than 0 occur (codes 3, 4, 5, 6, 9: closure failure, firing
  order, leak, mis-seated seals, open vent). Samples with a non-zero bottle flag
  are excluded.

- **Depth cut (Q&A #33, confirmed):** `depth ≤ 10 m`, shallowest bottle per
  station; both window and rule remain parameters, but these are the defaults
  and the values used for every reported result.
- **Daylight:** pre-dawn casts have no radiometry within the window and drop out
  at matchup (expected ≈ half the stations, §1.1). The attrition table reports
  them as a separate line ("no daylight radiometry"), not as a QC failure.
- **Sample time.** The AFC time is the cast start; the surface bottle closes at
  the end of the upcast. Against the cruise report's time-on-deck, the offset is
  0.72 h median (0.67–0.78 h, 16–84 %; two outlier casts). The FCM surface
  sample time is therefore **AFC time + 45 min**, and the binning window below
  is centred there (a parameter, recorded in provenance).
- **Binning (Q&A #32):** for each FCM sample, take the **median** of screened
  Rrs spectra within **±15 min**, per band; keep the 16–84 % spread as a
  per-band uncertainty and the n-in-bin as a weight/QC field. Window and
  estimator are parameters. *Claude Fun additions:* σ_bin(λ) = (p84 − p16)/2
  per band, combined with the floor in §5.2; the indices of the member 1-min
  spectra are stored so training can sample individual spectra (§7.2).
- **Lange-strict mode:** Lange matched by *exact date-hour-minute* to the
  1-min-selected spectrum. For the Tables 1–2 reproduction we support
  `window=0` (nearest 1-min spectrum) alongside our robust default, so the
  reproduction is apples-to-apples while our own training uses the sturdier
  estimator. With `window=0` there is no bin spread, so σ_Rrs is the §5.2 floor
  alone.
- **SST (Q&A #34, confirmed):** primary = `uway_sst` from the Jordan netCDF,
  linearly interpolated to the bin's mid-time; the HSAS ancillary SST column is
  a cross-check only (it is corrupt early in the cruise). The pipeline computes
  and stores both, plus their difference as a QC diagnostic.
- **Qualifier flags and detection limits (settled from the AFC files).** No
  numeric DL is published, but the files encode it. Abundances are event counts
  divided by the analysed volume V, so one event is the smallest non-zero
  abundance and sets the DL. AMT23's fractional values show the step directly:
  single-event values cluster at 1.31–1.40 cells mL⁻¹ (median 1.37). For the
  rounded files, a likelihood fit of `round(k · step)` to all values below
  20 cells mL⁻¹ recovers the step; the signature is that value 2 never occurs on
  AMT24 or AMT25 while 1 and 3 are common, exactly as rounding multiples of
  ≈ 1.37 predicts.

  | cruise | step = DL (cells mL⁻¹) | V (mL) | evidence |
  |---|---|---|---|
  | AMT23 | 1.37 (per-sample 1.31–1.40) | 0.73 | fractional values, direct |
  | AMT24 | 1.38 | 0.72 | fit; Δlog L ≈ 11 vs step = 1 |
  | AMT25 | 1.29–1.40 (two near-equal fits) | 0.72–0.78 | fit; Δlog L ≈ 6 vs step = 1 |
  | AMT28 | 0.67 | ≈ 1.5 | fit; Δlog L ≈ 6 vs step = 1; lower confidence |

  V ≈ 0.73 mL corresponds to about 12 min at the FACSort HI flow rate
  (60 µL min⁻¹). The reader stores these values in `DEFAULT_PIPELINE` and
  re-runs the fit as a regression test. Rules:
  1. Zero → left-censored at DL (§5.4). `N` → missing. Non-zero bottle flag →
     excluded.
  2. Counts `n = round(abundance / step)` feed the counting-noise term (§5.4).
  3. Censoring is active: among AMT24 daylight surface samples, 5 have Pro = 0
     (the high-latitude ends of the transect); among the matched transfer
     samples, 1 on AMT23 and 1 on AMT28. Syn and Peuk have no zeros in either
     set.
  4. AMT28's V is the least certain and affects one matched sample; it is a
     low-priority item on the PML follow-up list
     (`requests/PML_follow_up.md`).

**Known coverage gap:** the BODC deposit is CTD-bottle-only, while Lange's
training set (n = 73–78) also used the ~30-minute **underway** FCM samples —
the very samples the paper shows are what make *Synechococcus* retrievable
(CTD-only retraining degraded Syn MAE 1.27→1.37). Per Q&A #35 this is item 1
on the PML follow-up list (`requests/PML_follow_up.md`); until it is resolved,
our comparison has a known, quantified handicap — and Lange's own CTD-only rows
(their Table 2 sensitivity) are the fair comparison line.

## 4. Pipeline products

The pipeline emits one netCDF/parquet per run under
`<AMT_ROOT>/AMT24/derived/claude_fun/` (the repo stays data-free): the screened
1-min Rrs stream (with per-stage QC flags, ρ_sky, L_NIR), and the matched
training matrix (spectra × 141 native bands with σ_Rrs, member-spectrum
indices, FCM counts with qualifier flags, censoring flags and DLs, SST,
station metadata, provenance). Every stage records counts in/out, so the
attrition table (raw → tilt → geometry → glint-selected → QC → matched, with the
daylight drop-out on its own line) is a standard output — that table is the
first thing to compare against Lange's n. The run also records the NIR-residual σ_abs
diagnostic (§5.2) and the inferred V and DL per taxon (Stage 5).

## 5. Model input

### 5.1 Wavelength range

Default 400–700 nm: the native bands 400.5–698.0 nm (86 bands). Reasons: Stage 3
guarantees non-negative Rrs only there; 750–800 nm is consumed by the glint fit
(Lw ≈ 0 there by construction, so it carries glint-correction residuals rather
than water signal); the UV end is the least reliable part of above-water
radiometry. The range is a config parameter; 414–660 nm is run as an ablation
for MOANA parity.

### 5.2 Tokens

Each band becomes one token. σ_Rrs is the Rrs uncertainty; no bandwidth
feature (decision 2026-10-30). Per-token features:

- `r = Rrs / s_R`, where `s_R` is a single global scalar (training-fold median
  Rrs over all bands and samples). **Do not normalise per spectrum**: absolute
  magnitude carries biomass information.
- `u = σ_Rrs / s_R`, with

  `σ_Rrs(λ) = max( σ_bin(λ), σ_floor(λ) )`,
  `σ_floor(λ) = sqrt( (f_rel·Rrs(λ))² + σ_abs² )`, with
  **f_rel = 0.027** and **σ_abs = 8.2×10⁻⁵ sr⁻¹**.

  - Source: a robust fit of this two-term model to the per-band Rrs
    uncertainties in the Brewin et al. 2023 compilation (PML HyperSAS,
    AMT23/25/26/28; 101 stations × 151 bands, 400–700 nm @ 2 nm, after
    converting the file's signed percentages to absolute values). The two terms
    describe the data across the spectrum: median absolute uncertainty is
    ≈ 2.9×10⁻⁴ sr⁻¹ at 400–420 nm (relative term dominant, ≈ 3 %) and levels off
    at ≈ 7.6–7.9×10⁻⁵ sr⁻¹ from 600 nm onward (absolute term dominant, where
    Rrs ≈ 1–3×10⁻⁴ sr⁻¹). Per-cruise red-band medians range 6.2–9.0×10⁻⁵ sr⁻¹.
    The fit's log-residual scatter is about ±0.65 (16–84 %), so individual
    stations vary by roughly a factor of two around the floor.
  - Why these values: σ is a model input, so the training spectra (AMT24) and
    the transfer-test spectra (AMT23/25/28, which carry these Brewin
    uncertainties, §8.5) must have σ on the same scale. A floor fitted to the
    transfer data guarantees that. Lin et al. 2022 report larger clear-sky
    averages for continuous HyperSAS measurements on AMT28 (≈ 6 % below 490 nm,
    ≈ 12 % above 550 nm). Above 550 nm most of that is the absolute term divided
    by small Rrs, which the two-term model represents explicitly; the blue-band
    difference (≈ 3 % here vs ≈ 6 %) is not resolved by the files at hand, and
    σ_bin exceeds the floor whenever conditions are worse.
  - The max (not a quadrature sum) with σ_bin avoids double counting, because
    the bin spread already samples environmental variability.
  - Diagnostic only: the pipeline also computes the robust spread
    (1.4826 × MAD) of glint-corrected AMT24 Rrs over 750–796 nm, where
    Lw ≈ 0 in open-ocean water. If it differs from 8.2×10⁻⁵ sr⁻¹ by more than a
    factor of two, stop and report at checkpoint 1.
  - Sanity check: the Stage 2 cross-check channel (provider fixed-ρ Rrs minus
    our fitted Rrs) gives an independent measure of how much of the uncertainty
    comes from the glint correction; report its per-band spread next to
    σ_floor.
- λ encoding: Fourier features `sin(2π(λ−300 nm)/P_k)`, `cos(…)` for 16
  log-spaced periods `P_k` from 5 nm to 500 nm (32 features). Continuous in λ,
  so any grid is valid input.

A spectrum is a variable-length set of tokens; batches are padded with an
attention mask. Nothing in the model depends on the number of bands or their
order.

### 5.3 SST

SST is an input (decision 2026-10-30). It enters as a single global feature:
`uway_sst` in °C, z-scored with training-fold statistics, concatenated to the
pooled encoder output `h` before the Bayesian head (§6.2). Linear °C rather than
MOANA's log10(SST) keeps the input defined at SST ≤ 0 °C in later global data;
the head can learn any needed nonlinearity. An SST-off run is kept as an
ablation (§8.3), and §8.2 measures how much the model relies on SST.

For AMT23/25/28 (§8.5), SST is the Brewin column `Temperature above Secchi
depth [degC]` (valid for all 101 usable stations; CTD profiler for 97,
underway for 4). This is a near-surface layer mean rather than underway SST,
a known source difference that is recorded per sample. If the Jordan et al.
netCDFs for those cruises are added later (the Jordan compilation spans
AMT19 and AMT22–29), switch to their `uway_sst` for consistency with training.

### 5.4 Labels

`y_k = log10(cells mL⁻¹)` for k ∈ {Pro, Syn, Peuk}, z-scored per taxon using
training-fold statistics. How censored and missing taxa enter the joint likelihood is specified in §6.3.

Counting noise: V is recovered for all four cruises (Stage 5), so each label has a
Poisson counting uncertainty `σ_count = 0.434/√n` in log10, with
`n = abundance × V` events. It is added to the diagonal of the predicted covariance
inside the NLL (`Σ_eff = Σ + diag(σ_count²)`), so the model's Σ describes
retrieval uncertainty rather than counting noise. The term matters only near the DL: at
100 cells mL⁻¹ (≈ 73 events) it is ≈ 0.05 dex.

## 6. Model architecture

### 6.1 Set encoder (deterministic)

Token MLP (features → 64) → 2 transformer encoder blocks (d = 64, 4 heads,
feed-forward 128, pre-LayerNorm, dropout 0.1, masked attention) → attention
pooling (one learned seed query attending over the tokens) → pooled vector
`h ∈ ℝ⁶⁴`. No positional encoding beyond the λ features, so the encoder is
permutation-invariant by construction. The standardised SST is concatenated to
`h` (§5.3).

### 6.2 Bayesian head

Two variational dense layers (64 units, GELU) and a variational output layer,
each with a mean-field Gaussian posterior over weights and biases, sampled with
Flipout or the local reparameterisation trick. The output layer gives the mean and
covariance of the three taxa:

- `μ_k` and `σ_k = softplus(s_k) + σ_min` (σ_min = 0.01 in standardised units,
  config), k ∈ {Pro, Syn, Peuk};
- a 3×3 correlation matrix R, built from three canonical partial correlations
  `ρ₁₂ = tanh(z₁)`, `ρ₁₃ = tanh(z₂)`, `ρ₂₃|₁ = tanh(z₃)`, which guarantees a valid
  (positive-definite) R for any z;
- `Σ = D R D` with `D = diag(σ_k)`.

`DEFAULT_MODEL["corr_mode"]` sets where z comes from:

- `"shared"` (**default this stage**): z is three variational parameters of
  the output layer with no input dependence. Each sample has its own σ_k, but
  the correlation structure is global. Three extra parameters are what ~33
  training labels can constrain.
- `"input"`: z is predicted per sample (3 more outputs), so correlations can
  vary with water type. Run as an ablation now; the candidate default for the
  next stage, when labels increase (§14).
- `"diag"`: R = I. Ablation, and the independent-taxa reference.

Reported products are the marginals (μ_k, σ_k = √Σ_kk) plus R.

The encoder is deterministic because variational attention is costly and
unstable at this data size; consequence: epistemic uncertainty from encoder
weights is not captured. A cheap option is MC dropout in the encoder at
prediction time (config flag, off by default).

### 6.3 Loss

The training objective is the negative ELBO:

`L = (1/N) Σᵢ NLLᵢ + (β/N) · KL( q(w) ‖ p(w) )`

- `NLLᵢ`: the multivariate Gaussian NLL of the label vector under
  `N(μ, Σ_eff)` in standardised log space, with one weight sample per step.
  - Missing taxa: marginalise by dropping their rows and columns of μ and Σ.
  - Censored taxa (below DL): for observed set O and censored set C,
    `NLLᵢ = −log N(y_O; μ_O, Σ_OO) − log P(y_C < log10 DL_C | y_O)`, where the
    conditional is Gaussian with mean `μ_C + Σ_CO Σ_OO⁻¹ (y_O − μ_O)` and
    covariance `Σ_CC − Σ_CO Σ_OO⁻¹ Σ_OC`. With one censored taxon (the only
    case in the current data, always Pro, Stage 5) the last term is a
    univariate Φ. Samples with two or more censored taxa do not occur in the
    AMT24 or transfer sets; the code excludes them and counts them in the
    attrition table rather than approximating a bivariate CDF.
  - With `corr_mode="diag"` this reduces exactly to the sum of independent
    per-taxon terms.
- The KL is between the variational posterior over **head weights** `q(w)` and
  a **weight prior** `p(w) = N(0, σ_w² I)`, σ_w = 1. There is no prior on the
  cell counts (decision 2026-10-30).
- `N` is the number of **independent labels** in the training fold, not the
  number of augmented spectra. Replicating spectra (§7.2) adds no independent
  information; dividing the KL by the replicated count would make the
  posterior overconfident.
- `β = 1`, with optional linear warm-up over the first 10 % of epochs.

### 6.4 Prediction

T = 200 weight samples → `(μ_t, Σ_t)`. Report, in log10 units: predictive
mean `mean(μ_t)`, aleatoric covariance `mean(Σ_t)`, epistemic covariance
`cov(μ_t)`, and total covariance (their sum). From the total covariance, report
per-taxon σ_k (the MOANA-comparable product) and the correlation matrix.
Output is an `xarray.Dataset`.

### 6.5 Implementation notes

TensorFlow with Keras 3. Custom Keras layers must subclass `keras.layers.Layer`;
this is the one permitted exception to "methods, not classes" (§9). Implement
the variational dense layer directly (~50 lines) rather than depending on
`tfp.layers`, which is built for Keras 2 (`tf_keras`). All architecture
hyperparameters live in `DEFAULT_MODEL` in `algorithm.py`.

## 7. Training

### 7.1 Defaults (`DEFAULT_TRAIN` in `train.py`)

Adam, learning rate 1×10⁻³ with cosine decay; batch 16; at most 3000 epochs;
early stopping on inner-fold NLL with patience 200; fixed seeds recorded in
provenance; optional 5-seed ensemble.

### 7.2 Augmentation (each a config flag; on for training, off for evaluation)

- **Replicate sampling:** each epoch, for each label, draw one member 1-min
  spectrum from its ±15 min bin (with σ_Rrs as in §5.2) instead of always using
  the median.
- **Noise:** `Rrs ← Rrs + ε·σ_Rrs`, ε ~ N(0, 1).
- **Spectral degradation** (required for the "any resolution" goal): random
  band dropout (keep fraction ∈ [0.3, 1]); convolution with a Gaussian spectral
  response of FWHM ∈ [3.5, 15] nm and resampling to a random regular grid with
  spacing ∈ [2, 10] nm; σ propagated as `(Σ_j w_j² σ_j²)^½` (band-independent
  errors, which is optimistic, hence the σ floor). Without this, the model only
  ever sees one 86-band grid; the set encoder makes other grids *acceptable*
  but nothing teaches the model to handle them.

## 8. Validation

### 8.1 Protocol

`moana_design.md` defines two MOANA validation targets that matter here:
(i) retrain on AMT24 and reproduce Lange Table 1 (hyperspectral, in-situ) as
bias/MAE/R², which is an **in-sample** comparison, with the CTD-only rows as
the fair line; and (ii) apply the published model to the held-out cruises
AMT23/25/28. In-sample numbers are meaningful for MOANA's few-parameter linear
regressions but not for a flexible network, so Claude Fun is compared in three
ways, in order of weight:

1. **Held-out cruises (primary).** Claude Fun trained on all matched AMT24
   labels vs published MOANA (target ii) and our AMT24-retrained MOANA, all
   evaluated on AMT23/25/28 (§8.5). Everything is out-of-sample for every
   model.
2. **Cross-validation on AMT24, identical folds for all models.**
   Leave-one-station-out (one label per station, so equal to leave-one-out),
   and latitude-blocked 5-fold (contiguous blocks along the transect), because
   neighbouring stations are spatially autocorrelated and leave-one-out
   overstates skill on unseen water masses. All fitted statistics (`s_R`, label
   and SST standardisation) are computed inside each training fold.
3. **Lange Table 1 (reference only).** Reported next to Claude Fun's
   cross-validated numbers with the explicit note that Lange's are in-sample
   and Claude Fun's are not.

### 8.2 Metrics

Per taxon (from the marginals, so directly comparable with MOANA):


- Point skill: bias, MAE and R² as defined for Lange Table 1, using the same
  implementation as the MOANA reproduction (so the 1.27/1.37 Syn figures are
  comparable), plus RMSE in log10.
- Probabilistic: mean NLL, CRPS, empirical coverage of the 68 % and 95 %
  intervals, PIT histogram, reliability plot; aleatoric/epistemic split.
- SST reliance: permutation importance (shuffle SST across held-out samples,
  report ΔNLL and ΔMAE), alongside the SST-off ablation. Note on reading the
  weights directly: in a variational head, weights on an uninformative input
  revert towards the prior N(0, σ_w²) rather than to zero, which shows up as
  extra epistemic variance. Use these two diagnostics, not weight magnitudes,
  to judge whether SST is used.

Joint (do the correlations help?):

- Joint NLL, and coverage of the 68 % and 95 % joint ellipsoids (squared
  Mahalanobis distance against χ² with 3 degrees of freedom).
- Energy score, the multivariate analogue of CRPS.
- These are compared across `corr_mode` = shared / input / diag (§8.3); per-taxon
  metrics should be nearly unchanged between them, and the joint metrics show
  whether the correlation is real.

### 8.3 Baselines and ablations

1. Climatology null: training-fold mean and standard deviation.
2. `epft_up.moana` retrained on the same folds (CTD-only configuration).
3. Published MOANA (vendored LUTs) on AMT23/25/28.
4. Lange Table 1 and Table 2 CTD-only rows (reference, §8.1).
5. Ablations: `corr_mode` = diag and input (vs shared); deterministic head (no KL); no spectral-degradation
   augmentation; SST off; 414–660 nm range; counter-binning input (spectra
   interpolated to a fine fixed grid with σ inflated by the counter-binning
   factor, into an MLP head), which tests whether the attention encoder earns
   its complexity.

### 8.4 Resolution robustness

Re-run the held-out predictions with inputs degraded to (a) MOANA 414–660 nm @
2 nm, (b) an OCI-like grid (Gaussian 5 nm FWHM, 2.5 nm sampling; replace with
official OCI spectral response functions when wired in), (c) a multispectral
subset (412, 443, 490, 510, 555, 670 nm). Report `|Δμ|/σ_total` relative to the
native-grid prediction and the change in §8.2 skill. Proposed pass criterion:
median `|Δμ|/σ_total < 0.25` for (a) and (b).

### 8.5 Transfer to AMT23/25/28

Labels: the three BODC AFC files (§2), through the same Stage 5 rules (≤10 m,
shallowest bottle, qualifier flags, DL inference per cruise).

Radiometry: the Brewin et al. 2023 compilation. Its contents, read from the file:

- 127 rows, one optics station per day near local solar noon (solar time
  8.4–12.7 h, median 11.2 h), from AMT23 (25 rows), AMT25 (35), AMT26 (37) and
  AMT28 (30). **There are no AMT24 rows**, so no same-cruise cross-check against
  our pipeline is possible. AMT26 has no AFC file and is unused unless one is
  added.
- Rrs on a 400–700 nm, 2 nm grid (151 bands), in two versions: BRDF-corrected
  (Lee et al. 2011) and uncorrected. **Use the uncorrected columns**
  (`Rrs(λ) [1/sr]`): the AMT24 pipeline applies no BRDF correction, and mixing
  the two would add a systematic shift. The 2 nm grid differs from AMT24's
  3.5 nm grid, so the transfer test is also a real-data grid change.
- Per-band uncertainty `Rrs(λ) Uncertainty [%]`, relative and signed (negative
  where Rrs < 0, and very large where Rrs ≈ 0). Convert to absolute,
  `σ_file(λ) = |U(λ)/100 · Rrs(λ)|`, then `σ_Rrs = max(σ_file, σ_floor)`,
  mirroring training. Where U is missing, σ_floor alone.
- Missing-value sentinel −999 for both Rrs and uncertainty. Rows are
  all-or-nothing: 101 rows have all 151 bands, 26 have none. Usable rows on the
  transfer cruises: AMT23 21, AMT25 26, AMT28 25, so **at most 72 transfer
  samples** before the FCM matchup.
- Six usable transfer rows (AMT25 1, AMT28 5) have small negative red Rrs
  (down to about −2.7×10⁻⁴ sr⁻¹). Stage 3 QC (no negative Rrs in 400–700 nm;
  second-derivative filter, native to this 2 nm grid) is applied exactly as for
  AMT24, with its attrition reported.
- Ancillary QC columns (visual cloud index, visual glint index, pitch and roll
  standard deviation, Li(750)/Es(750)) are carried through for stratified
  reporting, not used as filters.
- SST: see §5.3.

Matchup: each Brewin row is paired with the shallowest ≤10 m AFC sample
(sample time = AFC time + 45 min, Stage 5) from a CTD cast within
|Δt| ≤ 1.5 h and ≤ 10 km (the ship is on station for the noon optics and CTD
work); if several casts qualify, take the nearest in time. Δt and distance are
stored per pair, and the attrition (Brewin rows → QC → matched) is reported per
cruise.

Result (settled from the files): every usable Brewin row but one per cruise
matches a distinct cast, with |Δt| ≤ 1.4 h and distance ≤ 5.6 km, and the counts
are identical for windows from ±1.5 h to ±6 h:

| cruise | usable Brewin rows | matched | negative red Rrs (fail Stage 3) | after negative-Rrs QC |
|---|---|---|---|---|
| AMT23 | 21 | 20 | 0 | 20 |
| AMT25 | 26 | 25 | 1 | 24 |
| AMT28 | 25 | 24 | 5 | 19 |
| **total** | 72 | **69** | 6 | **63** |

The second-derivative filter may remove a few more; its attrition is reported.

### 8.6 SeaFlow + PACE test (minimal spec; full processor in next document)

SeaFlow v2 Pro/Syn/picoeukaryote abundances from cruises after PACE launch
(2024), matched to PACE OCI L2 Rrs: ±3 h, 5×5 pixel box centred on the ship
position, at least half the box valid after standard L2 flags, box median as
the spectrum and box spread (combined with the L2 Rrs uncertainty product if
available) as σ_Rrs; SST from the SeaFlow underway record. Expected domain
shifts: satellite atmospheric correction vs above-water radiometry, and surface
underway sampling vs ≤10 m CTD bottles. Test only; no training on these pairs in
this stage.

### 8.7 Spectral attribution (activation maps)

Goal: show which wavelengths drive each output (μ_k and σ_k, per taxon), to
tell broad-shape or water-mass proxies apart from localised features such as
pigment bands. All maps are computed on held-out samples only (CV hold-outs,
transfer sets) and averaged over the T weight samples,
so each map carries an uncertainty band.

Methods:

1. **Band-window occlusion (primary).** Delete all tokens inside a sliding
   window (10 nm wide, 5 nm step, across the model's range) and record Δμ_k,
   Δσ_k and, where labels exist, ΔNLL. Because the encoder is a set model trained
   with band dropout (§7.2), deleting tokens is an in-distribution perturbation
   with no arbitrary fill value. Maps from different grids (AMT24 3.5 nm, Brewin
   2 nm, PACE 1.25/2.5 nm) land on the same window axis.
2. **Integrated gradients (secondary).** Per-token attribution of μ_k with
   respect to `r` (§5.2), from a baseline equal to the training-fold mean
   spectrum on the same λ grid, 64 integration steps. Check completeness: the
   attributions sum to μ(x) − μ(baseline) within 1 %. Divide by the local band
   spacing so maps from different grids are per-nm densities.
3. **Attention weights are not attribution.** The pooling weights are reported
   as a diagnostic only; a large attention weight is not evidence of influence.

References:

- **MOANA's exact sensitivity**: the Jacobian ∂ŷ_k/∂Rrs(λ) of published and
  retrained MOANA, computed analytically through standardise → project →
  evaluate (`moana_design.md` §3). MOANA standardises each spectrum (subtracts
  its mean, divides by its standard deviation), so it is invariant both to
  scaling the spectrum and to adding a constant to it. Its Jacobian therefore
  satisfies Σ_λ J(λ)·Rrs(λ) = 0 and Σ_λ J(λ) = 0 — two built-in correctness
  checks.
- **Pigment features as a guide, not a criterion**: chlorophyll-a absorption
  near 440 and 675 nm; phycoerythrin absorption in *Synechococcus* around
  495–565 nm.

Sanity checks:

- Model randomisation: maps from a model with re-initialised head weights must
  lose their structure; if they do not, the method is reflecting the input, not
  the model.
- Stability: report the spread of maps across CV folds and seeds. Features that
  do not survive across folds are not interpreted.
- Known answer: run both methods on retrained MOANA itself, treated as a
  black-box model, before applying them to Claude Fun. Integrated gradients
  must recover MOANA's analytic Jacobian (correlation > 0.9 over 414–660 nm).
  Occlusion measures a finite deletion (MOANA interpolates across the gap)
  rather than a derivative, so it must instead place its largest effects in
  the same windows as the Jacobian.

Interpretation limits: a map shows what the model uses, not why. A wavelength
can matter because it tracks a water mass that co-varies with abundance. With at
most 33 AMT24 labels, maps will be noisy; evidence that the network learned
cell-specific spectral features belongs to the next stage. Maps are stratified
by latitude band where n allows.

Outputs: per-taxon maps (mean and fold spread) for μ_k and σ_k, with MOANA's
Jacobian overlaid, as an `xarray.Dataset` plus report figures.

## 9. Code layout and conventions

```
epft_up/claude_fun/
    cf_design.md     # this document
    __init__.py      # public surface: run_cf(), train_cf(), predict_cf(), validate_cf()
    io.py            # .sav reader, FCM/BODC reader (flags, DL inference), Jordan-SST reader, Brewin reader
    pipeline.py      # stages 1–5: screening, glint, QC, QC-resampling, matchup; DEFAULT_PIPELINE
    algorithm.py     # tokenisation, standardisation, encoder, Bayesian head; DEFAULT_MODEL
    train.py         # ELBO training, augmentation, CV folds; DEFAULT_TRAIN
    validation.py    # metrics, baselines, robustness and transfer tests
    attribution.py   # occlusion, integrated gradients, MOANA Jacobian (§8.7)
    reference/lange2020_tables.csv   # published numbers
epft_up/data/moana/  # vendored LUTs (already in place, sha256-pinned; baseline only)
epft_up/tests/test_cf.py
```

Per the project conventions and Q&A decisions: **methods, not classes** — the
natural state objects are plain dicts/`xarray.Dataset`s passed explicitly
(exception: custom Keras layers, §6.5); imports at top; every method's
docstring lists inputs/outputs; inline comments explain the *why* (especially
every Lange threshold, each of which cites the paper section). All thresholds
and paths live in module-level config dicts (`DEFAULT_PIPELINE`, including
`AMT_ROOT`; `DEFAULT_MODEL`; `DEFAULT_TRAIN`), so a run's provenance is a
serialisable dict (all three, plus git hash, input file hashes and seeds),
matching the IOPtics provenance habit. Plotting (matplotlib, `Agg`) goes in the
report scripts, not the package.

## 10. Tests (`epft_up/tests/test_cf.py`)

- Pipeline: `.sav` column layout and wavelength axis; BODC code → taxon mapping
  equals the expected table in Stage 5 (`P700A90Z` → Syn, `P701A90Z` → Pro,
  `PYEUA00A` → Peuk); qualifier flags parsed into censored/missing; V recovered
  from synthetic quantised counts, and the fallback triggered on unquantised
  ones; tilt formula; Δϕ convention check; L1 glint fit recovers known ρ_sky,
  L_NIR on synthetic spectra and respects bounds; second-derivative filter on a
  synthetic noisy spectrum; median/16–84 % binning; σ_floor formula and the
  max with σ_bin; σ_abs diagnostic on synthetic NIR residuals; attrition counts
  non-increasing stage to stage, with the daylight line present.
- AFC reader: header-based mapping on all four cruises (including the paired
  `QV:BODC` columns); non-zero bottle flags excluded; `N` → missing; zero →
  censored; on the real files (`$OS_COLOR`-gated), the step fit reproduces the
  Stage 5 table (1.37 / 1.38 / 1.29–1.40 / 0.67) and the AMT24 daylight count
  (35 daylight, 33 on radiometry days).
- Brewin reader: −999 → missing for Rrs and uncertainty; signed % uncertainty
  converted to absolute; the uncorrected (not BRDF-corrected) Rrs columns
  selected; on the real file (`$OS_COLOR`-gated), 101 full-spectrum rows with
  21/26/29/25 for AMT23/25/26/28; the σ_floor fit reproduces f_rel ≈ 0.027 and
  σ_abs ≈ 8.2×10⁻⁵ sr⁻¹.
- Model: output invariant to token permutation; accepts variable band counts in
  one padded batch; `σ_k > 0`; R positive-definite for random z (including
  extreme values); the marginal of the joint Gaussian equals the univariate
  per-taxon distribution; with `corr_mode="diag"` the joint NLL equals the sum
  of univariate NLLs; the censored conditional term matches numerical
  integration; missing-taxon marginalisation; counting noise added to the
  diagonal only when enabled; KL scaled by
  independent-label N; predictive covariance decomposition sums correctly; loss
  decreases on a tiny synthetic dataset; fixed seed reproduces outputs.
- Attribution: occlusion with an empty window changes nothing; integrated
  gradients satisfy completeness within 1 %; the analytic MOANA Jacobian matches
  finite differences and satisfies Σ_λ J(λ)·Rrs(λ) = 0 and Σ_λ J(λ) = 0; on a toy model with a
  known single-band dependence, both methods put their mass on that band.

## 11. Dependencies

`numpy`, `scipy` (`io.readsav`, `optimize`), `pandas`, `pyarrow` (parquet),
`h5py`, `xarray`, `netCDF4` (xarray backend), `tensorflow` (≥ 2.16, ships
Keras 3), `earthaccess` (PACE, §8.6 only), `matplotlib` (report scripts only),
`pytest`.

## 12. Deliverables and implementation order

Implement and stop for review at each checkpoint:

1. Stages 0–5 with tests → **checkpoint:** attrition table vs Lange's n (with
   the daylight drop-out; at most 33), the matched AMT24 matrix, the DL fit
   reproducing the Stage 5 table, and the NIR-residual σ_abs diagnostic against 8.2×10⁻⁵ sr⁻¹.
2. Model and training with tests → **checkpoint:** leave-one-station-out
   results for Claude Fun and baselines 1–2 on one config.
3. Full AMT24 validation (§8.1–8.4, §8.7) → **checkpoint:** validation report
   (metrics table, calibration plots, SST reliance, robustness table, the
   attribution known-answer check on MOANA, and attribution maps with MOANA's
   Jacobian overlaid).
4. Transfer tests (§8.5–8.6), with attribution maps on the transfer sets.

Do not change any threshold or default silently; anything not settled here
goes back to Billy as a question.

## 13. Open items

None blocking. Settled 2026-10-30: data root (`$OS_COLOR/AMT`); Brewin
compilation contents and transfer matchup (§8.5); σ floor values (§5.2);
detection limits and analysed volumes (Stage 5); AMT24 matched-n upper bound
(§1.1). On the PML follow-up list, non-blocking: AMT24 underway FCM samples
(Q&A #35) and confirmation of AMT28's analysed volume.

## 14. Next stage (next design document)

- Train on AMT and SeaFlow+PACE together, with a data processor that ingests
  the different input forms (in-situ hyperspectral, satellite L2, multispectral)
  into the common (λ, Rrs, σ_Rrs) + SST representation.
- **Input-dependent correlations.** Decide whether `corr_mode="input"`
  becomes the default, based on this stage's joint metrics and the larger
  label set.
- Revisit model size once the number of independent labels is known.
- Add the AMT24 underway FCM samples if PML supplies them (Q&A #35).
