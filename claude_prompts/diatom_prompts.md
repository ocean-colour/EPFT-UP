# Diatom carbon biomass

## Goals

Explore the new diatom carbon biomass algorithm provided by Alison Chase.

## Context

The GitHub repo is `https://github.com/alisonchase/diatom_carbon_biomass`.
I have put a local copy on this machine at `Oceanography/python/diatom-carbon-biomass-model`

## Prompts

### Explore

1. Explore the code in the GitHub repo.  Then generate a report of the code.
Also, if you have any questions, ask them in the Q&A section.  Use Fable if you can.  Log your work.

## Reports

### Report 1 — Code review of `diatom-carbon-biomass-model` (2026-09-15)

#### 1. What this repository is

Five Python modules (3,287 lines total) implementing an empirical diatom carbon
biomass algorithm: a **Random Forest regressor that maps four HPLC pigment
concentrations to diatom carbon**, with an ONNX export path and a companion
uncertainty/domain bundle. Supporting the paper Chase et al. (preprint).

Actual remote: `https://github.com/alisonpchase/diatom-carbon-biomass-model`
(not `alisonchase/diatom_carbon_biomass` as recorded in the Context section
above). Two commits, both 2026-09-04, author "Ali Chase". Licence **CC0 1.0**
(public domain dedication) — notably permissive, and different from the BSD-3
used in this project's sibling repos.

| File | Lines | Role |
| --- | --- | --- |
| `model_config.py` | 114 | Single source of truth: features, ratios, log list, target, RF hyperparameters, plot colours, cruise names |
| `prepare_features.py` | 164 | Ratio derivation, log10 transform, validity masking. Shared by every other script |
| `train_evaluate_model.py` | 1,766 | Hold-out + cross-validation evaluation, metrics, 10 diagnostic figures |
| `train_model_final_export.py` | 447 | Fits the production model on all data; exports ONNX + joblib bundle + JSON metadata |
| `run_model_inference.py` | 555 | CSV-in/CSV-out prediction with per-tree uncertainty and domain flags |

**Not in the repository:** `trained_models/` and `training_data/` are both
git-ignored, and `*.onnx`, `*.joblib`, `*.csv` are routed to Git LFS. So the
README's "Contents" list advertises a `trained_models/` directory that is not
distributed, and there is no training data, no exported model, and no test
suite. Nothing here can be executed end-to-end as cloned.

#### 2. The algorithm

**Inputs** — four pigment concentrations (mg m⁻³): `chla`, `chlb`, `chlc`
(c1+c2), `ppc` (photoprotective carotenoids). Three ratios are derived
(`chlb_a`, `chlc_a`, `ppc_a`, each pigment over chl-a), then **all six model
features are log10-transformed**. Note `ppc` itself is not a feature — only its
ratio. `t_avg` and `s_avg` (temperature, salinity) are present but commented
out in `model_config.py:14-15`, annotated with importances of 0.2264 and 0.0930.

**Target** — `diatCarb`, log10-transformed, fit under squared-error loss.
Back-transformed with a plain `10**`.

**Model** — `RandomForestRegressor`, hyperparameters hardcoded in
`model_config.py:50-61`: `n_estimators=1500, max_depth=15, max_leaf_nodes=500,
min_samples_split=10, min_samples_leaf=5, max_features='sqrt', max_samples=0.8,
random_state=42`. The inline comments say "optimized based on hyperparameter
testing", but **no search code exists anywhere in the repository** — no
GridSearch, no Optuna, no param grid. The provenance of these values is
undocumented.

**Training data** — in-situ cruise matchups; the cruise map in
`model_config.py:99-115` names 16 cruises: NAAMES 01-04, EXPORTS 01-02,
PEACETIME, Tara Mission Microbiomes 1-6, SO-PACE KM2418/KM2419/TGT444. The
README's example filename (`ifcb_enviro_09Mar26.csv`) indicates the target is
IFCB-derived diatom carbon.

**There is no satellite code in this repository.** No xarray, netCDF, h5py, no
Rrs or wavelengths, no L2 flags, no geolocation, no chunking. Both the training
and inference paths are tabular CSV. Applying this to PACE would require a
pigment-retrieval step that lives elsewhere, and whose error is not represented
in any uncertainty this code emits.

#### 3. Evaluation design

**Hold-out split** (`train_evaluate_model.py:207-265`) is *cruise-blocked and
chl-a-stratified*: each cruise is assigned to one of five chl-a strata
(`0.001-0.03, 0.03-0.17, 0.17-1, 1-5, 5-30` mg m⁻³) by its **median** chl-a,
then within each stratum the cruise list is split 75/25 with seed 42. Splitting
on cruises rather than samples is the right call and avoids station-level
leakage.

**Cross-validation** (`:1318-1501`) is 5-fold `GroupKFold` grouped by cruise, run
on the training portion only. Per fold it refits and computes permutation
importance on the held-out fold. It is **not nested** — but since no
hyperparameter search happens in-script, there is nothing to nest.

**Metrics** (`:298-391`) are computed mostly in native units: MAE, MedAE, RMSE,
R², Spearman ρ, MAPE, median relative error, MBE, median bias, MPE, plus `r2_log`
in log space. A `learning_curve` diagnostic and a comparison against the Chase
et al. 2022 Eq. 6 chl-a power law (`1.5·chla^1.9`, `:1635`) round it out.

#### 4. Uncertainty and domain flagging

Three different things are called "uncertainty" across the repo, and **they are
not the same estimator**:

1. **Evaluation** (`:418-419`) uses `forestci.random_forest_error` — the
   Wager–Hastie–Efron infinitesimal jackknife. This is the *sampling variance of
   the ensemble mean*, i.e. how well-determined the RF's prediction is, not how
   far a new observation will fall from it.
2. **Inference** (`run_model_inference.py:395-400`) uses the **standard
   deviation and 5/50/95th percentiles across the 1500 trees**. This is ensemble
   disagreement — a different quantity again.
3. Both are converted to native units by a first-order delta method,
   `σ_y ≈ ln(10)·y·σ_log`, which is symmetric about a prediction whose true error
   distribution is strongly right-skewed. The quantile columns are the sounder
   product.

Neither is a calibrated predictive interval, and **no coverage check exists**
anywhere — nothing verifies what fraction of held-out observations actually fall
within ±1σ. For a project whose premise is rigorous uncertainty, this is the
central gap.

**Domain flagging** (a genuinely good feature) is computed at export time: the
per-feature training min/max, plus a Mahalanobis mean, inverse covariance, and a
threshold at the 95th percentile of training distances. Inference flags samples
outside the training range or beyond the Mahalanobis threshold. Two caveats: the
95th-percentile threshold flags 5% of the *training* data by construction, and
see defect D5 below.

#### 5. Defects found

Ranked by consequence. All verified directly against the source; the mechanical
ones are re-checked by `reports/scripts/diatom_code_checks.py`.

| # | Location | Defect |
| --- | --- | --- |
| **D1** | `run_model_inference.py:125-148, 304` | **The exported JSON metadata is never read.** The exporter writes `<stem>_metadata.json`; `load_metadata()` is always passed a *directory*, and the stem-specific filename is only added to the search list when the argument `is_file()` — never true. So only a literal `metadata.json` is sought, which nothing writes. It silently falls back to the ONNX `metadata_props`. |
| **D2** | `run_model_inference.py:313-314` | **Consequence of D1: domain flags never appear.** `has_domain_flags` is gated on `domain_flag_methods`, which exists only in the JSON, not in the ONNX props. Every domain column — `mahal_dist`, `domain_flag_range`, `domain_flag_mahal`, `domain_flag_any` — is silently dropped from the output, and the summary block never prints. The entire domain-flagging feature is dead in practice. |
| **D3** | `run_model_inference.py:310, 526` | **Consequence of D1: the log-target flag is parsed unsafely.** The fallback reads the ONNX string `target_log10`, and `bool("false")` is `True` in Python. It works today only because the flag happens to be true; a linear-target model would have `10**` applied to linear predictions. The dead `load_model_and_config` at `:239` does it correctly with `.lower() == "true"`. |
| **D4** | `run_model_inference.py:306-317` | **Inference ignores the model's embedded feature definitions.** `FEATURE_NAMES`, `VAR_TO_LOG`, `RATIO_DEFS` are taken from the *local* `model_config.py`, never from the ONNX metadata the exporter deliberately wrote. Correctness depends on the config file on the inference machine matching the one used at training. A width-changing edit fails loudly (ONNX shape mismatch); a same-width edit — reordering features, dropping one from `VAR_TO_LOG` — produces **silently wrong predictions**. No assertion guards this. |
| **D5** | `train_model_final_export.py:139-185` | **Mahalanobis distance on a rank-deficient space.** Because ratios are logged after division, `log(chlb_a) ≡ log(chlb) − log(chla)` exactly: the 6-column feature matrix has rank 4 (verified numerically — singular values drop to ~1e-13). The covariance is singular, and the `+1e-6·I` ridge is the only reason it inverts, giving the inverse eigenvalues of 1e6 along two null directions. Training distances are unaffected (their null component is exactly zero), so it is numerically benign today, but the metric is scale-arbitrary in two of six directions and the `pinv` fallback at `:169-171` is unreachable. Compute it on the four independent variables instead. |
| **D6** | `train_evaluate_model.py:1351, 1375` vs `:1739` | **Return-arity mismatch.** `cross_validation_analysis` returns a bare array on two guard paths but a 2-tuple on success; the caller always unpacks two values. Either guard firing raises instead of degrading gracefully. |
| **D7** | `train_evaluate_model.py:1481, 1494` | **Hardcoded `../figures/` paths** for the two cross-validation CSVs, ignoring `--figures-dir` and writing relative to the parent of the working directory. Fails outright if that directory does not exist. |
| **D8** | `train_model_final_export.py:73-87` | **Path resolution makes the README's own example wrong.** A non-existent directory path without a trailing slash is treated as a file, so `--output-model ./trained_models --model-name diatom_model.onnx` writes `./trained_models.onnx` and silently ignores `--model-name` unless `./trained_models/` already exists. |
| **D9** | `run_model_inference.py:481-501` | **Ambiguous positional CLI.** The usage string advertises an optional *middle* argument, which positionals cannot express. With `$DIATOM_INFERENCE_MODEL` set, `run_model_inference.py in.csv model.onnx` treats `model.onnx` as the output path and **overwrites the model with a predictions CSV**. |
| **D10** | `train_model_final_export.py:373-387` | **The ONNX export is undermined by its own companion.** The bundle pickles all 1500 fitted `DecisionTreeRegressor` objects, and inference loads it unconditionally. So portability and pickle-safety — the reasons to ship ONNX — are lost anyway: you still need a version-matched scikit-learn, and loading executes arbitrary code. |

Lesser issues: dead code in `run_model_inference.py` (`onnx_predict`,
`inverse_target_log10`, `load_model_and_config`, and `tree_mean`, which should
really be compared against the ONNX prediction as a runtime parity check);
`run_onnx_predict:153` lacks the `HAS_ORT` guard its dead twin has; the ratio
`if out_col in df.columns: pass` block in `prepare_features.py:60-62` is a no-op;
strata bins, ratio definitions and the 27% in-situ error constant are each
duplicated in two places; `"\\n"` appears in ~11 f-strings where `"\n"` was meant,
printing a literal backslash-n; MBE is `pred − obs` while MPE is `obs − pred`, so
the two bias metrics carry opposite signs; a NaN Spearman ρ is silently replaced
by 0.0 (`:341-343`); and `matplotlib`'s `boxplot(labels=...)` at `:772` is
deprecated in ≥3.9.

#### 6. Scientific concerns

**Retransformation bias.** The model minimises squared error on log10(diatCarb)
and inverts with `10**`, which returns a geometric-mean-like quantity — the
conditional *median*, not the mean. Yet the headline metrics include native-scale
MBE, MPE and MAPE, which compare against the arithmetic mean. The resulting low
bias is not small: at a residual spread of σ_log10 = 0.25 it is 15%, at 0.35 it
is 28%, at 0.50 it is 49% (computed in the verification script). Either apply a
Duan smearing correction, or state explicitly that the product is a median and
report `tree_q50` as the central estimate.

**Single split, 16 cruises.** Test metrics rest on whichever 3-5 cruises seed 42
happened to select. With cruise-blocked data there is no cheap way to make that
robust except repetition: repeated splits, or leave-one-cruise-out, would give a
distribution rather than a point estimate. Related — cruises are stratified by
their *median* chl-a, so a cruise spanning oligotrophic to bloom conditions
(NAAMES) gets one label, while the sample-level `stratum` column computed at
`:198-200` is used only for a bar chart. The plot and the split mean different
things by "stratum".

**Native-scale R² over five decades of biomass** is dominated by a handful of
high-biomass points. The final printout (`:1761-1762`) reports native R² and
RMSE; `r2_log` exists but is not the headline. For a log-distributed target the
log-space and relative metrics should carry the interpretive weight.

**Feature redundancy.** Six features spanning a rank-4 space, with
`max_features='sqrt'` giving 2 candidates per split, means importance is diluted
across collinear features — permutation importance on `chla` and `chlb_a` is
partly measuring the same information twice. Worth noting when the importances
are interpreted.

**Provenance gaps in the artifacts.** The ONNX carries a good set of metadata
(feature names, transforms, training min/max, library versions, platform,
free-text notes) — better than most. But nothing records the training file name
or hash, a git SHA, a timestamp, a model version, or the units of the inputs and
target. The dropped-row counts are printed and discarded. For EPFT-UP's
provenance goals, these are the obvious additions.

#### 7. Assessment

The code is clearly written, honestly structured, and does several things well:
one config file as the single source of truth, one shared feature-preparation
path used by training *and* inference (which eliminates the most common class of
train/serve skew), cruise-blocked validation, a real out-of-domain check, and
ONNX metadata with library versions. That is a better starting point than most
published ocean-colour algorithm code.

The weaknesses are concentrated in the part EPFT-UP cares most about: the
uncertainty is uncalibrated and inconsistent between the evaluation and
production paths, and the domain flagging — the best provenance feature in the
repo — is switched off by a filename bug (D1/D2). Those, plus the
retransformation bias, are the three things to address before building on this.

#### 8. Verification

`reports/scripts/diatom_code_checks.py` re-derives the mechanical claims above
(run with `conda run -n ocean14 python reports/scripts/diatom_code_checks.py`):
the rank-4 feature matrix and the singular covariance; the magnitude of the
retransformation bias; the metadata filename resolution that produces D1; the
ONNX-vs-JSON key sets that produce D2 and D3; and the return-arity and hardcoded
paths of D6 and D7.

## Q&A

Questions for Alison (or for J.X.P. to route), numbered for later reference.

1. **Is the domain flagging supposed to be active?** Because of the metadata
   filename mismatch (D1/D2), `run_model_inference.py` as committed emits no
   `domain_flag_*` or `mahal_dist` columns at all. Do the outputs shown in the
   preprint carry those columns — i.e. were they produced by a different (local)
   version of this script?

2. **Which uncertainty is the published one?** The evaluation figures use
   forestci's infinitesimal jackknife; the shipped inference uses per-tree
   spread. These are different quantities and will differ substantially in
   magnitude. Which one do the paper's uncertainty figures report, and should
   they agree?

3. **Has the uncertainty been checked for coverage?** Is there a figure or
   number showing what fraction of held-out observations fall within ±1σ or
   inside the q05-q95 band? If the ensemble spread is materially narrower than
   the CV residuals, the shipped `tree_std_log` will understate error.

4. **Is `diatCarb_predicted` intended as a median or a mean?** The `10**`
   back-transform with no smearing correction gives the former, but the reported
   MBE/MPE/MAPE implicitly test the latter.

5. **Where did the RF hyperparameters come from?** No search code is in the
   repository. What was the search space and CV objective — and, importantly, did
   the tuning see the hold-out cruises? That determines whether the reported test
   metrics are optimistic.

6. **Why were `t_avg` and `s_avg` dropped** despite the noted importances (0.23,
   0.09)? Was it to keep the algorithm applicable to satellite pigment retrievals
   without ancillary fields, or did they hurt generalisation across cruises?

7. **Was raw `ppc` tested as a feature?** `chla`, `chlb`, `chlc` appear both raw
   and as ratios, but `ppc` appears only as `ppc_a`. Deliberate?

8. **What are the units of `diatCarb`** — mg C m⁻³? They are not recorded in the
   ONNX metadata, the JSON, or the output CSV header.

9. **Is a single seed-42 split the basis of the reported skill,** or were
   repeated cruise-level splits / leave-one-cruise-out run and summarised
   elsewhere? With 16 cruises the test composition matters a great deal.

10. **Is the training CSV shareable?** Nothing in the repo can be run end-to-end
    without it (`training_data/` and `trained_models/` are both git-ignored, and
    no model artifact is distributed). Is there a DOI or archive for the matchup
    dataset, and would the exported ONNX + bundle be available for EPFT-UP to
    test against?

11. **Intended input domain for satellite application.** The model is trained on
    HPLC pigments. If it is to be driven by PACE pigment retrievals, is the
    retrieval uncertainty meant to be propagated, and if so where — since the
    current uncertainty assumes exact inputs?

12. **Repository URL** — the Context section of this prompt doc records
    `alisonchase/diatom_carbon_biomass`, but the local clone's remote is
    `alisonpchase/diatom-carbon-biomass-model`. Which is canonical?

## Logs

### 2026-09-15 (Explore prompt #1 — code review of the diatom carbon biomass model)

Executed the first prompt under *Prompts/Explore*: explored
`Oceanography/python/diatom-carbon-biomass-model`, wrote Report 1 above, and
raised twelve questions in the Q&A section.

What was done:

- Read all five Python modules (3,287 lines) plus README, `environment.yml`,
  `.gitignore`, `.gitattributes`, and the git history. Per the prompt's "Use
  Fable if you can", the two largest modules were analysed by Fable subagents
  (one on `train_evaluate_model.py`, one on the export/inference pair) running in
  parallel, while `model_config.py` and `prepare_features.py` were read directly.
  Every consequential claim from the subagents was then re-verified against the
  source before it entered the report; one subagent claim about which metadata
  keys live in which artifact was checked and found to need scoping, which is why
  check 5 of the verification script extracts the key sets by AST rather than by
  regex.
- Wrote `reports/scripts/diatom_code_checks.py`, which substantiates the
  mechanical findings: the rank-4 log feature matrix, the singular Mahalanobis
  covariance, the size of the retransformation bias, the metadata filename
  resolution behind D1, the ONNX-vs-JSON key sets behind D2/D3, and the
  return-arity and hardcoded paths of D6/D7. It runs against the local clone and
  needs no training data.

What I learned about the repository under review:

- It is a **pigment-to-carbon** model, not an optics model: four HPLC pigments in,
  diatom carbon out, via a 1500-tree Random Forest on log10 features and a log10
  target. There is no satellite, netCDF, or Rrs code anywhere in it — any PACE
  application depends on an upstream pigment retrieval that lives outside this
  repository, and whose error is not represented in the shipped uncertainty.
- Its best structural decision is that training and inference share one
  `prepare_features.py` code path driven by one `model_config.py`, which removes
  the usual train/serve skew. Its worst is that inference then re-reads that
  config locally instead of the metadata embedded in the model (D4), so the
  guarantee holds only if the two machines' config files match.
- The single highest-impact defect is a filename mismatch: the exporter writes
  `<stem>_metadata.json` but inference only ever looks for `metadata.json`, so
  the whole domain-flagging feature is silently inert (D1/D2) and the log-target
  flag is parsed via `bool("false") is True` (D3).
- Nothing in the repo is runnable as cloned: `trained_models/` and
  `training_data/` are git-ignored and no artifact is distributed, despite the
  README listing `trained_models/` under "Contents". There are no tests.
- Licence is CC0 1.0, unlike the BSD-3 used across this project's sibling repos —
  worth noting if any of this code is ever vendored into EPFT-UP.
- For EPFT-UP specifically: the uncertainty is uncalibrated (no coverage check
  anywhere) and uses *different estimators* in the evaluation and production
  paths, so it is not yet the "rigorous treatment of uncertainty" this project
  aims at. The domain flagging, once repaired, is the piece most worth adopting.

No git commands that change repository state were run.
