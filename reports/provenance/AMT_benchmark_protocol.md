# AMT benchmark protocol and provenance

**Purpose.** A record of what EPFT-UP has done with the pigment–hyperspectral
Rrs matchup data so far, and the evaluation protocol that makes results
comparable across groups. Anyone training a model on these data (e.g. Billy
Li's Bayesian neural networks) can then be scored with the same splits,
baselines and metrics as the EPFT-UP SDP study.

**Status:** 2026-10-24, repository commit `1b60500` + working tree.
Questions to Tom Jordan are outstanding (`requests/AMT_Rrs_HPLC_questions.md`).
Items marked *pending* are defined here but not yet produced. The prompts
"Test on new AMT data" in `claude_prompts/2nd_derivative_prompts.md` produce
them.

**Sharing.** The AMT file is for the AI_for_Ocean group only, and is shared
further only once Tom Jordan agrees (asked 2026-10-24). The Kramer deposit is
public (CC-BY-4.0).

---

## 1. Data sets

| Set | Role | N | Groups | Location | SHA-256 |
|---|---|---|---|---|---|
| **Kramer et al. (2022)**, PANGAEA doi:10.1594/PANGAEA.937536 | training set of all EPFT-UP SDP models | 145 | 8 campaigns | `$OS_COLOR/PANGAEA/Kramer2022/PANGAEA_937536.tab` | `63da2981596ef8c3d856e0ae12eab4cb715a9ec535e92d83458a7a7f5f910f3b` |
| **EXPORTS-NA 2021** (Kramer et al. 2024 test spectra + SeaBASS HPLC) | hold-out test 1 (done, report §9) | 17 | 1 campaign | `$OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA/exports_na_matchups.csv` | `90d8fd20b1231a98479c1d267f1ebc93f31479b32a7145cd2dc2755e2b9e5288` |
| **AMT_Rrs_HPLC** (compiled by Tom Jordan, received 2026-10-24) | hold-out test 2 (*pending*); later, extra training data | 210 rows → 190 unique spectra | 8 cruises | `$OS_COLOR/AMT/Rrs_HPLC/AMT_Rrs_HPLC.csv`; shared drive *AI_for_Ocean* `data/Color/AMT/Rrs_HPLC/` | `0290c9f49fae79dd52dbf920e77ff580c15af9115cd26ed734659d2c7827bca6` |

Supporting inputs, cloned at pinned commits and never copied into the repo:

| Item | Source | Pin |
|---|---|---|
| Kramer's MATLAB code, A(λ),B(λ) table, a_w table, EXPORTS-NA test spectra | `github.com/sashajane19/Rrs_pigments` | commit `b3e366225d16c4fe85c7406cd3d57e554b8cb57f` |
| Python port with Kramer's trained PCR coefficients | `github.com/max-danenhower/rrs-SDP-pigments` | commit `fb17c2f580dea9488b6ba11c690fe8c959e0a0fe` |
| A(λ),B(λ) table `aph_A_B_Coeffs_Sasha_RSE_paper.txt` | Rrs_pigments | `a6b0492df26f0cde9163f342f747853b448ae5967744f964716a2d42d6d3c64e` |
| a_w table `aw_mcf16_350_700_1nm.txt` | Rrs_pigments | `d3e8519d6485f905c3cc74afe6fc43ee67709310197010928062e2d5408f35ad` |
| EXPORTS-NA spectra `Kramer_rrs_testdata.mat` | Rrs_pigments | `b0141be8ddd3880be133b8e8847385dd71325ea711b62c6f9dac33ec5c200d2c` |
| GSM residual product for the 145 (`gsm_dRrs_insitu.npz`) | EPFT-UP `scripts/sdp/reproduce_gsm.py` | `a155b6ce40957f7276f01dc79d96ec901a13d8c5838e4eafa69cfbce4ce52e9b` |

## 2. What has been done (chronological)

| Date | Step | Script(s) | Output / record |
|---|---|---|---|
| 2026-10-02 | Fetch and verify the PANGAEA deposit and the reference repos | `scripts/sdp/fetch_kramer2022.py` | report §3.1 |
| 2026-10-02 | Ingest; Table 1 check; smoothing/rounding signature; Lange QC | `scripts/sdp/ingest_kramer2022.py` | `reports/figures/sdp/ingest_summary.json` |
| 2026-10-02 | GSM-like residual δRrs, reproduction of Kramer Fig. 4 | `scripts/sdp/reproduce_gsm.py`, `fetch_woa_ts.py` | `gsm_summary.json`; product above |
| 2026-10-02 | PCR reproduction (Table 2, Figs. 3, 6) | `scripts/sdp/reproduce_pcr.py` | `pcr_summary.json`; `products/pcr_rrsD2_1nm.npz` |
| 2026-10-02/03 | Tchla nulls, ratios, LOCO; source spaces | `diagnostics_null_loco.py`, `source_spaces.py` | `diag_*`, `source_spaces_*` |
| 2026-10-03 | Maths; learned weighting (shared W); Bayesian uncertainty | `maths_section.py`, `learned_weighting.py`, `uncertainty.py` | `maths_*`, `weighting_*`, `uncertainty_*` |
| 2026-10-03 | EXPORTS-NA matchups and hold-out 1 | `exports_na_matchups.py`, `exports_na_holdout.py` | `exports_na_*`; report §9 |
| 2026-10-04 | Report, slides, numbers check | `report_provenance.py`, `make_slides.py` | `reports/SDP_Claude_Report.md`, `reports/slides/` |
| 2026-10-24 | AMT inspection | `scripts/sdp/inspect_amt.py` | `amt_inspect_summary.json`, `amt_*.png`; report Appendix A |
| 2026-10-24 | AMT staging (laptop + shared drive), README, per-row flags | `scripts/sdp/stage_amt.py` | `README.md`, `AMT_Rrs_HPLC_flags.csv` beside the CSV |
| 2026-10-24 | Questions to the data provider | — | `requests/AMT_Rrs_HPLC_questions.md` |

All JSON outputs record the script, input checksums and run time. The code
lives in `epft_up/sdp/` (100 unit tests).

## 3. Harmonisation rules

The AMT rules are adopted defaults until Tom Jordan answers. *Pending:*
implemented in an AMT loader.

| Step | Kramer deposit (as published) | AMT (rule) |
|---|---|---|
| Spectral smoothing | already 5 nm moving mean, rounded to 1e-6 sr⁻¹ | apply a 5 nm moving mean (primary); unsmoothed as a sensitivity case |
| Grid | 400–700 nm, 1 nm | identical |
| One row per spectrum | yes | average the HPLC of the 19 shared-spectrum pairs; **drop Sample_ID 107/108** (one spectrum, Tchla 0.073 vs 0.165, likely mismatch) |
| Non-positive Rrs | none | exclude (9 rows: 34, 47, 48, 54, 60, 86, 107, 108, 198) |
| Salinity | in situ | in situ; Sample_ID 132/159 (S = 22.23) → WOA23 climatology |
| Below-detection pigments | zero | zero below the LOD proxy (next table); ½·LOD where a log is taken |
| MVchlb | measured | measured on AMT29/30; elsewhere Tchlb as a flagged proxy; MVchlb **scored** only on AMT29/30 |
| T/S for pure-water backscatter | in situ | in situ (except the salinity fix) |

**LOD proxy** [mg m⁻³] (the deposit's smallest non-zero value per pigment,
or 0.001 where it has no zeros; `epft_up/sdp/validate.py::lod_proxy`):

| Tchla | HexFuco | ButFuco | Allo | Fuco | Perid | Zea | DVchla | MVchlb | Chlc12 | Chlc3 | Neo | Viola |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 0.003 | 0.002 | 0.001 | 0.002 | 0.004 | 0.001 | 0.002 | 0.003 | 0.001 | 0.001 | 0.002 | 0.001 |

If the HPLC labs' real LODs arrive, they replace the proxy for AMT.

## 4. Splits (fold assignment)

The grouping columns already exist: `campaign` in the Kramer deposit (via
`epft_up.sdp.data.load_kramer2022().meta`) and `Cruise` in the AMT file.
*Pending:* a per-sample fold file (`sample_id, set, group, fold`) written
by the AMT loader.

| Scheme | Train | Test | Used for |
|---|---|---|---|
| R: random 75/25 × 100 (seed 1) | deposit | deposit | comparison with Kramer et al. 2022 only (leaks within cruises) |
| L8: leave-one-campaign-out | deposit, 7 campaigns | the 8th | the main skill estimate so far (report §6–§8) |
| H1: frozen hold-out | deposit (all 145) | EXPORTS-NA (17) | done (§9) |
| **H2: frozen hold-out** | deposit (all 145) | AMT, scored pooled **and per cruise** (8 groups) | *pending* (Test on new AMT data) |
| L16: leave-one-group-out | deposit + AMT | each of 8 campaigns + 8 AMT cruises | *later*, after Tom's answers |
| LA: leave-all-AMT-out | deposit | all AMT | *later* (= H2 when nothing is retrained) |

Rules:
- A shared-spectrum pair never straddles train and test.
- Hyperparameters are chosen inside the training set by an inner
  leave-group-out CV, never on the test groups.

## 5. Targets

- **13 absolute pigments** [mg m⁻³]: Tchla, HexFuco (AMT `Hex`), ButFuco
  (`But`), Allo, Fuco, Perid, Zea, DVchla, MVchlb, Chlc12, Chlc3, Neo,
  Viola.
- **12 log₁₀ pigment:Tchla ratios** for every pigment above except Tchla,
  with zeros replaced by ½·LOD before the log. This is the primary target
  for composition.

## 6. Baselines (every model is judged against these)

Fitted on the training set only:

| Baseline | Absolute pigment P | Ratio r = log₁₀(P/Tchla) |
|---|---|---|
| null_GSM | log₁₀ P = a + b·log₁₀ Tchla_GSM | r = a + b·log₁₀ Tchla_GSM |
| null_OC4 | the same with OC4v6 Tchla | the same with OC4v6 Tchla |
| const | — | the training-mean r |

- Tchla_GSM comes from the GSM-like inversion in `epft_up/sdp/gsm.py`
  (Kramer 2022 Eqs. 1–6, in-situ T/S).
- OC4v6 comes from `epft_up/sdp/validate.py::oc4v6`.
- *Pending:* both are provided per sample in the harmonised AMT table, so
  other groups need not run our code.
- For each target, a model is compared with the **best** of these baselines.

## 7. Metrics

Point predictions:
- **Absolute pigments:** R² in linear space (Kramer Table 2 convention),
  log-R², and bias, RMS and centred RMS of log₁₀(pred/obs) [dex].
- **Ratios:** log-R², bias, RMS, centred RMS [dex].
- **Skill above the null:**
  - ΔR² (random and LOCO schemes) with a paired bootstrap 95% CI over
    samples;
  - for hold-outs, ΔRMS = RMS(model) − RMS(best null) with a paired
    bootstrap CI.
  - "Beats null" means the CI excludes zero.
  - On a single hold-out campaign this CI ignores between-campaign
    variation, so **per-group results and their spread** are reported
    alongside (H2: 8 AMT cruises).

Probabilistic predictions (BNNs, the Bayesian shared W):
- coverage of the 68% and 95% central intervals;
- SD of z = (obs − μ)/σ;
- PIT histogram;
- CRPS (Gaussian closed form, or sample-based for a BNN);
- all reported pooled and per group.

**To be scored with our code**, return a CSV with
`set, sample_id, target, mu[, sigma or quantiles]`. Targets are named as
in §5, with ratios as `ratio:<pigment>`; units mg m⁻³ for absolute values
and dex for ratios.

## 8. Reference results to beat (from the EPFT-UP report)

| Scheme | Model | Result |
|---|---|---|
| L8 | Kramer PCR on δRrs″ | mean LOCO R² 0.19 (absolute), log-R² 0.17 (ratios); beats the null for 1 of 12 ratios |
| L8 | shared low-rank W (group-inner CV) | 0.39 / 0.24; beats the null for 6 of 12 ratios (Fuco, Zea, DVchla, Chlc12 strongest) |
| L8 | Bayesian shared W | 68/95% coverage 0.71/0.94 (ratios) |
| H1 | PCR / shared W on EXPORTS-NA | report §9 |

## 9. Open items

- [ ] Tom Jordan's answers: radiometry processing and resampling, the shared
  pairs, LODs, MVchlb, salinity, AMT24, sharing with Billy, credit.
- [ ] The AMT loader, the harmonised table (with Tchla_GSM, OC4 and flags),
  and the fold file: Test on new AMT data, #1.
- [ ] H2 results, report §11 and slides: Test on new AMT data, #2–#3.
- [ ] Once Tom answers: re-run H2 if the harmonisation changes; then L16/LA
  retraining; and give Billy the harmonised table.
