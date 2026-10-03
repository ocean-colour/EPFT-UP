# The 2nd-derivative (SDP) approach to phytoplankton pigments from hyperspectral Rrs

**Author:** Claude (Opus 5.5), for J. Xavier Prochaska
**Status:** in progress. §3 (Data), §4 (reproduction) and §6.1–6.2 (diagnostics) are written; the other sections are outlines
that later Execution prompts in `claude_prompts/2nd_derivative_prompts.md` fill
in.
**Scope:** a reproduction and critique of Kramer, Siegel, Maritorena & Catlett
(2022, *RSE* **270**, 112879), "SDP", and the alternatives to it, written for any
ocean-color scientist.

---

## 1. Introduction

*(Execution #10.)* Motivation; what SDP claims; what this report tests.

## 2. Literature

*(Execution #10; the reading notes are in the prompt doc's Logs, 2026-09-28/29.)*
Catlett & Siegel 2018; Kramer & Siegel 2019; Kramer, Siegel & Graff 2020;
Kramer et al. 2022; Kramer et al. 2024; El Hourany & Kramer 2026; and Lange
et al. 2020 as a brief non-derivative comparison.

## 3. Data

### 3.1 Source and provenance

All results use the PANGAEA deposit that accompanies the paper: Kramer et al.
(2021), doi:10.1594/PANGAEA.937536, CC-BY-4.0. It is the tab-delimited
"textfile" export, downloaded 2026-10-02 to
`$OS_COLOR/PANGAEA/Kramer2022/PANGAEA_937536.tab`, with
SHA-256 `63da2981…5f910f3b`, pinned in `epft_up/sdp/data.py`. Two code
repositories are used as references only and are cloned beside the data at
pinned commits. Nothing from them is copied into EPFT-UP:

| Repository | Commit | Role |
|---|---|---|
| `sashajane19/Rrs_pigments` (MATLAB, Kramer) | `b3e3662` (2025-11-03) | A(λ),B(λ) table, a_w, GSM inversion, `rrsModelTrain.m` |
| `max-danenhower/rrs-SDP-pigments` (Python port) | `fb17c2f` (2026-07-24) | trained coefficient sheets for a numerical cross-check |

Neither repository has a LICENSE file. However, the `Rrs_pigments` README says
"all code and data in this repository are freely available for use by anyone
for any and all applications". `scripts/sdp/fetch_kramer2022.py` fetches and
verifies everything. `epft_up.sdp.data.load_kramer2022()` returns the data
with a provenance record: source, checksum, the repo commits actually checked
out, loader version and load time.

### 3.2 Contents

There are **145 samples**: the quality-controlled set of the paper. Each has
Rrs at 1 nm from 400 to 700 nm (301 bands, all positive), 27 HPLC pigment
columns, in-situ temperature and salinity, time, position and sampling depth
(0–7 m). The 13 modelled pigments are Tchla, HexFuco, ButFuco, Allo, Fuco,
Perid, Zea, DVchla, MVchlb, Chlc12, Chlc3, Neo and Viola. The deposit also
carries Tchlb, Tchlc, αβ-carotene, Diadino, Diato, MVchla, chlorophyllide,
DVchlb, Lut, phaeophytin, phaeophorbide and Pras.

**The 33 spectra removed by Kramer's visual QC are not in the deposit**, so that
step cannot be reproduced and is inherited as-is.

### 3.3 Agreement with Kramer 2022 Table 1

The per-campaign counts match **exactly** in all eight campaigns, and the
Tchla range, median and mean match to the paper's quoted precision. The one
exception is the BIOSOPE mean: 0.325 here against 0.326 in the paper, a
rounding difference. The three Polarstern legs in the deposit (ANT-XXIV/4,
ANT-XXV/1, ANT-XXVI/4: 7 + 1 + 18 samples) are the paper's "ANT" (26). Over
all samples, Tchla runs from 0.019 to 4.151 mg m⁻³ with median 0.110,
against the paper's 0.019–4.15 and 0.110.

| Campaign | N | Tchla min–max | median | mean |
|---|---|---|---|---|
| ANT | 26 | 0.033–4.151 | 0.232 | 0.648 |
| NAAMES | 11 | 0.094–0.987 | 0.496 | 0.540 |
| RemSensPOC | 27 | 0.049–1.094 | 0.090 | 0.173 |
| SABOR | 9 | 0.070–1.312 | 0.252 | 0.471 |
| Tara Oceans | 16 | 0.021–0.950 | 0.168 | 0.194 |
| Tara Med | 29 | 0.026–0.170 | 0.055 | 0.064 |
| BIOSOPE | 23 | 0.019–1.471 | 0.069 | 0.325 |
| EXPORTS | 4 | 0.172–0.292 | 0.223 | 0.228 |

![Fig. 1 reproduction](figures/sdp/kramer_fig1_map.png)
![Fig. 2A reproduction](figures/sdp/kramer_fig2a_rrs.png)

*Reproductions of Kramer 2022 Fig. 1 (locations coloured by HPLC Tchla) and
Fig. 2A (measured Rrs coloured by source); compare with the paper's panels.*

### 3.4 The deposited spectra are already smoothed and rounded

Two properties of the deposit matter for everything downstream:

1. **Rrs values are rounded to 1×10⁻⁶ sr⁻¹** (six decimals; every value lies
   on that grid). Relative to the red Rrs (median ≈1.1×10⁻⁴ sr⁻¹ at
   650–700 nm), that is about 1% quantization.
2. **The spectra were smoothed by the paper's 5 nm moving mean *before*
   rounding.** A 5-point boxcar has transfer-function zeros at 0.2 and
   0.4 cycles nm⁻¹. At exactly those frequencies, the power of the deposited
   2nd difference falls to the white-noise floor expected from the rounding:
   P_data/P_quant = 1.07 and 0.96 overall, and 0.7–1.3 in every campaign. At
   the sidelobe (0.25–0.30 cycles nm⁻¹) it stays 3–4× above the floor. Only
   periods longer than ~7 nm (below ~0.15 cycles nm⁻¹) carry signal well
   above the floor: about 8× at 0.15, 350× at 0.10, 3×10⁴ at 0.05.

Consequences: the δRrs reproduction (Execution #2) must **not** apply the 5 nm
moving mean again. Kramer's MATLAB pipeline (`Kramer_hyperRrs.m`) also
applies none, consistent with this. The rounding floor is a real,
quantifiable noise term for the 1 nm 2nd derivative, with σ ≈ √6·10⁻⁶/√12 ≈
7×10⁻⁷ sr⁻¹ nm⁻² per band. It enters the noise and degrees-of-freedom
treatment in Execution #6. The script is `scripts/sdp/ingest_kramer2022.py`,
and the numbers are in `reports/figures/sdp/ingest_summary.json`.

### 3.5 Objective QC check

The Lange et al. (2020) screen rejects spectra with |Rrs''| > 2×10⁻⁴ sr⁻¹ nm⁻²
at 610–660 nm. It rejects **none** of the 145, on either the native 1 nm grid
or Lange's 2 nm grid. The largest in-band |Rrs''| is 1.1×10⁻⁵ (1 nm), 20×
below the threshold, so that threshold, tuned for HyperSAS data, cannot
discriminate on this already-smoothed deposit. That is consistent with these
being the spectra that passed Kramer's visual QC, but it is no test of that
QC.

### 3.6 Pigment zeros (below detection)

The authors set below-detection values to zero, and pigments are reported to
0.001 mg m⁻³. The zero fraction is substantial for several modelled pigments:

| Pigment | Neo | Allo | Perid | DVchla | MVchlb | Viola | Fuco | ButFuco | Chlc3 | HexFuco | Tchla, Zea, Chlc12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| zero fraction | 0.74 | 0.70 | 0.54 | 0.35 | 0.33 | 0.30 | 0.16 | 0.10 | 0.06 | 0.01 | 0 |

Two implications:
- For log-ratio targets (Q&A #13, #20), the ½·LOD replacement will dominate
  Neo, Allo and Perid, so their sensitivity analysis is essential.
- Kramer's stated exclusion rule ("below detection in >75% of samples") does
  not reproduce exactly from the deposit. DVchlb (0.68 zero) and Lut (0.68)
  were excluded while Neo (0.74) was kept, so the rule was presumably applied
  to the pre-QC 178 samples, or with "not measured" counted.

## 4. Reproduction of Kramer et al. (2022)

### 4.1 The residual δRrs (GSM-like model)

**Implementation.** `epft_up/sdp/gsm.py` is written from Kramer 2022 §2.3
(Eqs. 1–6) and the primary sources it cites, not ported:

- Gordon et al. (1988): rrs = 0.0949u + 0.0794u², with u = b_b/(a+b_b).
- The Lee et al. (2002) conversion between Rrs and rrs.
- Pure-seawater scattering from Zhang et al. (2009) at each sample's in-situ
  T and S, with b_bw = b_sw/2.
- a_ph = A(λ)·Tchla^B(λ), from the reference A,B table.
- a_dg = a_dg(443)·exp(−S_dg(λ−443)), with S_dg from Carder et al. (1999).
- b_bp = b_bp(443)·(443/λ)^η, with η from Lee et al. (2002).
- Three free parameters, fitted by unweighted least squares in rrs.

The deposited spectra are used as-is, with **no further smoothing** (§3.4).
The script is `scripts/sdp/reproduce_gsm.py`, and the numbers are in
`reports/figures/sdp/gsm_summary.json`. The δRrs product for §4.2 is saved
with its provenance at
`$OS_COLOR/PANGAEA/Kramer2022/products/gsm_dRrs_insitu.npz`.

**Where the paper, the MATLAB original and the Python port differ:**

| Item | Paper text | `Rrs_pigments` (MATLAB) | Python port | Used here |
|---|---|---|---|---|
| S_dg (Eq. 5) | S = −0.01447 + 0.00033·Rrs490/Rrs555 inside exp(S(λ−443)) | exp(−(0.01447 + 0.00033 r)(λ−443)) | as MATLAB | MATLAB = Carder (1999). The printed sign is a typo; used literally, GSM-Tchla R² collapses to 0.07 |
| η band ratio | "rrs(490)/rrs(555)" | rrs(440)/rrs(555) (Lee 2002) | **Rrs**(440)/Rrs(555) | rrs(440)/rrs(555); the 490 variant gives R² 0.62 |
| a_w | Mason et al. (2016) | `aw_mcf16_350_700_1nm` | same file | same file (identical to `asw_all.mat`) |
| b_bw T/S | WOA13 ¼° climatology | user-supplied | user-supplied | in-situ (deposit); WOA23 as a check |
| Optimizer | "non-linear fit" | `fminsearch` from (0.15, 0.01, 0.0029), TolX = TolFun = 1e-9, ≤2000 iterations, unbounded; NaN unless exitflag = 1 | `scipy.fmin`, tolerances 1e-6, no convergence check | scipy Nelder–Mead with MATLAB's settings and initial simplex; trials with Tchla ≤ 0 get infinite cost |
| Smoothing | 5 nm moving mean, trim 4 nm | none | none | none: already applied in the deposit (§3.4) |

**Fig. 4 is reproduced** (log₁₀ statistics, OLS of log modelled on log
HPLC):

| | paper | here, N=145 | here, without one fit (N=144) |
|---|---|---|---|
| OC4v6 | y = 0.87x − 0.14, R² = 0.75 | y = 0.873x − 0.138, R² = 0.746 | — |
| GSM-like | y = 0.96x − 0.093, R² = 0.86 | y = 0.942x − 0.126, R² = 0.717 | **y = 0.961x − 0.091, R² = 0.864** |

The OC4 match shows that the data, the band sampling and the log-space
statistics are the paper's. The GSM difference is **one spectrum**: SABOR,
2014-07-31 (index 69), HPLC Tchla 0.252. Its fit is degenerate: Tchla →
3.8×10⁻⁴ while a_dg(443) = 0.12 m⁻¹ absorbs the blue. With it removed, the
statistics are the paper's to the quoted digits. The paper's Fig. 4B has no
point anywhere near log₁₀(GSM) ≈ −3.4, so in Kramer's run this sample either
failed to converge (and was set to NaN by `gsm_invert.m`) or converged
elsewhere. We can't run MATLAB here, but the mechanism is plausible:

- 22 of this fit's 366 Nelder–Mead trial points have Tchla ≤ 0, where
  `Tchla^B` is complex in MATLAB.
- MATLAB's simplex sorting and comparisons on complex costs differ from our
  infinite-cost guard and from the port's NaNs.

The paper's other visible outlier is (log₁₀ HPLC, log₁₀ GSM) ≈ (−1.2, −2.3),
which is our Tara Med sample with HPLC 0.068 and GSM 0.005. It is reproduced.

![Fig. 4 reproduction](figures/sdp/kramer_fig4_chl.png)

*Kramer Fig. 4: (A) OC4v6, (B) the GSM-like fit. The circled point is the
degenerate SABOR fit discussed above.*

Two further features matter for downstream use:
- Five ANT (Polarstern) spectra are fitted with **b_bp(443) < 0**. The MATLAB
  fit is unbounded, so we keep this faithfully. A positivity-bounded fit
  moves those five to b_bp ≈ 0 and Tchla 0.4–3, and lowers the Fig. 4 R²
  to 0.66. It changes no other spectrum: δRrs differs by 6×10⁻⁸ relative
  RMS elsewhere.
- Per-campaign log₁₀(GSM/HPLC) biases range from −0.16 (Tara Med, NAAMES)
  to +0.04 (RemSensPOC), with SABOR the noisiest (SD 0.9, driven by
  index 69).

**Fig. 2B–C is reproduced** qualitatively. The modelled spectra track the
measured ones. The residuals have the paper's range, about ±3×10⁻⁴ sr⁻¹
with a larger excursion near 400 nm. They also show its features: the
black EXPORTS dip near 430 nm, the structure at 450–520 nm, and the ANT
bump near 680 nm, which is chlorophyll fluorescence the model lacks.

![Fig. 2B–C reproduction](figures/sdp/kramer_fig2bc_model_residual.png)

**Cross-checks and sensitivities of δRrs.** All are relative RMS over every
sample and band, against an RMS δRrs of ≈7×10⁻⁵ sr⁻¹:

| Change | Effect on δRrs | Comment |
|---|---|---|
| Python port, with the port's η convention | 4×10⁻⁵ (max 4.5×10⁻⁸ sr⁻¹) | agreement to below the deposit's rounding (2.9×10⁻⁷); b_sw identical to 8×10⁻¹⁶ |
| η from Rrs (port) vs rrs (MATLAB) | 8×10⁻⁴ (max 7×10⁻⁷) | negligible |
| WOA23 monthly vs in-situ T/S (Q&A #22) | 1.2×10⁻³ (max 1.2×10⁻⁶) | negligible. In-situ − WOA: T = +0.7 ± 1.4 °C, S = 0.0 ± 0.4. Tchla changes ≤0.4% except the degenerate SABOR fit (15%) |
| Bounded optimizer | 0.58 | entirely the five negative-b_bp ANT spectra |
| A second 5 nm moving mean | 0.04 (max 3×10⁻⁵) | 10–100× the T/S or η effects. Re-smoothing would be a real error, and it grows under differentiation (§4.2) |

![δRrs sensitivity](figures/sdp/gsm_sensitivity.png)

*Per-wavelength RMS of each perturbation compared with the RMS of δRrs. The
WOA, η and port differences sit at or below the deposit's rounding level.*

The WOA comparison uses WOA23, since WOA13 is no longer served: the ¼°
monthly objectively analysed surface fields at each sample's month and
nearest ocean cell (`scripts/sdp/fetch_woa_ts.py`; values in
`$OS_COLOR/PANGAEA/Kramer2022/woa23_surface_ts.csv`). The largest T
difference, +8.8 °C, is a RemSensPOC sample at 55°N, 49°W in August (in-situ
19.0 °C vs WOA 10.2 °C). That in-situ value looks suspect but has no effect
on δRrs.

**Open for §4.2:** keep the degenerate SABOR spectrum (faithful to our run)
or drop it (probably faithful to the paper's). Both will be reported.
### 4.2 The PCR pigment model

**Implementation.** `epft_up/sdp/models.py::train_pcr_kramer` re-implements
`rrsModelTrain.m` from its source, with Kramer's settings from
`Kramer_Rrs_pigments.m`:

- Predictors are `diff(δRrs, 2)`: plain second differences on the 1 nm
  grid, 299 bands at 401–699 nm.
- 100 permutations of a random 75/25 split.
- Inside each permutation, a 5-fold CV chooses the number of PCs
  (≤ 30, by validation MAE).
- Fold coefficients are averaged and back-transformed to A_m(λ) and C_m.
- Outputs are clipped at ≥ 0.

Two quirks of the original are kept deliberately:
- Each inner validation fold is z-scored with **its own** mean and SD.
- The fold coefficients are averaged in standardized units before the
  conversion to raw units.

Within a fold the PC scores are orthogonal and centred, so the OLS
coefficients of the first l PCs do not depend on l. All 30 candidate models
therefore come from one fit; a unit test checks this against explicit
refits. The script is `scripts/sdp/reproduce_pcr.py`; the numbers are in
`reports/figures/sdp/pcr_summary.json` and `kramer_table2.csv`, and the
model ensembles are in `$OS_COLOR/PANGAEA/Kramer2022/products/pcr_rrsD2_1nm.npz`.
MATLAB's random splits cannot be reproduced across languages, so agreement
is statistical (seed 1).

**What Table 2's "normalized MAD" is.** It is the validation MAE divided by
the mean *modelled* value of the validation set: Tchla 0.507 here vs 0.498
in the paper. Dividing by the mean observed value instead gives 0.516 for
Tchla and 1.03 for Fuco (paper 0.84), so that is not it.

**Table 2 is reproduced.** Every pigment's mean R² is within one quoted SD
of the paper's. The table uses all 145 spectra; the 144-spectrum run
(without the degenerate SABOR GSM fit of §4.1) differs by ≤ 0.03 in R².

| Pigment | R² here | R² paper | Δ / SD_paper | MADn here | MADn paper | median #PCs |
|---|---|---|---|---|---|---|
| Tchla | 0.73 ± 0.13 | 0.72 ± 0.15 | +0.1 | 0.51 ± 0.15 | 0.50 ± 0.13 | 22 |
| Chlc12 | 0.67 ± 0.13 | 0.70 ± 0.13 | −0.2 | 0.71 | 0.70 | 20 |
| Chlc3 | 0.62 ± 0.12 | 0.68 ± 0.13 | −0.4 | 0.70 | 0.64 | 17 |
| Fuco | 0.63 ± 0.13 | 0.65 ± 0.15 | −0.2 | 0.85 | 0.84 | 18 |
| ButFuco | 0.51 ± 0.15 | 0.62 ± 0.16 | −0.7 | 0.64 | 0.59 | 13 |
| DVchla | 0.49 ± 0.13 | 0.55 ± 0.12 | −0.5 | 0.64 | 0.59 | 25 |
| HexFuco | 0.48 ± 0.13 | 0.54 ± 0.16 | −0.4 | 0.73 | 0.69 | 13 |
| Perid | 0.55 ± 0.12 | 0.49 ± 0.13 | +0.4 | 0.74 | 0.78 | 18 |
| MVchlb | 0.41 ± 0.16 | 0.42 ± 0.19 | −0.1 | 0.98 | 0.98 | 14 |
| Neo | 0.40 ± 0.16 | 0.42 ± 0.21 | −0.1 | 1.11 | 1.13 | 15 |
| Allo | 0.38 ± 0.15 | 0.40 ± 0.19 | −0.1 | 1.28 | 1.22 | 12 |
| Viola | 0.35 ± 0.17 | 0.38 ± 0.18 | −0.2 | 1.14 | 1.10 | 13 |
| Zea | 0.37 ± 0.11 | 0.37 ± 0.10 | 0.0 | 0.50 | 0.47 | 21 |

The ±σ is the spread over the 100 validation splits. The Monte-Carlo error
of each mean is ≈ σ/10 ≈ 0.015. Our R² is slightly lower on average (mean
Δ ≈ −0.2 SD), most for ButFuco and DVchla, but nowhere outside the quoted
spread.

**Kramer's own trained coefficients reproduce on our δRrs''.** The Python
port ships `original_a_coefs.xlsx` / `original_c_coefs.xlsx`: 100 trained
models per pigment on the same 401–699 nm grid. Applied to *our* δRrs'',
with no retraining, they give full-reconstruction R² within ±0.05 of our own
models for every pigment (Tchla 0.79 vs 0.82, Fuco 0.79 vs 0.80, Chlc12 0.81
vs 0.80). The median coefficient spectra correlate at r = 0.85–0.90 for 10
of 13 pigments, and 0.78 for Perid. Zea and DVchla correlate at only ≈0.5,
yet still predict as well: these pigments' coefficients are poorly
constrained, so many weight spectra fit equally well. Median intercepts
agree to within ≈0.01–0.06 mg m⁻³. This is an independent confirmation that
our δRrs and PCR are Kramer's, not just statistically similar.

![PCR coefficient spectra](figures/sdp/pcr_coefficients.png)

*Median and interquartile range of our 100 A(λ) for Tchla, Fuco and HexFuco,
compared with the median of Kramer's original coefficients. Note that the
largest weights sit at 520–700 nm, where δRrs'' is smallest. That is a
consequence of z-scoring each band before the PCA (taken up in §5).*

**Fig. 6 (median of the 100 models applied to every spectrum).** The
paper's panels are in **log₁₀** space, so samples that are zero on either
axis drop out. Modelled values below the detection limit are set to zero,
as in the paper. The deposit has no method LODs, so I use a proxy:
- for pigments with below-LOD zeros, the smallest non-zero deposited value
  (0.001–0.004 mg m⁻³);
- for Tchla, Zea and Chlc12, which are never below LOD here, the 0.001
  reporting resolution.

| Pigment | here (N=144): slope, intercept, R² (N in plot) | paper: slope, intercept, R² |
|---|---|---|
| Tchla | 0.88, −0.02, 0.71 (127) | 0.94, −0.004, 0.73 |
| Fuco | 0.65, −0.20, 0.71 (89) | 0.80, −0.08, 0.60 |
| Perid | 0.55, −0.83, 0.60 (56) | 0.78, −0.46, 0.51 |
| HexFuco | 0.69, −0.27, 0.64 (125) | 0.76, −0.20, 0.64 |
| MVchlb | 0.62, −0.39, 0.63 (78) | 0.74, −0.26, 0.51 |
| Zea | 0.51, −0.63, 0.50 (144) | 0.53, −0.59, 0.55 |

R² agrees, but **our slopes are lower by 0.02–0.23.** This is not our
training: Kramer's original coefficients on our δRrs'' give the same low
slopes (Tchla 0.84, Fuco 0.61, HexFuco 0.72). Nor is it the regression
convention: ordinary y-on-x, reduced major axis and inverted x-on-y each
fail to match all six panels. One visible difference is that 17 low-Tchla
samples have a modelled Tchla ≤ 0 (the linear model goes negative) and drop
out of the log plot. The paper's panel A seems to keep more low points. With
no quoted uncertainty for Fig. 6, I record this as an **open discrepancy**.
It does not affect Table 2, which agrees.

![Fig. 6 reproduction](figures/sdp/kramer_fig6_pigments.png)

**Fig. 3 (clustering of the 12 accessory-pigment:Tchla ratios; Ward linkage
on 1 − R).** The *measured* dendrogram (cut at 0.65) recovers exactly the
paper's five groups:
- haptophytes {HexFuco, ButFuco};
- green algae {Allo, MVchlb, Neo, Viola};
- diatoms {Fuco, Chlc12, Chlc3};
- dinoflagellates {Perid};
- cyanobacteria {Zea, DVchla}.

The *modelled* dendrogram (cut at 0.80) is fragile:
- The green-algal and diatom groups survive.
- With the detection-limit handling above, {HexFuco, Perid, Zea, DVchla}
  merge and ButFuco stands alone.
- With a coarser proxy (Tchla LOD = 0.019, i.e. 17 more zeroed Tchla
  values), five groups appeared, with Chlc3 moving to Perid.

Modelled ratios exist only where modelled Tchla > 0 (127 of 144 samples), and
clipped zeros dominate the low-concentration ratios. The paper's statement
that "the same five groups emerge" from the modelled pigments therefore
holds only loosely, and it depends on how zeros are handled.

![Fig. 3 reproduction](figures/sdp/kramer_fig3_dendrograms.png)

**Supplementary variants** (mean R² over 100 splits, N = 144):

| Pigment | δRrs'' 1 nm (paper) | Rrs' + Rrs'' 1 nm | δRrs'' 5 nm | δRrs'' 10 nm |
|---|---|---|---|---|
| Tchla | 0.73 | 0.76 | 0.72 | 0.71 |
| Fuco | 0.63 | 0.72 | 0.67 | 0.71 |
| HexFuco | 0.50 | 0.59 | 0.62 | 0.62 |
| ButFuco | 0.53 | 0.61 | 0.64 | 0.63 |
| Chlc12 | 0.66 | 0.75 | 0.69 | 0.71 |
| Chlc3 | 0.64 | 0.71 | 0.70 | 0.70 |
| DVchla | 0.46 | 0.54 | 0.54 | 0.54 |
| Perid | 0.55 | 0.54 | 0.52 | 0.45 |
| Zea | 0.36 | 0.36 | 0.43 | 0.42 |
| MVchlb | 0.42 | 0.43 | 0.47 | 0.50 |
| Neo | 0.42 | 0.43 | 0.48 | 0.55 |
| Allo | 0.39 | 0.43 | 0.47 | 0.47 |
| Viola | 0.37 | 0.39 | 0.45 | 0.50 |

Two results:
- **Rrs' + Rrs'' of the *measured* reflectance (no GSM residual) does as
  well as or better than δRrs''** for 11 of 13 pigments. That is consistent
  with the paper's "comparable results", and it is the first hint that the
  residual step adds little; §6.3 tests this properly.
- **Coarser sampling does not hurt here.** Here, "5 nm" and "10 nm" mean
  every 5th or 10th band of the already-smoothed δRrs, then `diff(·, 2)`.
  At 5 nm, skill is flat or higher, in line with the paper. At 10 nm,
  skill is *not* "notably worse" as the paper reports: only Perid drops,
  0.55 → 0.45. The paper does not say how it degraded the resolution (the
  supplement tables S2–S3 are not available to us), so this may be a recipe
  difference. With 29 rather than 299 predictors, the PCR has far fewer
  noisy dimensions to fit, which helps the weaker pigments. Either way, the
  claim that hyperspectral 1 nm resolution is needed is not supported by
  this dataset under this recipe. §5 treats the effective degrees of
  freedom.

## 5. Maths and statistics — *Execution #6*

## 6. Diagnostics — *Execution #4, #5*

The code is in `epft_up/sdp/validate.py` and `scripts/sdp/diagnostics_null_loco.py`.
The numbers are in `reports/figures/sdp/diag_summary.json` and
`diag_skill_table.csv`. Everything here uses all 145 spectra, as the paper
does. Dropping the degenerate SABOR fit changes no conclusion; those numbers
are in the JSON.

**Design.** Every model is trained and scored on *identical* splits, so the
comparisons are paired. There are two validation schemes:
- Kramer's 100 random 75/25 splits (the same split sizes; seed 1);
- **leave-one-campaign-out (LOCO)**: each of the 8 campaigns held out in
  turn, with one prediction per sample pooled over the folds.

The models compared:

| Model | Predictor | What it measures |
|---|---|---|
| PCR | δRrs'' (Kramer 2022) | the paper's method |
| null_GSM | log₁₀ P = a + b log₁₀ Tchla_GSM, fitted per split | what a hyperspectral Tchla retrieval alone gives |
| null_OC4 | the same with OC4v6 band-ratio Tchla | what a *multispectral* Tchla gives |
| oracle | the same with HPLC Tchla | the ceiling of Tchla covariance (not achievable) |
| const | training mean | no skill (ratios) |

The two null models have only 2 parameters, fitted on the training split.
Zeros are replaced by ½·LOD wherever a log is taken, through
`validate.replace_zeros`, with the LOD proxy of §4.2.

### 6.1 Skill beyond Tchla

**Absolute concentrations: SDP has no skill beyond Tchla, except for the
cyanobacterial pigments, and that exception disappears out of campaign.**
R² below is Table 2's linear-space statistic. "PCR > null_GSM" is the
fraction of the 100 paired splits in which PCR has the higher R².

| Pigment | PCR | null_GSM | null_OC4 | oracle | PCR > null_GSM |
|---|---|---|---|---|---|
| Tchla | 0.72 | **0.81** | 0.65 | — | 22% |
| Chlc12 | 0.64 | **0.78** | 0.55 | 0.92 | 18% |
| Chlc3 | 0.62 | **0.73** | 0.47 | 0.82 | 31% |
| Fuco | 0.61 | **0.72** | 0.57 | 0.81 | 21% |
| MVchlb | 0.44 | **0.62** | 0.45 | 0.74 | 18% |
| HexFuco | 0.48 | **0.61** | 0.39 | 0.73 | 30% |
| Perid | 0.55 | 0.56 | **0.60** | 0.53 | 55% |
| ButFuco | 0.50 | **0.52** | 0.34 | 0.60 | 53% |
| Neo | 0.43 | **0.52** | 0.40 | 0.66 | 34% |
| Viola | 0.38 | **0.52** | 0.42 | 0.68 | 28% |
| Allo | 0.38 | **0.47** | 0.38 | 0.75 | 33% |
| Zea | **0.37** | 0.20 | 0.08 | 0.21 | 87% |
| DVchla | **0.48** | 0.06 | 0.02 | 0.06 | 100% |

What the table shows:
- A two-parameter power law in the **GSM-retrieved Tchla beats PCR for 11
  of 13 pigments**, Tchla itself included. That GSM Tchla is a by-product of
  the very residual step SDP uses.
- **Even the multispectral OC4 Tchla matches or beats PCR for 6 of 13**
  pigments.
- PCR beats the Tchla nulls only for **DVchla and Zea**, the two pigments
  whose concentration does *not* scale with Tchla (Prochlorococcus and
  Synechococcus dominate the oligotrophic end, so the null R² is ≈0).
- The oracle sits above PCR for every pigment except Perid and Zea. Most of
  the information in the HPLC table is Tchla covariance, and SDP recovers
  less of it than a direct Tchla retrieval does.
- Log-space scoring gives the same picture. There, PCR is further penalized
  by its clipped zeros: Tchla log-R² is 0.43 for PCR vs 0.72 for null_GSM.

**Composition (log₁₀ pigment:Tchla).** The table compares three versions of
the ratio:
- PCR trained directly on the log ratio;
- the ratio of PCR's *absolute* predictions ("derived"), which is what a
  user of SDP pigment products would compute;
- the better of the two Tchla nulls (a ratio-on-Tchla power law).

| Ratio | random: PCR_log | derived | best null | PCR > null | LOCO: PCR_log | best null | ΔR² 95% CI |
|---|---|---|---|---|---|---|---|
| DVchla:Tchla | 0.69 | 0.29 | 0.42 | 100% | 0.36 | 0.32 | [−0.05, 0.12] |
| Fuco:Tchla | 0.61 | 0.03 | 0.38 | 100% | 0.40 | 0.31 | [−0.01, 0.17] |
| Zea:Tchla | 0.57 | 0.22 | 0.42 | 90% | 0.34 | 0.35 | [−0.10, 0.09] |
| Chlc3:Tchla | 0.55 | 0.03 | 0.34 | 98% | 0.17 | 0.20 | [−0.14, 0.07] |
| Chlc12:Tchla | 0.44 | 0.03 | 0.34 | 94% | 0.28 | 0.27 | [−0.05, 0.08] |
| Neo:Tchla | 0.38 | 0.09 | 0.20 | 91% | 0.20 | 0.04 | [0.02, 0.28]* |
| Allo:Tchla | 0.31 | 0.02 | 0.05 | 99% | 0.08 | 0.20 | [−0.24, 0.02] |
| Perid:Tchla | 0.29 | 0.08 | 0.19 | 83% | 0.09 | 0.08 | [−0.08, 0.08] |
| MVchlb:Tchla | 0.24 | 0.02 | 0.17 | 74% | 0.03 | 0.04 | [−0.06, 0.03] |
| HexFuco:Tchla | 0.18 | 0.02 | 0.06 | 89% | 0.03 | 0.21 | [−0.30, −0.06] |
| ButFuco:Tchla | 0.18 | 0.02 | 0.06 | 80% | 0.00 | 0.08 | [−0.21, 0.01] |
| Viola:Tchla | 0.17 | 0.02 | 0.04 | 92% | 0.01 | 0.29 | [−0.44, −0.13] |

\* Neo is 74% below detection, and its result is an artifact of the zero
replacement. Under LOCO its PCR log-ratio R² is 0.03, 0.20 and 0.47 for
zeros at 0.1, 0.5 and 1.0 × LOD. The ratios with few zeros are insensitive
to that choice: Fuco 0.39/0.42/0.41, Zea 0.33 throughout, Chlc12 0.28
throughout.

Three conclusions:
1. **The derived ratios are empty.** The pigment ratios implied by SDP's
   absolute products carry essentially no compositional information:
   R² 0.02–0.09, rising to 0.2–0.3 only for Zea and DVchla. SDP's absolute
   pigments are, to first order, Tchla times a constant.
2. **Trained on the ratio itself, the spectra do carry composition in
   random splits.** Fuco:Tchla, Chlc3:Tchla, DVchla:Tchla and Zea:Tchla
   reach R² 0.55–0.69 and beat the nulls in 90–100% of splits. This is
   where δRrs'' has genuine information beyond Tchla.
3. **That information does not transfer between campaigns (§6.2).**

### 6.2 Leakage: random splits vs leave-one-campaign-out

Kramer's random splits put samples from the same cruise in both training and
validation. Holding out a whole campaign shows how much of the skill is
cruise-specific.

**Absolute concentrations, LOCO pooled R²:**

| | Tchla | Fuco | Chlc12 | Chlc3 | HexFuco | MVchlb | Perid | DVchla | Zea |
|---|---|---|---|---|---|---|---|---|---|
| PCR, random | 0.72 | 0.61 | 0.64 | 0.62 | 0.48 | 0.44 | 0.55 | 0.48 | 0.37 |
| **PCR, LOCO** | **0.44** | **0.32** | **0.29** | **0.24** | **0.12** | **0.15** | **0.45** | **0.04** | **0.01** |
| null_GSM, LOCO | 0.66 | 0.56 | 0.54 | 0.58 | 0.39 | 0.37 | 0.48 | 0.21 | 0.02 |
| null_OC4, LOCO | 0.46 | 0.39 | 0.28 | 0.20 | 0.16 | 0.13 | 0.50 | 0.12 | 0.02 |

- PCR loses **roughly half its R²** when a campaign is held out: Tchla 0.72
  → 0.44, HexFuco 0.48 → 0.12. Its advantage on the cyanobacterial
  pigments vanishes: DVchla goes from 0.48 to 0.04.
- The Tchla nulls lose much less (null_GSM Tchla 0.81 → 0.66).
- Under LOCO, PCR is **significantly worse than null_GSM** (the paired
  bootstrap 95% interval of ΔR² lies below 0) for Tchla, HexFuco, Allo,
  DVchla, MVchlb, Chlc12 and Chlc3. It is better for none of the 13
  pigments.

![Absolute-concentration skill](figures/sdp/diag_skill_absolute.png)

*R² of PCR and the Tchla-only null models for each pigment: random splits
(left) vs LOCO (right), linear (top) vs log (bottom) scoring. The open
diamonds (oracle, HPLC Tchla) show how much of each pigment is explained by
Tchla covariance alone.*

**Composition under LOCO.**
- **No ratio retains significant skill beyond the Tchla nulls** when a
  campaign is held out. Neo appears to, but only through the zero
  replacement (see above).
- The best candidates are **Fuco:Tchla** (0.40 vs 0.31; CI [−0.01, 0.17]),
  then DVchla:Tchla and Zea:Tchla (≈ the null).
- For HexFuco, ButFuco, MVchlb and Viola, PCR's held-out log-ratio RMS error
  is **no better than predicting the training-mean ratio**: 0.28 vs 0.25,
  0.31 vs 0.30, 0.35 vs 0.35, and 0.35 vs 0.33 dex. Only Fuco (0.29 vs
  0.37), Zea (0.44 vs 0.54) and DVchla (0.77 vs 0.98) clearly beat a constant
  ratio out of campaign.

![Ratio skill](figures/sdp/diag_skill_ratio.png)

**Why LOCO hurts.** Held-out-campaign predictions are systematically offset
by campaign. When it is held out, BIOSOPE (South Pacific gyre, Tchla ≈ 0.02
mg m⁻³) is predicted at 0.2–0.5 mg m⁻³. Many RemSensPOC and Tara Med
spectra are predicted at ≤ 0. The PCR weights the fine spectral structure
of δRrs'', and part of that structure is specific to the instrument, the
processing or the region of each campaign (§3.4, §4.1): RAMSES for ANT,
HyperPro for the others, and different processing chains. A two-parameter
Tchla model does not see that structure and transfers better.

![LOCO Tchla](figures/sdp/diag_loco_tchla.png)

**Bottom line for §6.1–6.2.** On this dataset:
1. SDP's absolute pigment concentrations are matched or beaten by a
   Tchla-only power law, using the GSM Tchla that SDP already computes,
   under both random and campaign-held-out validation.
2. δRrs'' does carry compositional information within campaigns, most
   clearly for Fuco:Tchla and the cyanobacterial ratios.
3. **None of it is demonstrated to transfer to an unseen campaign.** At most
   a weak Fuco:Tchla signal survives (not significant at 95%).

Kramer's Table 2 statistics are real, but they measure mostly Tchla and
within-cruise similarity. The question for §6.3 and §7 is whether a
different predictor space or a regularized, noise-aware weighting can make
the compositional signal transfer.
### 6.3 What the residual buys (source-space comparison)

## 7. Alternatives: a learned wavelength weighting — *Execution #7*

## 8. Uncertainty — *Execution #8*

## 9. Hold-out test: EXPORTS North Atlantic 2021 — *Execution #9*

## 10. Conclusions and recommendations — *Execution #10*
