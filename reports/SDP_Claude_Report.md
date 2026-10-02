# The 2nd-derivative (SDP) approach to phytoplankton pigments from hyperspectral Rrs

**Author:** Claude (Opus 5.5), for J. Xavier Prochaska
**Status:** in progress. §3 (Data) is written; the other sections are outlines
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

### 4.1 The residual δRrs (GSM-like model) — *Execution #2*
### 4.2 The PCR pigment model — *Execution #3*

## 5. Maths and statistics — *Execution #6*

## 6. Diagnostics — *Execution #4, #5*

### 6.1 Skill beyond Tchla (null model; pigment:Tchla ratios)
### 6.2 Leakage: random splits vs leave-one-campaign-out
### 6.3 What the residual buys (source-space comparison)

## 7. Alternatives: a learned wavelength weighting — *Execution #7*

## 8. Uncertainty — *Execution #8*

## 9. Hold-out test: EXPORTS North Atlantic 2021 — *Execution #9*

## 10. Conclusions and recommendations — *Execution #10*
