# Exploring the 2nd derivative

## Goals

Explore the 2nd derivative approach to estimating PFTs

## Context

We are going to explore the series of papers in ocean color that attempt to
infer this and that including phytoplankton functional types (PFTs) using 
the 2nd derivative of Rrs.  

I have put these papers in `context/papers`:

- catlett2017.pdf : The original paper in the Kramer+ series
- kramer2020.pddf : EOFs on PFTs
- kramer2022.pddf : The detailed manuscript describing the Kramer+ approach to PFT recovery
- El_Hournay2026.pdf : A new paper by Hournay & Kramer that further recovers PFTs

## Prompts

### Setup

1. Read the papers in `context/papers` described above.  We are going to have a discussion on the 2nd derivative approach and then write a report on it.  
Our tasks will include:
    - Read the literature
    - Discuss the maths/stats behind the methodology
    - Reproduce the Kramer+2022 paper with their data
    - Explore alternatives
Begin by asking me a series of questions in the Q&A section below.
Use Opus 5.5.  Log your work.

## Q&A

### Setup

Questions from Claude after reading the papers (2026-09-28). Answer inline under
each question (or just the ones you care about); unanswered ones get the default
noted in *[default: ...]*.

#### A. Scope of the literature and the report

1. **Paper set.** `context/papers/` also contains `lange2020.pdf` (PCR on
   *standardized* Rrs → Pro/Syn/picoeuk cell counts, AMT24). It is not a
   derivative paper. Should it be in scope as a "PCR-on-Rrs without derivatives"
   comparison point, or ignored here (it's already covered by the MOANA work)?
   *[default: include as a brief comparison point only]*
2. **Missing links in the chain.** El Hourany & Kramer (2026) lean heavily on
   Kramer et al. (2024), which uses the δRrs residual to define *optical*
   communities, and every paper leans on Kramer & Siegel (2019) for the global
   HPLC grouping. Also, Kramer 2022's key details (Table S6 A(λ),B(λ)
   coefficients; the 1st+2nd-derivative and 5/10 nm variants; Fig. S5
   coefficient spectra) live in its **Supplementary Material**. Can you add
   Kramer+2024, Kramer & Siegel 2019, and the Kramer 2022 (and ideally Catlett
   2018) supplements to `context/papers/`?
3. **Report audience & form.** Who is the report for (you / the group / an
   external paper or white paper)? Markdown in `reports/` like
   `MOANA_Claude_Report.md`? *[default: Markdown in `reports/`, written for
   ocean-optics scientists]*

#### B. Reproducing Kramer+2022

4. **Data source.** The paired dataset is on PANGAEA
   (doi:10.1594/PANGAEA.937536). Is it OK to fetch it into `epft_up/data/`
   (or wherever you keep external data), and do you already have a local copy
   or know whether it holds the 178 raw or the 145 QC'd spectra?
5. **Reference code.** El Hourany+2026 list three implementations:
   `sashajane19/Rrs_pigments` (original, likely MATLAB),
   `max-danenhower/rrs-SDP-pigments` (a Python port?), and NASA's
   `oci_sdp` notebook (suggesting this is now a PACE/OCI product). Which is
   the reference: (a) write an independent Python version from the paper and
   use the repos only to check; (b) port/wrap one of them; (c) treat the NASA
   SDP notebook as the benchmark? *[default: (a) independent re-derivation,
   checked against the repos, as we did for MOANA]*
6. **How close is "reproduced"?** Several steps are under-specified or not
   reproducible as written:
   - QC of 33 spectra was by **visual inspection** (noise at 610–660 nm).
   - The GSM-like fit (Eqs. 1–6) needs the NOMAD A(λ),B(λ) (Table S6) and
     WOA T/S for b_bw.
   - The PCR text has inconsistencies: "divided by the variance" (surely the
     std), MAD "N = 25% of 145 = 36 … the training set" (36 is the *test*
     set), and "following Eq. (8)" when Eq. (7) is meant. How the number of PCs
     is chosen (MAD on the training set? an inner CV?) is not stated.
   Is the target (i) matching Table 2 (mean R², normalized MAD) within the
   quoted SDs, (ii) matching Fig. 6/Fig. 3 qualitatively, or (iii) both plus
   coefficient spectra (Fig. S5)? And should I try to reproduce the visual QC
   with an objective rule (e.g. the Lange-style |Rrs''| threshold at 610–660 nm)?
   *[default: (i)+(ii), with an objective QC rule that recovers the 145 as
   closely as possible]*
7. **Catlett & Siegel 2018 too?** That paper works on a_ph (PnB, 491 samples),
   with 1st+2nd derivatives, a 15 nm Hamming filter, and 500 permutations.
   Reproduce it as well (if the PnB data are accessible), or treat it as
   background only? *[default: background only]*

#### C. Maths/stats: what should the discussion probe?

These are the issues I think matter most. Tell me which to emphasize, drop, or add.

8. **Skill beyond Tchla.** Accessory pigments co-vary strongly with Tchla
   (Kramer 2020 normalizes by Tchla for exactly this reason), and Kramer 2022
   models *absolute* concentrations with R² computed in *linear* space. My
   suspicion is that much of the headline R² (0.4–0.7) is Tchla retrieval in
   disguise. Should a central test be a **null model**, i.e. pigment predicted
   from retrieved Tchla alone (e.g. log–log from the GSM fit or OC4), plus
   evaluation on pigment:Tchla ratios? *[default: yes, central]*
9. **Validation design.** The "100-fold cross-validation" is really repeated
   random 75/25 subsampling across 8 cruises whose samples are spatially and
   temporally correlated, so it likely leaks. Should we add
   **leave-one-campaign-out** CV as the honest estimate? (El Hourany use
   temporal blocking; Lange use 80/20 bootstraps.) *[default: yes]*
10. **What the residual actually removes.** δRrs = Rrs − Rrs_mod, where
    Rrs_mod has a_ph = A(λ)·Tchla^B(λ) with NOMAD-mean *spectral features*. So
    δRrs'' ≈ Rrs'' − (curvature predicted by a mean-community a_ph at that
    Tchla): it is a compositional *anomaly* relative to NOMAD, not
    "phytoplankton signal minus CDOM/NAP" (the second derivative already
    kills most of the smooth adg/bbp shapes). Kramer found Rrs'+Rrs''
    performs comparably. Should I quantify how much the residual step buys
    over plain Rrs'' (and over El Hourany's spline residual M1)?
    *[default: yes]*
11. **Noise and resolution.** A finite-difference 2nd derivative amplifies
    white noise by ~√6/Δλ², and a 5 nm moving mean followed by a 1 nm
    derivative leaves highly correlated channels (the Cael+2020 information
    content argument). Do you want a formal treatment (effective degrees of
    freedom, noise propagation with PACE-OCI-like per-band uncertainties) as
    part of the "maths" section? *[default: yes, including PACE OCI noise]*
12. **n ≪ p regularization.** With 145 samples and ~290 correlated predictors,
    PCR is one of several shrinkage estimators. For the discussion, compare
    PCR with ridge/PLS (and possibly a Bayesian linear model that gives
    predictive uncertainties natively)? *[default: yes]*

#### D. Alternatives and project goals (uncertainty + provenance)

13. **Target variable.** Which should the alternatives aim at: (a) absolute
    pigment concentrations (Kramer 2022), (b) pigment:Tchla ratios, possibly
    treated as compositional data (log-ratios), (c) EOF amplitudes (Catlett),
    or (d) discrete assemblage classes (El Hourany SOM + Random Forest)?
    *[default: (b) as primary with (a) for comparison to the paper]*
14. **Alternatives you already have in mind?** Candidates I see: inverting
    a_ph(λ) first (e.g. via an IOP inversion from IOPtics) and doing the
    derivative analysis on a_ph (closer to Catlett); Gaussian-band
    decomposition of a_ph; PLS/ridge/GP; El Hourany's M1/M2/M3 source spaces;
    the MOANA approach. Should any of these be prioritized or excluded? Is the
    IOPtics machinery meant to be used here?
15. **Uncertainty deliverable.** For EPFT-UP, what counts as success on
    uncertainty: per-retrieval predictive intervals, propagation of the input
    Rrs uncertainty, and/or calibration checks (coverage) under
    leave-one-campaign-out? *[default: all three, at least for the
    best alternative]*
16. **PACE application.** Is applying the reproduced and alternative models to
    PACE OCI L2/L3 Rrs (and comparing with NASA's SDP product) in scope for
    this prompt series, or strictly in situ? *[default: in situ first, PACE
    later as its own prompt]*

## Logs

### 2026-09-28 (Setup #1: read the 2nd-derivative literature; posed Q&A)

Read all five PDFs in `context/papers/` (text extracted with `pdftotext` to the
session scratchpad, so nothing was added to the repo):

- **Catlett & Siegel 2018 (JGR, "catlett2017.pdf")**: a_ph(λ) from Plumes &
  Blooms (SBC, n=491, 350–700 nm, 1 nm). The 1st and 2nd derivatives come from
  centered finite differences after a 15 nm Hamming smoothing (chosen by
  optimizing pigment-to-2nd-derivative-peak correlations). Pigment communities
  are defined by hierarchical clustering (1−R, Ward) and EOFs of standardized
  pigments. PCR is fitted on z-scored a_ph' and a_ph'' with the coefficients
  back-transformed to A(λ), B(λ), plus 500-permutation CV. EOF modes 1, 2 and 4
  reach R² > 0.8. Fractional contributions (CHEMTAX/DP) are retrieved worse
  than pigment concentrations or EOF amplitudes. The core argument is that
  success comes from pigment covariance ("communities of absorption features").
- **Kramer, Siegel & Graff 2020 (Front. Mar. Sci.)**: HPLC only (NAAMES 1–4,
  n=229, 16 pigments normalized to Tchla). Clustering, EOFs and WGCNA network
  community detection give 5 groups (diatom, dino, hapto, green, cyano), and
  dinoflagellates separate only via the EOF mode 2 vs 4 plane. This is the
  basis for the "5 groups is the ceiling for HPLC" argument.
- **Kramer, Siegel, Maritorena & Catlett 2022 (RSE 270, 112879)**: the target
  paper. 145 QC'd global open-ocean HPLC+hyperspectral Rrs matchups (178 before
  visual QC; 8 campaigns; ±2 h), 400–700 nm, 1 nm, 5 nm moving mean. A
  GSM-like model (Gordon quadratic, a_ph = A·Tchla^B from NOMAD, S_dg from
  Rrs490/555, η from Lee 2002, and a nonlinear fit for Tchla, adg443, bbp443)
  gives the residual δRrs. The second derivative of δRrs feeds PCR for 13
  pigments, optimized on MAD, with 100 random 75/25 splits and the median of
  the 100 models used for the full reconstruction. Results: R² 0.37–0.72
  (Table 2); 1–5 nm resolution is fine, 10 nm degrades; the 5 HPLC groups are
  recovered from the modeled pigments. Data: PANGAEA 10.1594/PANGAEA.937536.
- **El Hourany & Kramer 2026 (SSRN preprint)**: 237 HPLC–Rrs stations on a
  2.5 nm grid. Four source spaces: raw Rrs, M1 (spline residual), M2
  (Gordon-type semi-analytical residual, i.e. Kramer's δ), and M3
  (Savitzky–Golay 2nd derivative of rrs). SOM classes of standardized pigment
  ratios (K = 2…50) are predicted with a Random Forest. At K=12 the balanced
  accuracy is ~0.67 for M1–M3 versus 0.42 for raw Rrs. The transformed spaces
  degrade faster under MODIS-like noise, and concatenating raw Rrs stabilizes
  them. psbO data provide a taxonomic cross-check (only 9 full triplets), and
  there is an exploratory PACE L3 application. Code: sashajane19/Rrs_pigments,
  max-danenhower/rrs-SDP-pigments, NASA oceandata `oci_sdp` notebook.
- **Lange et al. 2020 (Opt. Express)**: not listed in this prompt but present
  in the folder. PCR on spectrum-standardized Rrs (plus SST) gives Pro, Syn
  and picoeukaryote counts (AMT24). QC removed spectra with
  |Rrs''| > 2e-4 at 610–660 nm, which is a candidate objective replacement for
  Kramer's visual QC.

Repository notes: no Kramer data or code exists yet in EPFT-UP or IOPtics.
`context/papers/README` is empty. The PDF is named `El_Hournay2026.pdf`
(author is El Hour**any**). `lange2020.pdf` is duplicated in `papers/` and
`context/papers/`.

Issues flagged for discussion (see Q&A §C): (1) skill beyond Tchla; (2) random
splits leaking across correlated cruise samples; (3) what the NOMAD-based
residual actually removes; (4) noise amplification and effective degrees of
freedom; (5) PCR vs other shrinkage estimators; (6) several textual
inconsistencies in the Kramer 2022 PCR description. 16 questions posed in
Q&A → Setup; awaiting answers before any calculations.