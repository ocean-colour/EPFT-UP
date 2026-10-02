# Exploring the 2nd derivative

## Goals

Explore the 2nd derivative approach to estimating PFTs

## Context

We are going to explore the series of papers in ocean color that attempt to
infer this and that including phytoplankton functional types (PFTs) using 
the 2nd derivative of Rrs.  

I have put these papers in `context/papers`:

- catlett2018.pdf : The original paper in the Kramer+ series
- Kramer2019.pdf : 
- kramer2020.pdf : EOFs on PFTs
- kramer2022.pdf : The detailed manuscript describing the Kramer+ approach to PFT recovery
- Kramer_etal_2024.pdf : Kramer's latest
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

2. I have answered your questions in the Q&A section below.  Please review them and let me know if you have any additional questions or comments.
Use Opus 5.5.  Log your work.

3. I have answered your second round of questions in the Q&A section below.  Please review them and let me know if you have any additional questions or comments.
Use Opus 5.5.  Log your work.

4. I have answered your third round of questions in the Q&A section below.  Please review them and react accordingly.  Use Opus 5.5.  Log your work.

### Execution

Conventions for every prompt below (from the Setup Q&A): code lives in
`epft_up/sdp/` (tests in `epft_up/tests/`, run scripts in `scripts/sdp/`,
figures in `reports/figures/sdp/`); external data and reference code go under
`$OS_COLOR/PANGAEA/Kramer2022/` and are never committed; every calculation is a
script on disk run with `conda run -n ocean14`; every product records its
provenance (input files + checksums, reference-repo commit SHAs, code version,
parameters, RNG seeds). Each prompt ends by appending its findings to the
relevant section of `reports/SDP_Claude_Report.md` (created in #1) and a dated
entry under Logs. Stop and ask if a result contradicts the paper by more than
the quoted uncertainty and the cause isn't obvious.

1. **Data ingest and provenance.** Read this file. Fetch PANGAEA
   10.1594/PANGAEA.937536 into `$OS_COLOR/PANGAEA/Kramer2022/` and clone
   `sashajane19/Rrs_pigments` and `max-danenhower/rrs-SDP-pigments` into
   `.../Kramer2022/ref/` (record commit SHAs; vendor nothing). Write
   `epft_up/sdp/data.py` to load the samples into a tidy structure (Rrs on the
   1 nm 400–700 grid, the 13 pigments plus the extra pigments, T/S, campaign,
   lat/lon/time) with a provenance record, plus tests. Verify N=145 and the
   per-campaign counts against Kramer 2022 Table 1. Apply the Lange-style
   |Rrs''| check at 610–660 nm and report how many of the 145 it would reject.
   Reproduce Kramer Figs. 1 and 2A as sanity plots. Start
   `reports/SDP_Claude_Report.md` with an outline (Introduction, Literature,
   Data, Reproduction, Maths/Stats, Diagnostics, Alternatives, Uncertainty,
   Hold-out test, Conclusions).
   Use Opus 5.5. Log your work.

2. **Reproduce δRrs (the GSM-like residual).** Read this file. Implement
   `epft_up/sdp/gsm.py` from Kramer 2022 Eqs. 1–6 *independently* (Gordon
   quadratic; a_w from Mason+2016; b_bw from Zhang+2009 using in-situ T/S;
   a_ph = A(λ)·Tchla^B(λ) from the reference A,B table; S_dg and η from the
   stated band-ratio relations; nonlinear fit for Tchla, a_dg(443), b_bp(443)).
   Include the 5 nm moving-mean smoothing and trimming as in the paper. Then
   cross-check against the Python port and `Kramer_hyperRrs.m`, documenting
   any differences in equations, constants or fit bounds. Deliverables:
   Fig. 4B (GSM Tchla vs HPLC; target R²≈0.86, slope≈0.96) and Fig. 2B–C
   (Rrs_mod, δRrs), plus the WOA-vs-in-situ T/S effect on δRrs (Q22).
   Tests for the forward model and the fit.
   Use Opus 5.5. Log your work.

3. **Reproduce the PCR model (Kramer 2022 §2.5, Table 2, Figs. 3 & 6).**
   Read this file. Implement `epft_up/sdp/spectral.py` (2nd-order finite
   differences, smoothing options) and the PCR in `epft_up/sdp/models.py`,
   following `rrsModelTrain.m`: z-scored δRrs'' predictors; for each of 100
   permutations an outer random 75/25 split, an inner 5-fold CV selecting the
   number of PCs (≤ max_pcs) by MAE, coefficients back-transformed to A_m(λ),
   C_m, and outputs clipped ≥0. Report mean±SD R² and normalized MAD per
   pigment against Table 2. Reconstruct all samples from the median of the
   100 models (with below-LOD values set to zero) for Fig. 6 and for the
   hierarchical clustering of modeled pigment ratios (Fig. 3B). Compare our
   coefficient spectra with the port's trained coefficients. Also run the
   paper's variants (Rrs'+Rrs''; 5 and 10 nm) at least for Tchla, Fuco and
   HexFuco.
   Use Opus 5.5. Log your work.

4. **Honest diagnostics: skill beyond Tchla and leakage.** Read this file.
   In `epft_up/sdp/validate.py`, add:
   (a) the **Tchla null model**: each pigment predicted from retrieved Tchla
   alone (log–log fit, using GSM Tchla and, separately, OC4/OCI-like
   band-ratio Tchla), scored exactly like Table 2;
   (b) evaluation of PCR and null on **pigment:Tchla ratios**, in linear and
   log space (zeros → ½·LOD via one swappable function, with a sensitivity
   check);
   (c) **leave-one-campaign-out** CV next to the paper's random splits.
   Summarize as a table of skill *above* the null model, per pigment and per
   validation scheme, and state plainly which pigments retain skill beyond
   Tchla.
   Use Opus 5.5. Log your work.

5. **What does the residual buy? Source-space comparison.** Read this file.
   Under identical PCR + LOCO settings, compare predictors: raw Rrs, Rrs'',
   δRrs, δRrs'', El Hourany M1 (spline residual) and M3 (Savitzky–Golay rrs''),
   and δRrs ⊕ {GSM Tchla, a_dg(443), b_bp(443)}, on both absolute pigments and
   ratios. Then test the sensitivity of δRrs to its baseline by swapping the
   NOMAD A,B for a featureless/alternative a_ph shape. Interpret in terms of
   what δRrs removes (Q&A §C-10, §E).
   Use Opus 5.5. Log your work.

6. **Maths/stats section.** Read this file. Write the Maths/Stats section of
   the report, with supporting scripts and figures:
   (i) the linear-operator argument, p̂ = wᵀDδRrs = (Dᵀw)ᵀδRrs, showing
   numerically that the effective weight spectra on δRrs from PCR-on-δRrs''
   can be expressed and compared directly, and what prior the
   derivative + z-score + PC truncation implies;
   (ii) noise propagation through smoothing + finite differences (white and
   correlated), the effective degrees of freedom of δRrs and δRrs'' at
   1/2.5/5/10 nm, using the IOPtics `'pct:X'` (+ floor) and `'pace'` noise
   models;
   (iii) PCR as a shrinkage estimator vs ridge/PLS (spectral-filter view).
   Use Opus 5.5. Log your work.

7. **Alternatives I: learned wavelength weighting.** Read this file.
   Implement and compare under LOCO, with ratios (log-ratio) as the primary
   target and absolute pigments for comparison to the paper:
   ridge and PLS per pigment; the **shared low-rank W** (reduced-rank
   multi-output regression on δRrs ⊕ GSM parameters) with noise whitening
   Σ_n^{-1/2} and a smoothness prior on the weight spectra; rank selected by
   CV. Plot the learned weight spectra and relate the selected rank to the
   4–5 group ceiling (Kramer & Siegel 2019; Kramer+2020). Report which
   wavelengths carry the weight and whether the gains over PCR survive LOCO.
   Use Opus 5.5. Log your work.

8. **Alternatives II: uncertainty.** Read this file. For the best model(s)
   from #7, build a Bayesian version (e.g. Gaussian priors matching the
   smoothness/low-rank structure; closed-form or sampled) that gives
   **per-retrieval predictive intervals**; **propagate input Rrs uncertainty**
   (in-situ-like and PACE OCI `'pace'` noise) through δRrs to the
   predictions; and run **coverage/calibration checks** under LOCO (e.g. the
   fraction of observations inside the 68/95% intervals, PIT histograms).
   Report how skill and coverage degrade from in-situ to OCI-level noise
   and at 5 nm sampling.
   Use Opus 5.5. Log your work.

9. **Independent hold-out: EXPORTS North Atlantic 2021.** Read this file.
   Fetch the 17 EXPORTS-NA (May 2021) HPLC + hyperspectral Rrs matchups used
   by Kramer+2024 from SeaBASS into `$OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA/`,
   processing them as in #1–#2 (stop and ask if the matchups can't be
   reconstructed unambiguously). Without retraining, apply the reproduced PCR,
   the Tchla null model and the best alternatives (with their intervals), and
   report skill and coverage on this campaign.
   Use Opus 5.5. Log your work.

10. **Write the report.** Read this file. Finish
    `reports/SDP_Claude_Report.md` for J. X. Prochaska as primary reader,
    written so any ocean-color scientist can follow it: literature summary
    (incl. Lange 2020 as a brief comparison), data and provenance,
    reproduction (with any deviations from the paper and their causes), maths
    and statistics, diagnostics, alternatives, uncertainty, hold-out results,
    limitations, and recommendations (incl. what PACE application would
    require). Every number must trace to a script in the repo.
    Use Opus 5.5. Log your work.

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
>A. Use your default
2. **Missing links in the chain.** El Hourany & Kramer (2026) lean heavily on
   Kramer et al. (2024), which uses the δRrs residual to define *optical*
   communities, and every paper leans on Kramer & Siegel (2019) for the global
   HPLC grouping. Also, Kramer 2022's key details (Table S6 A(λ),B(λ)
   coefficients; the 1st+2nd-derivative and 5/10 nm variants; Fig. S5
   coefficient spectra) live in its **Supplementary Material**. Can you add
   Kramer+2024, Kramer & Siegel 2019, and the Kramer 2022 (and ideally Catlett
   2018) supplements to `context/papers/`?
>A. I have added the Kramer2019.pdf and Kramer2020.pdf Kramer_etal_2024.pdf to the context/papers/ directory.  Do you mean this DOI for the Kramer 2022 Supp: https://doi.pangaea.de/10.1594/PANGAEA.937536 ?
3. **Report audience & form.** Who is the report for (you / the group / an
   external paper or white paper)? Markdown in `reports/` like
   `MOANA_Claude_Report.md`? *[default: Markdown in `reports/`, written for
   ocean-optics scientists]*
>A. I am the primary audience.  But, it should be understandable by any ocean color scientist.

#### B. Reproducing Kramer+2022

4. **Data source.** The paired dataset is on PANGAEA
   (doi:10.1594/PANGAEA.937536). Is it OK to fetch it into `epft_up/data/`
   (or wherever you keep external data), and do you already have a local copy
   or know whether it holds the 178 raw or the 145 QC'd spectra?
>A. Yes, and put the data in `$OS_COLOR/data/PANAGEA/Kramer2022/`
5. **Reference code.** El Hourany+2026 list three implementations:
   `sashajane19/Rrs_pigments` (original, likely MATLAB),
   `max-danenhower/rrs-SDP-pigments` (a Python port?), and NASA's
   `oci_sdp` notebook (suggesting this is now a PACE/OCI product). Which is
   the reference: (a) write an independent Python version from the paper and
   use the repos only to check; (b) port/wrap one of them; (c) treat the NASA
   SDP notebook as the benchmark? *[default: (a) independent re-derivation,
   checked against the repos, as we did for MOANA]*
>A. Use your default
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
>A. Use your default
7. **Catlett & Siegel 2018 too?** That paper works on a_ph (PnB, 491 samples),
   with 1st+2nd derivatives, a 15 nm Hamming filter, and 500 permutations.
   Reproduce it as well (if the PnB data are accessible), or treat it as
   background only? *[default: background only]*
>A. Use your default

#### C. Maths/stats: what should the discussion probe?

These are the issues I think matter most. Tell me which to emphasize, drop, or add.

8. **Skill beyond Tchla.** Accessory pigments co-vary strongly with Tchla
   (Kramer 2020 normalizes by Tchla for exactly this reason), and Kramer 2022
   models *absolute* concentrations with R² computed in *linear* space. My
   suspicion is that much of the headline R² (0.4–0.7) is Tchla retrieval in
   disguise. Should a central test be a **null model**, i.e. pigment predicted
   from retrieved Tchla alone (e.g. log–log from the GSM fit or OC4), plus
   evaluation on pigment:Tchla ratios? *[default: yes, central]*
>A. Yes, that is a very good idea
9. **Validation design.** The "100-fold cross-validation" is really repeated
   random 75/25 subsampling across 8 cruises whose samples are spatially and
   temporally correlated, so it likely leaks. Should we add
   **leave-one-campaign-out** CV as the honest estimate? (El Hourany use
   temporal blocking; Lange use 80/20 bootstraps.) *[default: yes]*
>A. use your default
10. **What the residual actually removes.** δRrs = Rrs − Rrs_mod, where
    Rrs_mod has a_ph = A(λ)·Tchla^B(λ) with NOMAD-mean *spectral features*. So
    δRrs'' ≈ Rrs'' − (curvature predicted by a mean-community a_ph at that
    Tchla): it is a compositional *anomaly* relative to NOMAD, not
    "phytoplankton signal minus CDOM/NAP" (the second derivative already
    kills most of the smooth adg/bbp shapes). Kramer found Rrs'+Rrs''
    performs comparably. Should I quantify how much the residual step buys
    over plain Rrs'' (and over El Hourany's spline residual M1)?
    *[default: yes]*
>A. Yes. 
11. **Noise and resolution.** A finite-difference 2nd derivative amplifies
    white noise by ~√6/Δλ², and a 5 nm moving mean followed by a 1 nm
    derivative leaves highly correlated channels (the Cael+2020 information
    content argument). Do you want a formal treatment (effective degrees of
    freedom, noise propagation with PACE-OCI-like per-band uncertainties) as
    part of the "maths" section? *[default: yes, including PACE OCI noise]*
>A. Yes, the default
12. **n ≪ p regularization.** With 145 samples and ~290 correlated predictors,
    PCR is one of several shrinkage estimators. For the discussion, compare
    PCR with ridge/PLS (and possibly a Bayesian linear model that gives
    predictive uncertainties natively)? *[default: yes]*
>A. Yes, the default

#### D. Alternatives and project goals (uncertainty + provenance)

13. **Target variable.** Which should the alternatives aim at: (a) absolute
    pigment concentrations (Kramer 2022), (b) pigment:Tchla ratios, possibly
    treated as compositional data (log-ratios), (c) EOF amplitudes (Catlett),
    or (d) discrete assemblage classes (El Hourany SOM + Random Forest)?
    *[default: (b) as primary with (a) for comparison to the paper]*
>A. Let's start with your default
14. **Alternatives you already have in mind?** Candidates I see: inverting
    a_ph(λ) first (e.g. via an IOP inversion from IOPtics) and doing the
    derivative analysis on a_ph (closer to Catlett); Gaussian-band
    decomposition of a_ph; PLS/ridge/GP; El Hourany's M1/M2/M3 source spaces;
    the MOANA approach. Should any of these be prioritized or excluded? Is the
    IOPtics machinery meant to be used here?
>A. My intuition is that the majority of signal is in δRrs = Rrs − Rrs_mod.  That is, what else can there be?  But, and this is key, we are going to want to learn a weighting matrix that gives preference as a function of wavelength.  
15. **Uncertainty deliverable.** For EPFT-UP, what counts as success on
    uncertainty: per-retrieval predictive intervals, propagation of the input
    Rrs uncertainty, and/or calibration checks (coverage) under
    leave-one-campaign-out? *[default: all three, at least for the
    best alternative]*
>A. Use your default
16. **PACE application.** Is applying the reproduced and alternative models to
    PACE OCI L2/L3 Rrs (and comparing with NASA's SDP product) in scope for
    this prompt series, or strictly in situ? *[default: in situ first, PACE
    later as its own prompt]*
>A. Use your default

#### E. Follow-up from Claude (2026-09-29)

I reviewed your answers, read the two new papers (Kramer & Siegel 2019;
Kramer et al. 2024), pulled the PANGAEA *metadata* (no data yet), and skimmed
the two public code repos. First, what I learned that changes the plan, then
the new questions.

**What changed**

- **Re Q2 (the PANGAEA DOI):** 10.1594/PANGAEA.937536 is the *data*, not the
  supplement. The supplement is Elsevier's "Appendix A" at
  doi:10.1016/j.rse.2021.112879 (usually `mmc1.docx`). It is now only
  nice-to-have, because the repos include the key inputs:
  `sashajane19/Rrs_pigments` (MATLAB) ships `aph_A_B_Coeffs_Sasha_RSE_paper.txt`
  (the Table S6 A(λ),B(λ)), `aw_mcf16_350_700_1nm.txt` (Mason+2016 water
  absorption), `betasw_ZHH2009.m`, and `gsm_invert.m`.
  `max-danenhower/rrs-SDP-pigments` is a Python port that also ships trained
  coefficients (`original_a_coefs.xlsx`, `original_c_coefs.xlsx`), which gives
  us a direct numerical check of our PCR.
- **The PCR recipe is settled by `rrsModelTrain.m`** (Catlett's code, edited by
  Kramer). For each of n_permutations: a random 75/25 outer split; inside the
  75%, a **k=5-fold CV picks the number of PCs** (≤ max_pcs) by the chosen
  metric (MAE for Kramer 2022); pigment outputs are clipped at ≥0; `rng(1)`.
  That resolves the Q6 ambiguities. We can't reproduce MATLAB's random splits
  exactly, only statistically.
- **PANGAEA holds only the 145 QC'd samples.** The metadata show 336 columns
  and ~48k data points, about 145 rows × (301 Rrs bands + ~35 metadata/pigment
  columns), and Kramer+2024 confirms 145. So the visual QC (178 → 145)
  **cannot** be reproduced from this dataset. I will instead check that an
  objective rule (the Lange-style |Rrs''| cut at 610–660 nm) keeps all 145, and
  record the QC step as "inherited, unverifiable". The file also carries
  in-situ T and S and extra pigments (Tchlb, Tchlc, αβ-carotene, Diadino, …).
- **Kramer+2024** adds 17 EXPORTS North Atlantic (May 2021) samples (N=162),
  processed the same way, and finds the δRrs-based communities agree with
  HPLC communities 74% of the time (only 3 communities; identical at 5 nm).
- **Kramer & Siegel 2019**: 4 groups globally (diatoms+dinos merged), up to 6
  locally. That caps what any global in-situ training set can resolve.

**On your Q14 comment (δRrs plus a learned wavelength weighting).** I agree,
and the maths makes it the natural framing:

- For *any* linear estimator, p̂ = wᵀ D δRrs = (Dᵀw)ᵀ δRrs, where D is the
  (fixed) smoothing + 2nd-difference operator. **The derivative adds no
  information for a linear model**: it only changes the implicit prior on
  the effective weights, through the PCA basis, the per-wavelength z-scoring
  and the PC truncation, and it changes how noise is weighted. So the object
  to learn is the weight spectrum on δRrs itself. "2nd derivative + PCR" is
  one hand-picked prior on it, and ridge/PLS/GP/Bayesian priors are others we
  can compare directly.
- "What else can there be?" Two things the residual throws away or mixes in:
  (1) the GSM fit parameters (Tchla, a_dg(443), b_bp(443)) themselves carry
  most of the Tchla-covariant skill. That matters for absolute
  concentrations, less for ratios. (2) δRrs also contains *structured
  non-phytoplankton* signal: GSM mis-specification (b_bp spectral shape,
  Raman, Chl fluorescence near 683 nm, a_w T/S dependence) and instrument
  artifacts. I propose predictors = δRrs ⊕ {GSM parameters}, with the
  weighting learned on δRrs.
- Note also that δRrs depends on the NOMAD A(λ),B(λ) baseline. Changing
  that baseline changes the "anomaly". It could itself be learned or
  swapped, e.g. for an a_ph model from IOPtics.

**New questions**

17. **Data path.** `$OS_COLOR` = `/Users/xavier/Projects/Oceanography/data/Color/`,
    so `$OS_COLOR/data/PANAGEA/Kramer2022/` expands to
    `.../data/Color/data/PANAGEA/Kramer2022/`, with "data" twice and
    "PANAGEA" misspelled. Use that exact path, or `$OS_COLOR/PANGAEA/Kramer2022/`?
    *[default: `$OS_COLOR/PANGAEA/Kramer2022/`]*
>A. Use your default
18. **What does the "weighting matrix" mean to you?** Options:
    (a) one weight spectrum w_m(λ) per pigment (standard regression);
    (b) a **shared, low-rank W** (pigments × λ) fitted jointly to all pigments
    (reduced-rank / multi-output regression), which exploits pigment
    covariance and is where the 4–5 group ceiling shows up naturally;
    (c) a metric M(λ,λ′) used to compare spectra (for classification /
    community detection à la Kramer+2024 / El Hourany);
    (d) wavelength weights *inside the GSM fit* (i.e. which bands define
    Rrs_mod and hence δRrs). *[default: (b), with a noise-whitening
    Σ_n^{-1/2} and a smoothness prior, compared against (a)]*
>A. Use your default
19. **Independent test set.** Should we pull the 17 EXPORTS-NA 2021 samples
    (SeaBASS) as a genuinely held-out campaign for the reproduced and
    alternative models? It's the cleanest out-of-sample check available.
    *[default: yes, after the reproduction]*
>A. Use your default
20. **Zeros in pigment ratios.** Kramer sets below-LOD values to zero, which
    breaks log-ratios for Q13(b). *[default: replace zeros with ½·LOD
    (multiplicative replacement) and report sensitivity to that choice]*
>A. Use your default.  Tom Jordan is fitting a model to a log-normal distribution for the concentrations.  We will do this eventually as well, but I don't want to steal his thunder.
21. **Rrs uncertainty model.** PANGAEA has no per-sample Rrs uncertainties.
    IOPtics already has a sweep-level noise-model framework (and a GSM). For
    Q11/Q15, should I reuse the IOPtics noise models (which ones?) and add a
    PACE-OCI per-band model, or define a simple HyperPro-like model here
    (relative % plus an absolute floor)? *[default: reuse IOPtics'
    framework; add an OCI model if one isn't there]*
>A. Use your default
22. **b_bw temperature/salinity.** The paper used WOA climatology T/S, while
    PANGAEA includes in-situ T/S. *[default: in-situ T/S; verify the effect
    on δRrs is negligible (it should be ≪ the residual)]*
>A. Use your default

#### F. Follow-up from Claude (2026-09-29, round 3)

All of Q17–Q22 take the defaults, so the design is essentially settled.
Remarks:

- **Q20 / Tom Jordan's log-normal model.** Understood; we won't build a
  distributional model of the concentrations. For the record, the principled
  version of what Tom is doing is a *censored* (Tobit-style) log-normal
  likelihood, where below-LOD values contribute P(c < LOD) rather than a
  number. I'll keep the zero handling as a single swappable function
  (½·LOD by default), so his model can drop in later without touching
  anything else.
- **Q21 is concrete now.** IOPtics `ioptics/noise.py` already provides
  `'pace'` (per-band OCI σ from `ocpy.satellites.pace.gen_noise_vector` on the
  native grid), `'pct:X'`, and `'insitu'`. So there is no OCI model to add. For
  the in-situ Kramer spectra I'll use `'pct:X'` with an absolute floor as the
  HyperPro-like model, and `'pace'` for the "what survives at OCI noise"
  experiments. That means EPFT-UP will import from IOPtics/ocpy; flag it if
  you'd rather keep EPFT-UP independent of IOPtics.

**Remaining questions (last round, I think)**

23. **Licensing of the reference code.** Neither `sashajane19/Rrs_pigments`
    nor `max-danenhower/rrs-SDP-pigments` has a license. So we shouldn't copy
    their files (e.g. the A(λ),B(λ) table, trained coefficients) into this
    repo. *[default: fetch them at run time into
    `$OS_COLOR/PANGAEA/Kramer2022/ref/`, record the commit SHA in provenance,
    commit nothing from them, and cite them in the report]* If you'd like to
    redistribute the A,B table, it's worth asking Sasha Kramer to add a
    license.
>A. Yes, that is fine
24. **Code layout.** *[default: a subpackage `epft_up/sdp/` with `data.py`
    (PANGAEA ingest + provenance), `gsm.py` (the GSM-like δRrs model),
    `spectral.py` (smoothing, finite differences, alternative transforms incl.
    El Hourany M1), `models.py` (PCR, ridge/PLS, reduced-rank W, Bayesian),
    `validate.py` (random-split, LOCO, Tchla null model); tests in
    `epft_up/tests/`; run scripts in `scripts/sdp/`; report at
    `reports/SDP_Claude_Report.md`]*
>A. Use your default
25. **Proposed order of work** (for you to turn into prompts; I won't start
    until you do):
    1. Data: fetch PANGAEA 937536 + reference repos, ingest, provenance, and
       sanity plots (Kramer Fig. 1–2).
    2. Reproduce the δRrs: GSM fit (Tchla vs. HPLC, Fig. 4B) and δRrs
       (Fig. 2C), cross-checked against the Python port.
    3. Reproduce the PCR: Table 2, Figs. 3 & 6, and coefficient spectra; check
       against the port's trained coefficients.
    4. Honest diagnostics: Tchla null model, pigment:Tchla ratios, LOCO CV,
       residual vs. plain Rrs'' vs. M1.
    5. Maths section: the linear-operator argument, noise propagation /
       effective DoF, OCI noise, 1/5/10 nm.
    6. Alternatives: ridge/PLS, reduced-rank W with noise whitening +
       smoothness prior, Bayesian with predictive intervals + coverage;
       predictors δRrs ⊕ GSM parameters.
    7. EXPORTS-NA 2021 holdout; then write the report.
    Is that the order you want, or should anything move (e.g. start the maths
    section alongside step 2)?
>A. That looks good.  Please generate a series of prompts for this in the Prompts/Execution section above.

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

### 2026-09-29 (Setup #2: reviewed answers; follow-up comments and questions)

Reviewed the answers to Q1–Q16 (mostly defaults). Accepted: Lange as a brief
comparison; report in `reports/` for J.X.P. but readable by any ocean-color
scientist; independent Python re-derivation checked against the repos; target
Table 2 + Fig. 3/6 with objective QC; Catlett as background; Tchla null model
central; leave-one-campaign-out CV; quantify residual vs plain Rrs'' vs M1;
formal noise/DoF treatment incl. PACE OCI; PCR vs ridge/PLS/Bayesian; ratios
primary; all three uncertainty deliverables; in situ first.

New reading: **Kramer & Siegel 2019** (JGR; 4,480 global surface HPLC samples,
17 pigments normalized to Tchla; clustering + EOFs → 4 global groups with
diatoms+dinos merged, 4–6 at six time-series sites) and **Kramer et al. 2024**
(Opt. Express 32, 34482; the SDP δRrs from Kramer 2022 on 145 + 17 EXPORTS-NA
2021 samples, N=162; WGCNA (β=12) + modularity community detection on δRrs
and on 12 pigment ratios; 3 communities each, 74% agreement, F1 0.62–0.79,
Mantel ρ=0.49; unchanged at 5 nm).

Learned about resources (read-only checks; nothing downloaded to the repo):
- PANGAEA 937536 metadata (fetched to scratchpad only) show 336 parameters and
  47,995 data points, i.e. the 145 QC'd samples with 1 nm Rrs 400–700,
  in-situ T/S and an extended pigment list. The 33 visually rejected spectra
  are not included, so the QC step can't be reproduced from this source.
- `sashajane19/Rrs_pigments` (MATLAB): `Kramer_hyperRrs.m`, `gsm_invert.m`,
  `gsm_cost.m`, `rrsModelTrain.m`, `aph_A_B_Coeffs_Sasha_RSE_paper.txt`
  (Table S6), `aw_mcf16_350_700_1nm.txt`, `betasw_ZHH2009.m`, test data.
  `rrsModelTrain.m` settles the PCR procedure: outer random 75/25 splits ×
  n_permutations; an inner 5-fold CV chooses the number of PCs ≤ max_pcs by
  the selected metric; clip ≥0 for pigments; rng(1).
- `max-danenhower/rrs-SDP-pigments`: Python port with trained A/C coefficient
  spreadsheets, useful for a numerical cross-check.
- The Kramer 2022 supplement is Elsevier Appendix A (not PANGAEA). It is no
  longer critical.
- IOPtics has a sweep-level noise-model framework and a registered GSM
  algorithm; both are candidates for reuse.

Key comment recorded in Q&A §E: for any linear estimator the 2nd derivative
is a fixed linear operator folded into the weights, (Dᵀw)ᵀδRrs. So it adds no
information and acts only as an implicit prior and noise weighting. This
supports the user's "learn a wavelength weighting on δRrs" framing. Posed
Q17–Q22 (data path, meaning of the weighting matrix, EXPORTS-NA holdout, zeros
in log-ratios, Rrs noise model, T/S for b_bw).

### 2026-09-29 (Setup #3: reviewed round-2 answers; final clarifications)

All of Q17–Q22 take the defaults: data at `$OS_COLOR/PANGAEA/Kramer2022/`;
weighting = shared low-rank W with noise whitening + smoothness prior, compared
to per-pigment w; EXPORTS-NA 2021 as a held-out campaign after the
reproduction; zeros → ½·LOD with sensitivity analysis (a distributional
log-normal model is deferred to Tom Jordan's work); IOPtics noise framework;
in-situ T/S for b_bw.

Checked the repos (read-only):
- IOPtics `ioptics/noise.py` already offers `'pace'` (OCI per-band σ via
  `ocpy.satellites.pace.gen_noise_vector`), `'pct:X'` and `'insitu'`, so no new
  OCI noise model is needed. Using it makes EPFT-UP depend on IOPtics/ocpy.
- Neither `sashajane19/Rrs_pigments` nor `max-danenhower/rrs-SDP-pigments`
  has a license (GitHub API: none). Their files shouldn't be vendored; the
  default is to fetch them at run time and pin the commit SHA.

Posed Q23–Q25 in Q&A §F: licensing/vendoring of the reference files, the
proposed `epft_up/sdp/` layout, and the proposed seven-step order of work.
No calculations or downloads yet.

### 2026-10-02 (Setup #4: round-3 answers; wrote the Execution prompts)

Round-3 answers: Q23 accepted (fetch the unlicensed reference repos at run
time into `$OS_COLOR/PANGAEA/Kramer2022/ref/`, pin SHAs, vendor nothing); Q24
default code layout (`epft_up/sdp/` with data/gsm/spectral/models/validate,
tests in `epft_up/tests/`, scripts in `scripts/sdp/`, report at
`reports/SDP_Claude_Report.md`); Q25 order of work approved, with a request
to turn it into prompts.

Wrote ten prompts under Prompts → Execution, plus a shared conventions
paragraph (code and data locations, ocean14, provenance, report section and
log updates, and a stop-and-ask rule for unexplained disagreements with the
paper). The seven proposed steps were expanded to ten:
1. data ingest
2. δRrs/GSM reproduction
3. PCR reproduction
4. Tchla null model + ratios + leave-one-campaign-out
5. source-space comparison (what the residual buys)
6. maths/stats section
7. learned weighting (ridge/PLS/low-rank W)
8. Bayesian uncertainty, propagation and coverage
9. EXPORTS-NA 2021 hold-out
10. final report

Step 4 of the original plan was split into #4 and #5, step 6 into #7 and
#8, and step 7 into #9 and #10. No calculations or downloads yet.

### 2026-10-02 (Execution #1: data ingest and provenance)

Fetched and pinned the inputs under `$OS_COLOR/PANGAEA/Kramer2022/`:
`PANGAEA_937536.tab` (PANGAEA textfile export, SHA-256 `63da2981…5f910f3b`)
and, in `ref/`, `sashajane19/Rrs_pigments` @ `b3e3662` and
`max-danenhower/rrs-SDP-pigments` @ `fb17c2f`. Nothing is vendored. Note that
`$OS_COLOR` already contains a differently spelled `PANAGEA/V3/` (another
dataset). I followed the agreed `PANGAEA/` spelling; the location is one
function, `epft_up.sdp.data.kramer2022_dir()`, if you want it moved.

New code:
- `epft_up/sdp/__init__.py`
- `epft_up/sdp/data.py`: pins (`PANGAEA_SHA256`, `REF_REPOS`); pigment and
  campaign maps; transcribed `TABLE1`; `load_kramer2022()` returns a
  `KramerData` (wave, Rrs, meta, pigments, provenance); `table1_comparison()`;
  `second_derivative_qc()` (Lange-style screen)
- `scripts/sdp/fetch_kramer2022.py`: fetch and verify, `--verify`
- `scripts/sdp/ingest_kramer2022.py`: Table 1 check, QC, the pre-smoothing
  test, Figs. 1/2A, `reports/figures/sdp/ingest_summary.json`
- `epft_up/tests/test_sdp_data.py`: 6 synthetic tests plus 1 real-data test,
  with a new `needs_kramer2022` marker in `conftest.py`. Full suite: 45 passed.
- `reports/SDP_Claude_Report.md`: outline, with §3 Data written

Findings:
- N=145; per-campaign counts match Table 1 exactly. The three Polarstern legs
  (7+1+18) are "ANT". Tchla statistics match to the paper's precision (overall
  0.019–4.151, median 0.110).
- **The deposited Rrs are already 5 nm moving-mean smoothed and then rounded
  to 1e-6 sr⁻¹.** At the boxcar zeros (0.2, 0.4 cycles/nm) the power of the
  2nd difference sits at the rounding floor (ratio 0.7–1.3 in every
  campaign); at the sidelobes it is 3–4× above. So #2 must not re-smooth.
  Kramer's MATLAB code doesn't smooth either. The rounding is ~1% of red Rrs
  and is a quantifiable noise term for Rrs'' (σ≈7e-7 sr⁻¹ nm⁻² per band),
  which feeds into #6. My first version of this test was mis-normalized
  (the 2nd difference of white noise isn't white). It was replaced by the
  explicit comparison with the quantization floor.
- The Lange |Rrs''| > 2e-4 screen rejects 0 of 145 (largest 1.1e-5), so it is
  uninformative on this deposit and cannot stand in for the visual QC (the
  33 rejected spectra aren't deposited).
- Many pigments are mostly zeros: Neo 0.74, Allo 0.70, Perid 0.54, DVchla
  0.35, MVchlb 0.33. The ½·LOD log-ratio sensitivity (Q20) will matter a lot.
  Kramer's ">75% below detection" exclusion rule doesn't reproduce from the
  deposit (DVchlb and Lut, 0.68 each, were dropped while Neo, 0.74, was kept).
- The `Rrs_pigments` README grants free use of its code and data. This is
  relevant to Q23, but nothing has been vendored.

Not done or open: `cartopy` (used only by the ingest script for Fig. 1) is not
in `requirements.txt`.

### 2026-10-02 (Execution #1 follow-up: directory spelling, cartopy)

Following the user's answers to the #1 open items:
- **Spelling.** Kept `$OS_COLOR/PANGAEA/` and moved the legacy
  `$OS_COLOR/PANAGEA/V3/` (the Valente et al. 2022 OC-CCI v3 compilation,
  PANGAEA.941318) to `$OS_COLOR/PANGAEA/V3/`, then removed the empty
  `PANAGEA/`. This is safe because `ocpy.insitu.pangaea.pangaea_path()`
  already searches `PANGAEA/V3` first and `PANAGEA/V3` only as a legacy
  fallback; it now resolves to the new location, which I verified. No other
  Python code in `~/Oceanography/python` references the old spelling. The
  IOPtics design doc and PAB `docs/context.md` still mention it in prose
  (not changed).
- **Licensing.** The user accepted that the `Rrs_pigments` README grant is
  fine. Still nothing is vendored.
- **cartopy** added to `requirements.txt` and to `install_requires` in
  `setup.py` (kept in sync, as `requirements.txt` asks).

### 2026-10-02 (Execution #2: reproduce δRrs, the GSM-like residual)

New code:
- `epft_up/sdp/gsm.py`: an independent implementation of Kramer 2022
  Eqs. 1–6, comprising Lee rrs conversions; Zhang et al. (2009) b_sw written
  from its equations (ocpy's `betasw_ZHH2009` can't be used: it raises "not
  successfully converted" and has a Boltzmann-constant typo, 1.38e-22);
  Carder S_dg; Lee η; the Gordon quadratic; `load_ref_tables()` (A,B,a_w from
  the pinned checkout, with checksums); `fit_gsm()` with `method='kramer'`
  (Nelder–Mead with MATLAB fminsearch settings) or `'lsq'` (bounded
  log-parameters); a `GSMFit` result with δRrs.
- `epft_up/sdp/spectral.py`: `moving_mean`, `trim_edges`, `kramer_preprocess`
  (the paper's §2.2 recipe for new, unsmoothed spectra).
- `scripts/sdp/fetch_woa_ts.py`: WOA23 ¼° monthly surface T/S at each sample
  via OPeNDAP → `$OS_COLOR/PANGAEA/Kramer2022/woa23_surface_ts.csv` (+ JSON).
- `scripts/sdp/reproduce_gsm.py`: fit, Fig. 4, Fig. 2B–C, port cross-check,
  sensitivities; product `$OS_COLOR/.../products/gsm_dRrs_insitu.npz`.
- `epft_up/tests/test_sdp_gsm.py`: 9 tier-1 tests (synthetic fixed-point
  spectra recovered by both optimizers, Zhang physics checks, smoothing) and
  2 tier-2 tests (real tables; Fig. 4B statistics). Full suite: 56 passed.
- Report §4.1 written.

Findings:
- **Fig. 4 reproduced.** OC4v6: y = 0.873x − 0.138, R² 0.746 (paper 0.87x −
  0.14, 0.75). GSM: R² 0.717 and slope 0.942 on all 145. Dropping one
  degenerate fit (SABOR 2014-07-31, idx 69: Tchla → 4e-4, a_dg(443) = 0.12)
  gives y = 0.961x − 0.091, R² 0.864, the paper's 0.96x − 0.093, 0.86. That
  point is absent from the paper's Fig. 4B. The likely cause is MATLAB's
  complex arithmetic when Nelder–Mead trials Tchla ≤ 0 (22 of 366 trials do,
  for this sample). It cannot be verified without MATLAB. Not a
  stop-and-ask: the cause is identified to a single sample.
- **The printed Eq. 5 has a sign typo.** Both codes use Carder (1999);
  taken literally, Eq. 5 gives R² 0.07. The paper's text says η uses
  rrs490/555, but the code uses rrs440/555 (Lee 2002). The 490 variant gives
  R² 0.62.
- Python port vs ours (same η convention): max |ΔδRrs| = 4.5e-8 sr⁻¹, below
  the deposit rounding; b_sw identical. The port uses above-water Rrs for η
  (0.08% effect).
- WOA23 vs in-situ T/S: 0.12% relative RMS on δRrs (max 1.2e-6). Q22 is
  closed: in-situ is fine.
- Five ANT spectra are fitted with b_bp(443) < 0 (unbounded fit, faithful to
  MATLAB). A bounded fit changes only those five.
- An accidental second 5 nm smoothing changes δRrs by 4% RMS, 10–100× the
  other choices. This confirms #1's warning not to re-smooth.
- `Kramer_rrs_testdata.mat` in `Rrs_pigments` holds 17 EXPORTS-NA spectra
  with T, S and chl (no pigments), which is likely the #9 hold-out. Its
  pigments must still come from SeaBASS.

Open for #3: run the PCR both with and without the degenerate SABOR
spectrum.

### 2026-10-02 (Execution #3: reproduce the PCR model, Table 2, Figs. 3 & 6)

New code:
- `epft_up/sdp/models.py`:
  - `train_pcr_kramer()` is a re-implementation of `rrsModelTrain.m`
    (100 × 75/25 splits; inner 5-fold CV picks ≤ 30 PCs by MAE; outputs
    ≥ 0; the original's quirks are kept, namely own-statistics scaling of
    the inner validation folds and fold-averaging in standardized units).
  - An exact one-fit shortcut computes all nested PC models.
  - `predict_ensemble()` returns the median of the ensemble with clipping
    and an LOD.
- `epft_up/sdp/spectral.py`: `difference()` (MATLAB `diff` semantics,
  orders 1 and 2, subsampling `step`) and `derivative_features()`.
- `scripts/sdp/reproduce_pcr.py`: 13 pigments × {145, 144}; Table 2; Figs.
  6 and 3; a cross-check with the port's trained coefficients; supplementary
  variants. Products go to `$OS_COLOR/.../products/pcr_rrsD2_1nm.npz`.
- `epft_up/tests/test_sdp_models.py`: 7 tier-1 tests (including the
  nested-PC shortcut vs explicit OLS) and 1 tier-2 test (the Tchla row of
  Table 2). Full suite: 63 passed.
- Report §4.2 written.

Settings confirmed from `Kramer_Rrs_pigments.m`: max_pcs = 30, MAE
selection, k = 5, 100 permutations, predictors `diff(Rrs residual, 2, 2)`
(plain second difference, 299 bands at 401–699 nm). The 1st derivative is a
*forward* `diff`, not Catlett's centred Eq. 1.

Findings:
- **Table 2 reproduced:** all 13 pigments' mean R² lie within one quoted SD
  (mean Δ ≈ −0.2 SD; worst ButFuco −0.7 SD and DVchla −0.5 SD). Table 2's
  "normalized MAD" is MAE / mean *modelled* value (Tchla 0.507 vs 0.498).
  Results are the same with and without the degenerate SABOR spectrum
  (|ΔR²| ≤ 0.03).
- **Kramer's original trained coefficients** (port `original_*_coefs.xlsx`)
  applied to our δRrs'' reproduce our full-reconstruction R² to ±0.05 for
  every pigment. The median A(λ) correlate at r = 0.85–0.90 (Zea and
  DVchla ≈ 0.5, poorly constrained). Strong confirmation that our δRrs and
  PCR are Kramer's.
- **Fig. 6 is in log₁₀ space**, not linear (my first version was linear,
  which was wrong). R² matches; slopes are 0.02–0.23 lower than the paper's
  (Tchla 0.88 vs 0.94; Perid 0.55 vs 0.78). Kramer's own coefficients give
  the same low slopes, and no regression convention matches all panels.
  This is an **open discrepancy**, with no uncertainty quoted in the paper
  and Table 2 unaffected, so I did not stop. 17 low-Tchla samples have
  modelled Tchla ≤ 0.
- **LOD proxy:** the smallest non-zero deposited value for pigments with
  zeros; 0.001 (reporting resolution) for Tchla, Zea and Chlc12. The first
  version used Tchla's sample minimum (0.019), which wrongly zeroed
  modelled values. That has been corrected.
- **Fig. 3:** the measured-ratio dendrogram gives exactly the paper's five
  groups. The modelled-ratio dendrogram is fragile and depends on zero
  handling (4 groups now, 5 with the coarser LOD).
- **Variants:** measured Rrs' + Rrs'' ≥ δRrs'' for 11 of 13 pigments. 5 nm
  sampling is flat or better. 10 nm is **not** worse except for Perid,
  contradicting the paper's "notably worse at 10 nm". The paper's
  degradation recipe is unknown (I used subsampling of the smoothed δRrs).
- The z-scoring puts the largest PCR weights at 520–700 nm, where δRrs'' is
  smallest. This is for #6.
- The port's `run_sdp` re-smooths new spectra with a **trailing**
  `rolling(5, min_periods=1)` window (a 2 nm shift). Relevant for PACE
  application.
