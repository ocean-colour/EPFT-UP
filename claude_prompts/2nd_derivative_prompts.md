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

11. **Slides.** Generate a set of slides for the report.  
    Provide context on the what, why, how, and what we found.
    Use Figures as often as you can.  Include ones from the 
    Kramer papers too.
    Avoid using fonts with a size less than 20pt.
    Place the file in the `reports/slides/` directory and call it 
    `SDP_Claude_Slides.pptx`.
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

### 2026-10-02 (Execution #4: honest diagnostics: skill beyond Tchla, ratios, LOCO)

New code:
- `epft_up/sdp/validate.py`:
  - `lod_proxy` (moved here as the single definition) and `replace_zeros`
    (the one swappable ½·LOD rule);
  - `oc4v6`;
  - `random_splits`, `loco_splits`;
  - fitters (`PCRFitter`, `LogLinearFitter` (the Tchla null),
    `ConstantFitter`) and `cross_validate`;
  - `score_linear` (Table 2 style), `score_log`, `per_split_scores`,
    `pooled_predictions`, `paired`, `bootstrap_paired`.
- `epft_up/sdp/models.py`: factored out `fit_pcr_one()` (one model from one
  training set). `train_pcr_kramer` is unchanged numerically (Tchla R²
  0.7342 before and after).
- `scripts/sdp/diagnostics_null_loco.py`:
  - 13 pigments × {PCR, null_GSM, null_OC4, oracle};
  - 12 ratios × {PCR_log, PCR_lin, PCR_derived, nulls, const};
  - random (100) and LOCO (8) schemes, linear and log scoring;
  - a zero-replacement sensitivity (0.1, 0.5, 1.0 × LOD);
  - the N=144 repeat;
  - figures and `diag_summary.json` / `diag_skill_table.csv`.
- `scripts/sdp/reproduce_pcr.py` now takes `lod_proxy` from `validate`
  (Table 2 output byte-identical).
- `epft_up/tests/test_sdp_validate.py`: 7 tests. Full suite: 70 passed.
- Report §6.1–6.2 written.

Findings (N=145; N=144 gives the same verdicts):
- **Absolute concentrations: no skill beyond Tchla.** A two-parameter
  log–log fit on the GSM-retrieved Tchla beats PCR for 11 of 13 pigments
  under random splits (Tchla 0.81 vs 0.72, Fuco 0.72 vs 0.61, Chlc12 0.78 vs
  0.64). Even the OC4 Tchla matches or beats PCR for 6 of 13. PCR wins only
  for DVchla (0.48 vs 0.06) and Zea (0.37 vs 0.20), the pigments that don't
  scale with Tchla.
- **LOCO halves PCR's skill** (Tchla 0.72 → 0.44, HexFuco 0.48 → 0.12,
  DVchla 0.48 → 0.04). Under LOCO, PCR is significantly worse than null_GSM
  for 7 pigments and better for none.
- **Ratios:** ratios derived from SDP's absolute products have R² ≈ 0.02–0.09
  (0.2–0.3 for Zea and DVchla), so they carry no compositional information.
  PCR trained on log ratios has real within-campaign skill (Fuco:Tchla 0.61,
  DVchla 0.69, Zea 0.57, Chlc3 0.55; beating the nulls in ≥ 90% of splits).
  Under LOCO nothing survives significantly. Neo:Tchla appears to, but it is
  an artifact of zero replacement (LOCO R² 0.03/0.20/0.47 for 0.1/0.5/1.0 ×
  LOD). The best candidate is Fuco:Tchla (0.40 vs 0.31; ΔR² CI
  [−0.01, 0.17]). For HexFuco, ButFuco, MVchlb and Viola, LOCO log-ratio RMS
  is no better than a constant ratio.
- Mechanism: held-out campaigns are offset as a block (BIOSOPE Tchla
  predicted 10× too high; RemSensPOC and Tara Med often ≤ 0). The PCR uses
  campaign-specific fine spectral structure (instrument and processing).

Not a stop-and-ask: these are new tests that the paper did not run, not
contradictions of a quoted number. They are, however, the most consequential
result so far, so they are worth discussing before #7.

### 2026-10-03 (Execution #5: source-space comparison; standard Benchmark)

User instruction (2026-10-02): from here on, models must be judged on
**both** held-out-campaign skill and skill above the Tchla null. This is
implemented as a reusable scorecard:
- `validate.Benchmark`: 25 targets (13 absolute pigments plus 12 log10
  ratios), 100 random and LOCO splits, cached baselines (null_GSM,
  null_OC4, const). `evaluate(fitter_factory, label)` returns, per target,
  R², ΔR² vs the best null, the paired-win fraction (random), the bootstrap
  CI (LOCO), `beats_null_random` / `beats_null_loco`, and ratio RMS vs the
  constant. `pcr_factory(X)` gives Kramer's PCR. It reproduces #4's numbers
  exactly.
- `spectral.py`: `spline_residual` (El Hourany M1; csaps p = 0.0005 with the
  penalty rescaled 2.5 nm → 1 nm) and `savgol_second_derivative` (M3; the
  window is not given in the paper, so 7, 11 and 21 nm are tried).
- `gsm.py`: `smoothed_tables(tables, fwhm)` (featureless a_ph baseline).
- `scripts/sdp/source_spaces.py`: 16 predictor spaces through the
  Benchmark, plus a curvature budget; `--replot`. Figures and tables are in
  `reports/figures/sdp/source_spaces_*`.
- Tests: 2 Benchmark tests and 3 transform tests. Full suite: 75 passed. I
  relaxed one of my new tests: M1's heavy spline leaves ~0.4% of a smooth
  exponential, which is expected behaviour.
- Report §6.3 written.

Findings:
- **The 2nd derivative is the harmful step.** All finite-difference spaces
  have the worst LOCO skill (absolute ≈ 0.19–0.21, ratio ≈ 0.16–0.18). The
  undifferentiated δRrs (0.38) and M1 (0.40) do about twice as well. M3's
  LOCO skill rises monotonically with the SG window (7 → 21 nm: abs 0.25 →
  0.37, ratio 0.19 → 0.27).
- The residual helps only without the derivative (δRrs 0.38 vs Rrs 0.32;
  δRrs'' 0.19 ≈ Rrs'' 0.20). The purely empirical M1 matches or beats δRrs,
  so it is envelope removal, not the GSM physics, that helps.
- **The NOMAD a_ph features are irrelevant:** 30/80 nm-smoothed A,B change
  δRrs'' by < 0.5% and leave all scores unchanged.
- **A smoothing-mismatch artifact:** the deposited Rrs are 5 nm smoothed but
  the 1 nm a_w table is not. Its fine structure enters Rrs_mod'', so the
  model curvature has 1.30× the variance of Rrs'' and var(δRrs'') =
  1.52 × var(Rrs''). With a_w smoothed consistently: 0.64× and 0.99×, and
  r = 0.61 with the published δRrs''. This is present in Kramer's pipeline
  too (we match the port's coefficients). PCR is only mildly affected (LOCO
  0.21 vs 0.19).
- Out of campaign, the compositional signals that beat the Tchla null are
  Fuco:Tchla (M1 0.56, M3 21 nm 0.55, δRrs 0.50 vs null 0.31; not with
  δRrs''), DVchla:Tchla (M1 0.49, M3 21 nm 0.45 vs 0.32), Zea:Tchla (M3
  21 nm 0.48 vs 0.35) and Chlc12:Tchla (M1, δRrs_flat80 0.40 vs 0.27).
  Neo:Tchla is a zero artifact. No space beats the null on absolute Tchla.
- Best spaces for #7: δRrs (consistent smoothing), M1, and wide-window M3,
  with log ratios as the target.

### 2026-10-03 (Execution #6: Maths/Statistics section)

New code:
- `epft_up/sdp/theory.py`: `difference_matrix`, `boxcar_matrix`,
  `subsample_matrix`, `effective_weights` (w = Lᵀ A),
  `implied_prior_covariance` (Lᵀ S⁻² L), `propagate_covariance`,
  white/correlated covariances, `dof_signal` (Rodgers d_s via generalized
  eigenvalues, pseudo-whitening), `filter_factors`, `ridge_coefficients`,
  `pcr_coefficients`, `frequency_power`.
- `epft_up/sdp/noise.py`: `pace_sigma` (ocpy `gen_noise_vector`, the
  IOPtics 'pace' model; ioptics itself is not importable in ocean14),
  `pct_floor_sigma` (the IOPtics two-part floor, reproduced), and
  `rounding_sigma`.
- `scripts/sdp/maths_section.py`: parts (i)–(iii), five figures,
  `maths_summary.json`.
- `epft_up/tests/test_sdp_theory.py`: 9 tests. Full suite: 83 passed.
- Report §5 written.

Two fixes along the way:
- The first run's "PACE correlated" model (purely Gaussian-correlated) gave
  an unphysical d_s(δRrs'') > d_s(δRrs). A pure smooth covariance is
  singular under differencing. It is now smooth plus 10% white.
- Instrument noise is now passed through the paper's 5 nm mean, and the
  PACE σ is rescaled to 1 nm. ocpy's PACE_error.csv is at ≈2 nm spacing, so
  the rescaling factor is √2, not √2.5; one of my tests had assumed 2.5.

Findings:
- (i) w_eff = Dᵀ A holds to 1e-14. The PCR-on-δRrs'' effective weights are
  99.7–99.8% high-frequency (periods < 10 nm) and uncorrelated (r ≈ 0.002)
  with the PCR-on-δRrs weights, which are ~1% high-frequency. Random-split
  skill is comparable or better on δRrs (HexFuco 0.63 vs 0.48). The
  implied prior Dᵀ S⁻² D is blind to offset and tilt and anti-smooth
  (≈5300× more prior power at f ≥ 0.3 than at f < 0.05; ≈1 without D), with
  red prior variance 7.4× blue.
- (ii) Per-band median SNR of δRrs'': ≈2 at deposit rounding, 0.04 at 2%
  in-situ noise, 0.02 at PACE white noise. Degrees of freedom:
  d_s(δRrs'') ≤ d_s(δRrs) always (a loss of 0.1–2). At 1/2.5/5/10 nm:
  - in-situ-like: 10.4, 8.6, 7.5, 6.1;
  - PACE white: 6.0, 5.3, 4.5, 3.5;
  - PACE smooth: 24, 16, 11, 8;
  - deposit rounding: 118, 93, 59, 31 (sample-limited).

  Prediction noise from δRrs''-PCR weights is 10–15× that of δRrs-PCR
  weights under white noise. Tchla's noise SD reaches ~12× its natural
  range at the PACE white level. The deposit's effective band-to-band noise
  must be ≲ 0.1% for the published PCR to work as it does.
- (iii) Filter factors on z-scored δRrs'': PCR (k ≈ 20) is a sharp cut;
  ridge (CV) has effective dof 68–103; PLS has 5–6 components, effective
  dof 54–74 and fᵢ up to 1.55. All share the same basis and its prior, so
  changing the estimator alone shouldn't fix transfer. That sets up #7.

### 2026-10-03 (Execution #7: learned wavelength weighting)

New code:
- `epft_up/sdp/weighting.py`:
  - `RidgeGCVFitter` (closed-form GCV) and `PLSFitter` (inner k-fold CV);
  - `PenalizedRRR`: a multi-output regression on raw-unit spectra ⊕
    standardized auxiliaries, with penalty λ_n wᵀΣ_n w + λ_s‖D₂w‖² + ε‖w‖²
    on the spectral block and a ridge on the auxiliaries, then reduced-rank
    projection; (λ_n, λ_s, r) by inner CV;
  - `SharedRRRFitter`: a Benchmark adaptor with one cached joint fit per
    training set, and inner CV that is random k-fold or **leave-campaign-out
    within the training set**.
- `scripts/sdp/learned_weighting.py`: 9 models through the Benchmark; rank
  and hyperparameter tallies; a full-data shared W (group inner) for the
  figures; band importance. About 25 min.
- `epft_up/tests/test_sdp_weighting.py`: 6 tests (rank recovery, smoothing,
  noise down-weighting, caching). Full suite: 89 passed.
- Report §7 written.

The first run had λ_s at the top of its grid (10) and λ_n close to it. I
re-ran with λ_n up to 100 and λ_s up to 1000. The choices moved to
100–1000, but the LOCO scores changed by < 0.01 (a flat CV surface), so I
did not extend further. My interim remark that very large λ_s means
"nearly linear weights" was wrong: the fitted weights keep 10–40 nm
structure (no power below 10 nm).

Findings (LOCO, against the best Tchla null; ratios are primary):
- Paper PCR on δRrs'': ratio log-R² 0.17, absolute R² 0.19, 1 ratio beats
  the null (the Neo artifact). Ridge or PLS on δRrs'' (estimator change):
  no better. Ridge or PLS on δRrs ⊕ GSM (basis change): absolute up to
  0.29–0.33, ratios unchanged. Per-pigment smooth+noise penalty: absolute
  0.39, ratios 0.18.
- **Shared low-rank W:** ratios 0.24–0.25, absolute 0.39. With group inner
  CV, 6 ratios beat the null: Fuco 0.45, Zea 0.48, DVchla 0.47, Chlc12
  0.43, Perid 0.20, plus Neo (artifact). 10 of 12 beat a constant out of
  campaign. Absolute ButFuco, Zea and Neo beat the null; Tchla never does.
  M1 ⊕ GSM performs about the same (Fuco:Tchla 0.55, Zea:Tchla 0.51).
- **Rank:** random inner CV chooses 5–7 (≈ the 4–5 group ceiling + local);
  group inner chooses **2 for ratios (91/108 splits)** and 3 for absolute
  (2 + Tchla). Latent 1 = diatom/Chlc12 vs Zea/DVchla (≈ Kramer & Siegel
  2019 EOF mode 1); latent 2 = Neo/Allo/Viola/HexFuco (≈ modes 2–3). Only
  about two compositional axes transfer between campaigns.
- **Wavelengths:** the shared-W weights have 0% power below 10 nm periods
  (PCR 99.7%) and ~50% at 10–40 nm. For the Fuco and cyanobacterial ratios,
  ≈ 46% of the weight contribution is at 600–700 nm (chlorophyll red
  absorption/fluorescence near 675–683; the chl c / phycobilin region near
  620–655) and only ≈ 8% at 400–450 nm.
- Gains over PCR survive LOCO for Fuco, Zea, DVchla and Chlc12 ratios and
  for absolute skill overall. They come from the smooth noise-aware prior,
  sharing across pigments and transfer-aware hyperparameter selection, not
  from the derivative.

### 2026-10-03 (Execution #8: Bayesian shared W, noise propagation, calibration)

New code:
- `epft_up/sdp/bayes.py`: `BayesianSharedW`, i.e. #7's shared W (inner
  leave-campaign-out) plus a per-target `BayesianRidge`. Default
  `calibration='oof'`: cross-fitted out-of-campaign predictions within the
  training set, then y = a + b·ŷ_oof. `'insample'` is kept for contrast.
  Also `gaussian_scores` (coverage, PIT, closed-form CRPS, SD of z).
- `epft_up/sdp/noise.py`: `band_sigma`, `covariance`, `draw` for 'insitu',
  'pace_white' and 'pace_corr', on any grid (5 nm mean at 1 nm, identity at
  5 nm).
- `epft_up/sdp/weighting.py`: `PenalizedRRR(Sigma_abs=...)` adds the untuned
  n·wᵀΣ_test w (errors-in-variables) term.
- `epft_up/sdp/gsm.py`: a `_safe_ratio` guard on the η and S_dg band ratios
  (negative or zero Rrs(555) in noisy draws overflowed η). It is inactive on
  the deposit: the clean δRrs matches #2's product to 1.2e-11.
- `scripts/sdp/uncertainty.py`: 6 scenarios × K=30 Monte Carlo draws
  through the full GSM chain; LOCO; shared W (unaware, unaware-insample,
  aware), PCR and nulls on the same noisy inputs; figures, CSVs, JSON.
  About 2 min.
- Tests: `test_sdp_bayes.py` (6), plus 1 GSM guard test and 1 Sigma_abs
  test. Full suite: 98 passed.
- Report §8 written.

Findings:
- In-sample Bayesian calibration is over-confident (deposit LOCO coverage
  64/90% of 68/95%; synthetic 48/80%), because the latent directions are fit
  to the same targets. Cross-fitted calibration fixes it: 68–71% / 93–94%
  coverage in all six scenarios, SD of z 1.05–1.10, near-flat PIT (a mild
  hump at 0.6–0.7 means a small negative bias). Without input-noise
  propagation, coverage falls to 59–63% / 88–90% under PACE noise.
  Concentrations are badly served by Gaussian intervals (SD of z 1.6–1.8);
  a log-normal or censored likelihood is the deferred item.
- **PCR on δRrs'' (the paper) collapses to LOCO R² ≈ 0** under any realistic
  noise: in-situ 2%, PACE white and PACE correlated, for ratios and absolute
  concentrations alike. It only works on deposit-clean spectra.
- **The noise-aware shared W degrades gracefully.** For the 4 transferable
  ratios, R² goes from 0.42–0.44 (clean) to 0.32–0.36 (PACE white), still
  beating the OC4-Tchla null by 0.03–0.09. Under smooth AC-like errors it is
  nearly unaffected (0.39–0.46), while the GSM-Tchla null collapses; the
  aware W then beats the null for 9/13 absolute pigments and 6/12 ratios.
  Uncertainty budget for Fuco:Tchla: model 0.27 dex, total 0.29 (in situ),
  0.33 (PACE white), 0.35 (PACE correlated).
- **5 nm sampling costs nothing** (clean and PACE white alike), in line with
  #6's degrees-of-freedom result.

### 2026-10-03 (Execution #9: EXPORTS-NA hold-out; partial, STOP-AND-ASK)

**Stopped per the prompt's instruction: the matchups cannot be
reconstructed.**
- SeaBASS programmatic access is refused: `cgi-bin/file_search.cgi` returns
  403 (also with a browser user agent and as a POST) and 444 with the
  Earthdata `~/.netrc` credentials. The `/archive/EXPORTS/...` pages return
  an empty table filled by JavaScript; `/search/archive/...` gives 404.
- PANGAEA has no EXPORTS-NA 2021 HPLC dataset (its Kramer & Siegel 2019
  compilation predates the cruise). There is no local copy.
- The Python port's `generate_coefficients` reads an
  `HPLC_Rrs_forAli_2025.xlsx` that isn't shipped.

What was done (provisional, `scripts/sdp/exports_na_holdout.py`): a
Tchla-only hold-out on the 17 EXPORTS-NA spectra in `Rrs_pigments`'
`Kramer_rrs_testdata.mat` (Rrs, T, S, lat/lon, HPLC chl; **no sample
times**; pinned commit), with every model trained on the 145.
- These spectra are processed differently from the deposit: not rounded,
  much smoother at high frequency, one clipped to 0 at 697–700 nm.
- PCR on δRrs'' gives the best Tchla here: bias +0.02 dex, RMS 0.09 (ours
  and Kramer's coefficients agree). The GSM and OC4 nulls calibrated on the
  145 are biased −0.18 and −0.25 dex. The shared W is biased −0.17, with
  68/95% coverage 0.59/1.00. N = 17 and the range 0.53–1.15 mg m⁻³ make
  R² uninformative.
- Report §9 holds a provisional Tchla-only section marked incomplete.

Follow-up (2026-10-03), after the user's answers to the stop-and-ask: (1) the
user will download the files and asked for a HOWTO; (2) use Kramer's 17
spectra as-is and flag the processing difference.
- Found the file layout from the SeaBASS cruise page
  (`/cruise/EXPORTSNA`), which lists every archive path in plain HTML.
  - **HPLC:** UCSB/CRSEO rosette files for JC214, DY131 and DY130
    (`*_rosette_HPLC_20230202_R1.sb`), plus Bowdoin/Roesler DY131 inline and
    pump files.
  - **Radiometry:** 99 UCSB/CRSEO `exports_na_{dy131,jc214}_*_C-OPS_CAST_R0.sb`
    casts (2–29 May 2021). Kramer's test file has positions but no times,
    so these casts are needed to time-match the HPLC.
- Direct archive URLs and the "Download All" endpoint return the
  JavaScript page or HTML, not data, so downloading needs a logged-in
  browser.
- Wrote `docs/HOWTO_SeaBASS_EXPORTS_NA.md`: exactly which 3 (+2) HPLC and
  99 C-OPS files, two browser routes (File Search with wildcards, or the
  archive browser), the target folder `$OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA/`,
  quick checks, and the matching plan (position → cast time → ±2 h surface
  HPLC; stop and ask on ambiguity; cross-check against Kramer's chl).
- #9 resumes once the files are in place.

HOWTO revision 2 (2026-10-03), after the user reported that the File Search
for `exports_na_*C-OPS_CAST*` returned 0 files (their screenshot is
`context/SeaBASS/Search - SeaBASS.pdf`).
- Cause: SeaBASS's "Keyword Search Filters" box matches only affiliation,
  investigator, experiment or cruise names, not file names. My first HOWTO
  was wrong on this.
- Re-examined the EXPORTSNA cruise file list. Kramer's 17 test spectra are
  1 nm hyperspectral, so the right radiometry for timing is NASA GSFC's DY131
  **HyperSAS** L2 Rrs: 34 per-station files `*_L2_Rrs_STATION_<n>_0_R1.sb`,
  13 stations, 2–29 May 2021. Their smoothness (≈3.3 nm native,
  interpolated to 1 nm) would also explain the test spectra's smoothness.
  The UCSB C-OPS casts are multispectral and are now optional.
- Rewrote `docs/HOWTO_SeaBASS_EXPORTS_NA.md`:
  - Route 1 (recommended): the archive browser, with the table filter box,
    in three directories.
  - Route 2: the File Search with the correct settings (keyword `EXPORTSNA`,
    dates, box, Data Use Warnings = Yes, Experimental ticked, HPLC / AOP +
    Hyperspectral products).
  - Where to put the files, quick checks (including the HTML-instead-of-data
    pitfall), the revised matching plan (HyperSAS station → time → ±2 h
    surface rosette HPLC; C-OPS as fallback), and an appendix of the 37
    exact archive URLs (3 HPLC + 34 HyperSAS).

### 2026-10-03 (Execution #9, completed: EXPORTS-NA hold-out)

Matchups. The user placed the SeaBASS files in `$OS_COLOR/SeaBASS/EXPORTS`.
- New `epft_up/sdp/seabass.py`:
  - `read_sb` handles −9999 → NaN, a below-detection mask and a
    `datetime` column;
  - `index_dir` strips the 10-hex-character bundle prefix and de-duplicates
    byte-identical copies. The user's download held 34 duplicated files.
- **Learned:** the test spectra are *not* the DY131 HyperSAS Rrs (no match;
  best relative RMS 5–13%). The HOWTO's inference was wrong, and a note was
  added to it. The working key is that the test positions are HPLC positions
  and the test chl is the replicate-mean surface HPLC Tchla.
- New `scripts/sdp/exports_na_matchups.py` matches the shallowest bottle at
  ≤ 12 m, with replicates averaged:
  - within 1 km, with |ΔTchla| ≤ 6e-4;
  - 17/17 unique, maximum 0.28 km, |Δ| ≤ 1e-16;
  - 13 DY131 + 4 JC214, depths 4–9 m, quality flags all 0.
  
  Output: `$OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA/exports_na_matchups.csv`
  plus a JSON with input SHA-256s.

Hold-out. `scripts/sdp/exports_na_holdout.py` was rewritten to cover 13
pigments and 12 log ratios, with no retraining:
- models: PCR (ours and Kramer's coefficients), the GSM/OC4/constant nulls,
  and the Bayesian shared W (abs hp 0.1/100/rank 3; ratio 1/100/rank 2);
- a paired bootstrap (2000 draws) of ΔRMS against the best null per target;
- outputs: `exports_na_scores.csv`, `exports_na_vs_null.csv`,
  `exports_na_tchla.png`, `exports_na_ratios.png` and
  `exports_na_summary.json`.

Findings (report §9 rewritten, no longer provisional):
- The campaign is narrow: Tchla 0.53–1.15, DVchla below LOD in all 17, Zea
  ≤ 0.07.
- PCR beats the best null on 7 of 13 absolute targets: Tchla, Fuco, Chlc12,
  Chlc3, HexFuco, ButFuco and Perid. Its ratios beat the null on Fuco,
  Chlc12, Chlc3 and Zea.
  - The gain is almost all campaign-level bias. The nulls under-predict
    Fuco by 0.6 dex.
  - PCR is worse by 0.5–1.1 dex on the trace pigments (Allo, MVchlb, Neo,
    Viola, DVchla).
- Kramer's coefficients and ours agree to ≤ 0.03 dex.
- The shared W sits in between.
  - Pooled coverage is 0.76/0.90 for abs and 0.28/0.89 for ratios: the
    errors are dominated by campaign bias, which a per-sample SD cannot
    represent.
  - Perid intervals miss; the trace-pigment intervals are too wide.
- This is consistent with LOCO (§6.2): one favourable campaign drawn from a
  wide between-campaign distribution. No contradiction with Kramer+2024,
  so no stop.
- Processing difference flagged: the test spectra are off the 1e-6 grid,
  smoother than the deposit, and have 4 non-positive values.

- Added `epft_up/tests/test_sdp_seabass.py` (reader + de-duplication); suite: 100 passed.

### 2026-10-04 (Execution #10: final report)

Finished `reports/SDP_Claude_Report.md`:
- **Header:** status set to complete.
- **§1 Introduction:** motivation, what SDP claims, the six questions, the
  two-test standard (Tchla null + LOCO), and the findings in brief.
- **§2 Literature:**
  - a summary table and a paragraph per paper: Catlett & Siegel 2018; Kramer
    & Siegel 2019; Kramer, Siegel & Graff 2020; Kramer 2022; Kramer 2024; El
    Hourany & Kramer 2026;
  - Lange 2020 as a brief non-derivative comparison;
  - where this report fits.
- **§10:**
  - 10.1 conclusions (7 points plus a bottom line);
  - 10.2 limitations;
  - 10.3 recommendations for SDP users, for EPFT-UP, and what PACE OCI
    application would require;
  - 10.4 provenance and run order of the 13 scripts.

Corrections made while writing:
- §4.2 said its N=145 Table 2 was in `kramer_table2.csv`. That CSV holds the
  N=144 run; the N=145 values are `pcr_summary.json` → `table2.n145`. The
  wording is fixed.
- §9.4: I had written that Kramer+2024 report EXPORTS-NA pigment retrievals.
  Rechecked the paper: they used the 17 samples only for community
  detection, and report no retrievals. They also call the processing
  "consistent with" the 145, which §9.1 shows is not quite so. The point is
  rewritten.

Traceability. New `scripts/sdp/report_provenance.py` writes
`reports/figures/sdp/report_manifest.json`:
- section → scripts → outputs, with SHA-256;
- each JSON's recorded `script` checked;
- every referenced figure exists, every script is cited, every output is
  mapped;
- 34 headline numbers quoted in §1/§10 re-read from the outputs and compared
  to the printed precision.

Result: 34/34 ok, 0 problems. One derived number was added through the
script: ⌈11.6²⌉ = 135 PACE pixels needed to bring SDP Tchla white-noise
below the natural SD of Tchla.

Learned:
- The report numbers are consistent with the saved outputs.
- `kramer_table2.csv` (N=144) vs `pcr_summary.json` (N=145) was the only
  mislabelled pointer found.
- The test suite is unchanged (100 passed).
