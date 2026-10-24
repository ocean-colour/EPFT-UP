# Email to Tom Jordan: questions on `AMT_Rrs_HPLC.csv`

**Prepared:** 2026-10-24 (Claude, for JXP), from EPFT-UP
`claude_prompts/2nd_derivative_prompts.md`, Q&A → New AMT data, Q3–Q8.
Evidence: EPFT-UP `reports/figures/sdp/amt_inspect_summary.json` and the
report's Appendix A.

---

**Subject:** AMT Rrs + HPLC compilation: thanks, and a few questions

Hi Tom,

Thanks again for the AMT Rrs + HPLC compilation. It is exactly what our
pigment-from-hyperspectral work has been missing.

We've started with a careful inspection, and it looks very good:
- radiometry and HPLC agree closely: a GSM-type inversion gives Tchla vs
  HPLC with slope 1.03 and R² 0.87; OC4 gives R² 0.88;
- the data sit on the same 400–700 nm, 1 nm grid as the Kramer et al.
  (2022) set we've been working with;
- with 8 cruises it roughly doubles both the sample size and the number of
  independent campaigns we can validate across.

Before we use it for testing and training, a few questions. Short answers
are fine, and "don't know" is useful too.

1. **Radiometry processing.**
   - Which instrument(s) and processing were used for each cruise? In
     particular the sky-glint correction (ρ) and any NIR/residual offset.
   - What time window was used to pair Rrs with HPLC?
   - Were the spectra resampled from the instruments' native grid
     (~3–4 nm?) to 1 nm, and was any spectral smoothing applied? Several
     cruises show a weak 4 nm periodicity in the second derivative, which
     looks like an interpolation signature. That matters for us because
     we work with spectral derivatives.

2. **Spectra that appear twice.** 20 spectra are each paired with two
   different HPLC samples (list below). Some are clear: replicate bottles
   (AMT25 016-2A/2B) or CTD vs underway (AMT26). Two others need your
   view:
   - **AMT30:** 37 rows have no `Original_ID` (Sample_ID 146–189). Twelve
     of them share a spectrum with a numbered station a few minutes
     earlier or later. What are these rows, and which HPLC should go with
     each spectrum?
   - **AMT29 Sample_ID 107/108** (DY110_CTD190 / CTD189) share a spectrum
     but differ by more than a factor of two in Tchla (0.073 vs 0.165).
     Is one of them mis-matched?

3. **Below-detection values.** The file has almost no zeros, and values go
   down to ~1e-4 mg m⁻³. Were values below the detection limit kept as
   measured, or replaced? Are the per-lab LODs (PML, DHI, NASA GSFC, FURG)
   available? We need to treat sub-LOD values the same way as in the
   Kramer set, which zeroes them.

4. **MVchlb / DVchlb** are only filled for AMT29 and AMT30 (NASA GSFC). Do
   you have them, or the divinyl/monovinyl split, for the other cruises?

5. **Salinity.** Sample_ID 132 and 159 (AMT30, 2023-03-01, 33.3°S 33.3°W)
   have S = 22.23, which looks like a bad value for the open South
   Atlantic. We'll substitute climatology unless you know better.

6. **AMT24** isn't included. Was that deliberate, e.g. a processing or
   matchup issue?

7. **Sharing and credit.**
   - We'd like to share the file with Billy Li, who would try it for
     training his Bayesian neural-network pigment models, within our AI
     for Ocean group only. Is that OK with you?
   - How would you like the compilation cited or credited (and the
     underlying cruise PIs and HPLC labs)?
   - Are updated versions likely? If so, we'll version our copy.

Thanks very much. Happy to share what we find, including the
per-cruise comparisons.

Best,
Xavier

---

**List for question 2: Sample_ID pairs that share one Rrs spectrum**

| Pair (Sample_ID) | Original_ID | Cruise | Time (file) | Tchla [mg m⁻³] |
|---|---|---|---|---|
| 22 / 23 | AMT25-016-2B / 016-2A | AMT25 | 2015-09-27 06:00 / 06:00 | 0.065 / 0.066 |
| 37 / 38 | AMT26-006-5A / UND-270-1022 | AMT26 | 2016-09-26 10:03 / 10:22 | 0.212 / 0.198 |
| 104 / 105 | DY110_UW169 / UW168 | AMT29 | 2019-11-01 16:48 | 0.222 / 0.221 |
| **107 / 108** | **DY110_CTD190 / CTD189** | AMT29 | 2019-11-03 13:15 | **0.073 / 0.165** |
| 110 / 111 | DY110_CTD231 / CTD230 | AMT29 | 2019-11-07 13:20 | 0.077 / 0.076 |
| 122 / 123 | DY110_UW348 / UW349 | AMT29 | 2019-11-19 14:36 | 0.530 / 0.518 |
| 124 / 125 | DY110_UW355 / UW354 | AMT29 | 2019-11-20 14:36 | 1.410 / 1.463 |
| 127 / 128 | DY110_UW360 / UW361 | AMT29 | 2019-11-21 15:48 | 0.739 / 0.754 |
| 130 / 156 | 9 / — | AMT30 | 2023-02-27 15:03 / 15:09 | 0.103 / 0.122 |
| 132 / 159 | 13 / — | AMT30 | 2023-03-01 14:01 / 14:33 | 0.071 / 0.075 |
| 134 / 163 | 17 / — | AMT30 | 2023-03-03 14:02 / 14:37 | 0.059 / 0.064 |
| 135 / 165 | 19 / — | AMT30 | 2023-03-04 14:01 / 14:48 | 0.053 / 0.055 |
| 136 / 166 | 21 / — | AMT30 | 2023-03-05 14:01 / 14:45 | 0.039 / 0.041 |
| 137 / 167 | 23 / — | AMT30 | 2023-03-06 14:02 / 14:37 | 0.030 / 0.031 |
| 138 / 168 | 24 / — | AMT30 | 2023-03-07 14:02 / 13:51 | 0.025 / 0.026 |
| 139 / 170 | 26 / — | AMT30 | 2023-03-08 14:03 / 14:35 | 0.030 / 0.032 |
| 141 / 173 | 30 / — | AMT30 | 2023-03-10 14:03 / 14:14 | 0.054 / 0.053 |
| 143 / 176 | 35 / — | AMT30 | 2023-03-13 13:59 / 13:51 | 0.177 / 0.178 |
| 144 / 179 | 37 / — | AMT30 | 2023-03-14 13:59 / 14:37 | 0.179 / 0.190 |
| 145 / 146 | 41 / — | AMT30 | 2023-03-16 13:59 / 12:57 | 0.125 / 0.097 |
