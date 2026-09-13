# MOANA re-derivation after the IOPtics → EPFT-UP move

**Run:** 2026-09-13 · `reports/scripts/moana_rederive.py` · conda env `ocean14` · source imported from IOPtics `moana` @ `3aa3b6e`.

Acceptance test for the Move (`claude_prompts/moana_prompts.md`, Q14/Q18): targets (i), (ii) and (iii-a) re-run from the raw inputs under `$OS_COLOR` with the migrated `epft_up.moana` code, and every number the report prints set beside the re-derived value at the same precision. Target (iii-b) has no in-situ counts on disk.

## Verdict: **PASS** — all 35 report-printed numbers match at their printed precision. Of the 5 values recorded only in the IOPtics prompt log, 1 differ: `pace.operational.syn.mad_dlog` log 0.007 vs 0.006.

Figures: all six pixel-identical to the IOPtics originals.

## Numbers

| # | quantity | report | re-derived | report section | |
|---|---|---|---|---|---|
| 1 | `amt24.n` | 30 | 30 | §12.2 | ✓ |
| 2 | `amt24.retrained.pro` | 1.01 / 1.22 / 0.77 | 1.01 / 1.22 / 0.77 | §12.2 table, full-fit | ✓ |
| 3 | `amt24.retrained.syn` | 1.00 / 1.36 / 0.91 | 1.00 / 1.36 / 0.91 | §12.2 table, full-fit | ✓ |
| 4 | `amt24.retrained.peuk` | 1.00 / 1.19 / 0.97 | 1.00 / 1.19 / 0.97 | §12.2 table, full-fit | ✓ |
| 5 | `amt24.cv.pro` | 1.25 / 1.78 / 0.25 | 1.25 / 1.78 / 0.25 | §12.2 table, 80/20 CV | ✓ |
| 6 | `amt24.cv.syn` | 0.68 / 3.34 / -6.5 | 0.68 / 3.34 / -6.5 | §12.2 table, 80/20 CV | ✓ |
| 7 | `amt24.cv.peuk` | 0.93 / 1.56 / 0.50 | 0.93 / 1.56 / 0.50 | §12.2 table, 80/20 CV | ✓ |
| 8 | `amt24.published.peuk` | 1.03 / 1.23 / 0.94 | 1.03 / 1.23 / 0.94 | §12.2 transfer experiment | ✓ |
| 9 | `amt24.published.pro.bias` | 1.96 | 1.96 | §12.2 transfer experiment (+96 %) | ✓ |
| 10 | `amt24.published.syn.bias` | 0.61 | 0.61 | §12.2 transfer experiment (−39 %) | ✓ |
| 11 | `amt24.basis.PC1` | 0.999 | 0.999 | §12.2 basis recovery | ✓ |
| 12 | `amt24.basis.PC2` | 0.998 | 0.998 | §12.2 basis recovery | ✓ |
| 13 | `amt24.basis.PC3` | 0.96–0.98 | 0.98 | §12.2 basis recovery (PC3–5 range) | ✓ |
| 14 | `amt24.basis.PC4` | 0.96–0.98 | 0.97 | §12.2 basis recovery (PC3–5 range) | ✓ |
| 15 | `amt24.basis.PC5` | 0.96–0.98 | 0.96 | §12.2 basis recovery (PC3–5 range) | ✓ |
| 16 | `heldout.n_stations` | 90 | 90 | prompt-14 log (stations retrieved) | ✓ |
| 17 | `heldout.23.pro` | 1.20/1.59/-1.8 | 1.20/1.59/-1.8 | §12.3 table | ✓ |
| 18 | `heldout.25.pro` | 2.47/2.90/0.15 | 2.47/2.90/0.15 | §12.3 table | ✓ |
| 19 | `heldout.28.pro` | 3.17/4.07/-0.4 | 3.17/4.07/-0.4 | §12.3 table | ✓ |
| 20 | `heldout.pooled.pro` | 2.19/2.75/-0.14 | 2.19/2.75/-0.14 | §12.3 table | ✓ |
| 21 | `heldout.23.syn` | 1.26/1.77/0.71 | 1.26/1.77/0.71 | §12.3 table | ✓ |
| 22 | `heldout.25.syn` | 1.08/2.05/0.54 | 1.08/2.05/0.54 | §12.3 table | ✓ |
| 23 | `heldout.28.syn` | 0.48/3.90/-10.6 | 0.48/3.90/-10.6 | §12.3 table | ✓ |
| 24 | `heldout.pooled.syn` | 0.85/2.46/-3.0 | 0.85/2.46/-3.0 | §12.3 table | ✓ |
| 25 | `heldout.23.peuk` | 0.79/1.40/0.76 | 0.79/1.40/0.76 | §12.3 table | ✓ |
| 26 | `heldout.25.peuk` | 1.16/1.58/0.77 | 1.16/1.58/0.77 | §12.3 table | ✓ |
| 27 | `heldout.28.peuk` | 1.24/1.43/0.80 | 1.24/1.43/0.80 | §12.3 table | ✓ |
| 28 | `heldout.pooled.peuk` | 1.07/1.47/0.78 | 1.07/1.47/0.78 | §12.3 table | ✓ |
| 29 | `heldout.pooled.n` | 66–71 | 66 | §12.3 (66–71 pooled matchups) | ✓ |
| 30 | `heldout.pooled.n` | 66–71 | 71 | §12.3 (66–71 pooled matchups) | ✓ |
| 31 | `pace.operational.syn.median_dlog` | -0.005 | -0.005 | prompt-14 log (bitexact_pace, 100k) | ✓ |
| 32 | `pace.operational.syn.mad_dlog` | 0.007 | 0.006 | prompt-14 log (bitexact_pace, 100k) | ✗ |
| 33 | `pace.atbd.syn.median_dlog` | 0.084 | 0.084 | prompt-14 log (bitexact_pace, 100k) | ✓ |
| 34 | `pace.atbd.syn.mad_dlog` | 0.090 | 0.090 | prompt-14 log (bitexact_pace, 100k) | ✓ |
| 35 | `pace.verdict` | operational | operational | §7.1 SETTLED box | ✓ |
| 36 | `pace.operational.syn.frac_exact` | 0.1–0.14 | 0.12 | §7.1 (~12 % of pixels exact) | ✓ |
| 37 | `figmethod.operational.syn.median_dlog` | -0.005 | -0.005 | §7.1 SETTLED box (figure annotation) | ✓ |
| 38 | `figmethod.operational.syn.mad_dlog` | 0.007 | 0.007 | §7.1 SETTLED box (figure annotation) | ✓ |
| 39 | `figmethod.atbd.syn.median_dlog` | 0.082 | 0.082 | §7.1 SETTLED box (figure annotation) | ✓ |
| 40 | `figmethod.atbd.syn.mad_dlog` | 0.090 | 0.090 | §7.1 SETTLED box (figure annotation) | ✓ |

Triples are bias / MAE / R² (Seegers et al. 2018 log₁₀-space bias and MAE; R² on the taxon's model scale). `amt24.basis.PCn` is |cos| between our retrained loading *n* and NASA's stored PC *n*; in-order match: yes.

## Figures (regenerated here vs the IOPtics originals, decoded pixels)

| figure | pixel-identical | pixels differing | max |Δ| (0–255) |
|---|---|---|---|
| `moana_pc_loadings.png` | yes | 0.000 % | 0 |
| `moana_product_masks.png` | yes | 0.000 % | 0 |
| `moana_clipping.png` | yes | 0.000 % | 0 |
| `moana_mapping_verdict.png` | yes | 0.000 % | 0 |
| `moana_mapping_consequence.png` | yes | 0.000 % | 0 |
| `moana_heldout_skill.png` | yes | 0.000 % | 0 |

## Run record

- target (i) `validate_amt24`: 5.9 s (37-day Level-2 chain + matchup + retraining + 100 bootstrap refits), n = 30
- target (ii) `validate_heldout_cruises`: 0.03 s, 90 stations retrieved, pooled n per taxon [66, 71, 71]
- target (iii-a) `bitexact_pace`: 6.6 s on 100,000 ocean pixels (seed 0), verdict `operational`
- target (iii-a), figure method (`plot_mapping_verdict` statistics on 60,000 pixels): operational -0.0049 / MAD 0.0066; ATBD +0.0824 / MAD 0.0905
- figures `moana_report_figs.main()`: 8.2 s
- machine-readable results and the AMT24 matchup table: `/Users/xavier/Projects/Oceanography/data/Color/EPFT-UP/moana/`
