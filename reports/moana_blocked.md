# What is blocking the MOANA reproduction

> **Provenance.** Imported from IOPtics `moana` @ `3aa3b6e` on 2026-09-13; original path `reports/moana_blocked.md`.

**Status as of 2026-09-12** · Author: Claude (Fable 5), for JXP
**Related:** `requests/PML_follow_up.md` §1 (the data ask) ·
`reports/MOANA_Claude_Report.md` §12.2 (target i), §13 (open items)

---

## Summary

Reproducing Lange et al. (2020) needs radiometry paired with flow-cytometry cell
counts. AMT24 collected those counts two ways — **CTD rosette casts** (twice
daily) and an **underway run sampling every ~30 minutes** across the South
Atlantic front. The BODC deposit contains **only the CTD casts**. The underway
samples are the ones that make the *Synechococcus* model work, they exist in no
public archive, and no amount of reprocessing on our side can substitute for
them. Two requests to PML have not yet produced them, for reasons that look like
a wording problem rather than a refusal.

The block is narrow: it stops the last step of one of three validation targets.
The other two are complete.

---

## 1. What target (i) requires

Target (i) is the *full* reproduction: retrain MOANA's PCA basis and its three
regressions from AMT24 data, and compare against the paper's Table 1. It is the
result that distinguishes **"we can rebuild this algorithm"** from **"we can run
NASA's coefficients"** — the latter we can already do, bit-for-bit.

Training requires matched pairs:

```
   hyperspectral Rrs spectrum   ⟷   flow-cytometry cell count
   (same place, same time, daylight)
```

We have the radiometry: PML delivered the AMT24 HyperSAS Level-2 set (37 daily
directories, ~1 M spectra). The constraint is entirely on the counts.

## 2. AMT24 collected cell counts two different ways

Lange et al. §2.1 is explicit that both were used:

1. **CTD rosette casts.** Two per day (~03–05 h and ~13–14 h ship time), bottles
   at each depth. This is what the BODC deposit holds: 814 bottles, 68 stations.
2. **Underway sampling.** From the ship's non-toxic seawater supply, roughly
   **every 30 minutes**, while crossing the front between the South Atlantic
   Gyre and temperate water — about **25–45° S**. Their Fig. 2b shows this as a
   dense run of dots where the CTD stations are sparse.

**The deposit contains only the first.** Every one of the 814 rows carries
`Gear=CTD`.

## 3. What the deposit actually contains — verified, not assumed

BODC DOI `10.5285/a2104adc-e98f-6789-e053-6c86abc0d557`
(`$OS_COLOR/AMT/AMT24/`, plus the 2026-09 re-export in
`AMT24_JR20140922_6param/`):

| property | value |
|---|---|
| rows | 814 bottles |
| stations | 68 |
| `Gear` | **`CTD` on every row** |
| `ODV_type` | `b` (bottle) on every row |
| taxon columns | 6 (Pro, Syn, picoeuk, cryptophytes, coccolithophores, nanoeukaryotes) |
| surface samples (≤ 10 m, shallowest per station) | **68** |
| sample times | bimodal, **03–06 h** (24) and **12–15 h** (24) UTC — the two daily casts |
| matched to daylight radiometry | **n = 30** |

Checks performed, so this need not be re-litigated:

- The 2026-09-07 re-delivery (`.out` + `.xlsx`) is the **same data**, not a
  superset: single sheet `botlist`, 814 × 29, and all six taxon columns
  **numerically identical** to the copy we already had (joined on `BODC_bot`,
  `np.allclose` per column).
- The workbook has **one sheet** — checked specifically because an underway
  table would plausibly have been a second one.
- The metadata report defines exactly **six** parameter codes and the file
  carries exactly those six columns. Nothing hidden, nothing unread.
- Our loader's code→taxon mapping matches the metadata document code for code
  and the values cell for cell (`epft_up/moana/io.py:FCM_CODES`).

## 4. Why the underway samples are the decisive ones

Not a nice-to-have. The reason is specific to *Synechococcus*: it lives in
**patches** at frontal boundaries. Sampling twice a day steps over the patches;
sampling every 30 minutes resolves them. A model trained without them never sees
the high-abundance tail it is supposed to predict.

**Lange et al. tested exactly this** (their Table 2). Retraining on CTD casts
alone degraded *Synechococcus*, hyperspectral:

| configuration | bias | MAE |
|---|---|---|
| CTD casts + underway (as published) | ~1.00 | **1.27** |
| CTD casts only | 0.84 (−16 %) | **1.37** |

**Our reproduction lands at 1.36.** We have reproduced their CTD-only result
almost exactly — which is simultaneously the reassuring part and the problem.
The pipeline is right; the sample set is not.

## 5. The gap, quantified

The decisive comparison is inside the front band, because that is where the
underway run happened:

| | this deposit | Lange et al. |
|---|---|---|
| samples in **25–45° S** | **16** | ~30 min apart |
| median spacing (time) | **14.9 hours** | ~0.5 hours |
| median spacing (latitude) | **1.11°** (~66 nM) | — |
| total usable training n | **30** | **73–78** |

Anyone holding the file can verify the 16 / 14.9 h row in a few lines; see §8.

## 6. Why there is no workaround

Three independent walls. Any one would be an obstacle; together they make this a
block rather than a difficulty.

**(a) The deposit structurally cannot contain them.** It is a *bottle* product —
its own metadata calls it a "Bottle Retrieval Data Report", and every row is
keyed by `BODC_bot` with `Rosette_Pos` and `Firing_Seq`. An underway sample has
no bottle, no rosette position and no firing sequence, so it is unrepresentable
in this series. **Re-requesting this DOI can never deliver them.**

**(b) Reprocessing on our side cannot reach Lange's n.** Even using all 68 CTD
surface samples, fewer than half the casts have daylight radiometry — many are
pre-dawn. **n = 30 is a ceiling set by the data**, not by our matchup window, our
screening, or any tunable choice. Widening the matchup window trades accuracy for
a handful of samples and still falls far short of 73–78.

**(c) They exist in no public archive.** Checked three times with different query
phrasings (2026-08-01, 2026-09-07, 2026-09-10): DataCite returns **exactly one**
AMT24 flow-cytometry DOI — the bottle one. Also absent from SeaBASS (the whole
`PML/AMT` tree holds only HPLC pigments and inline AC-S IOPs, on all nine
cruises) and from the published BODC series list for the cruise (all 167 are CTD
casts).

## 7. Why the ask has not landed yet

Two rounds, neither a refusal — both look like the request being read as already
satisfied.

| date | what happened | why it missed |
|---|---|---|
| 2026-09-07 | PML/BODC re-delivered the AMT24 AFC deposit as `AMT24_JR20140922_6param/` | A re-export of the same bottle data. Values identical column by column. |
| 2026-09-10 | Tom explained the column codes: *"the cell counts are the final set of columns … P700A90Z is Abundance of Synechococcus"* | Correct, and already what we do. Answers *"where are the cell counts?"*, not *"where are the **underway** cell counts?"* |

The diagnosis: **"we need the AMT24 flow-cytometry data" sounds satisfied**,
because the flow-cytometry deposit *has* been sent, twice. The bottle-vs-underway
distinction is invisible unless you already know to look for it.

The re-ask in `requests/PML_follow_up.md` §1 is therefore rewritten to open by
confirming we have the columns and read them correctly, *before* asking the real
question — and to state plainly that a spreadsheet suffices.

## 8. What would unblock it

**One file, from Glen Tarran (PML).** Specification:

| field | detail |
|---|---|
| samples | the ~30-minute underway flow-cytometry run, ~25–45° S, Sep–Nov 2014 |
| columns | UTC timestamp, latitude, longitude, and *Prochlorococcus* / *Synechococcus* / picoeukaryotes in cells mL⁻¹ |
| size | roughly 40–50 rows is enough |
| format | **plain CSV or spreadsheet — it does not need to be archived at BODC** |

That takes the reproduction from n = 30 to Lange's n = 73–78, directly tests the
1.36 → 1.27 *Synechococcus* gap, and turns the Table 1 comparison from an
inequality into an equality.

To re-verify the gap from the data on disk:

```python
from epft_up.moana.io import load_fcm, surface_fcm
s = surface_fcm(load_fcm(cruise=24))        # 68 surface samples
band = s[(s.lat <= -25) & (s.lat >= -45)]   # the front Lange sampled underway
len(band)                                    # -> 16
```

## 9. What is *not* blocked

Worth keeping in proportion — this stops the last step of one target out of
three:

- **Target (ii), held-out skill — done.** Published coefficients on AMT23/25/28
  in-situ hyperspectral Rrs, 66–71 matchups. Numbers that did not previously
  exist (report §12.3).
- **Target (iii-a), the PACE granule — done.** Settled empirically which PC
  mapping the shipping product uses, and established that exact bit-reproduction
  from L3M inputs is structurally impossible (report §7.1, §12.4).
- **Target (i) is partly done.** The **PCA basis already reproduces**: our
  retrained loadings recover NASA's stored ones *in order*, |cos| = 0.999 (PC1),
  0.998 (PC2), 0.96–0.98 (PC3–5) — from a different processing chain and less
  than half the sample count. It also settled that the training PCA was run
  **uncentred**, a convention neither the paper nor the ATBD states.
- **Two open items need no new data at all**: the out-of-domain flag experiment
  (report §9.3 / §13 item 4) and the subpixel-variability experiment
  (§9.9 / §13 item 7).

So the honest one-line status is: *the algorithm is reproduced except for the
sample set that makes its hardest taxon work, and obtaining that is a
correspondence problem rather than a technical one.*
