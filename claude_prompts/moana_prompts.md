# MOANA

## Goals

First, move the work that I've done in the IOPtics repo into this one.
Second, improve upon the original algorithm by using additional data from across the globe.
Third, publish it all to RTD with a DOI and all

## Context

The `moana` branch of `IOPtics` contains all of the work I have done on this to date.  See the prompt docs, reports, new code etc.

## Prompts

### Move

1. Read all the files in Context.  Then discuss with me a plan to migrate the docs and code from `IOPtics` into this new repository.  Put your questions in Q&A/Move below and I will answer them.
Use Fable if you can. Log your work.

2. I have answered your questions in Q&A/Move.  Please read them and ask me another round if needed.  I have started a new Claude session, so you may need to re-read the Context files.
Use Fable if you can. Log your work.
 
3. I have answered your questions in Q&A/Move.  Please read them and ask me another round if needed.  I have started a new Claude session, so you may need to re-read the Context files.
Use Fable if you can. Log your work.

4. I have answered your questions in Q&A/Move.  Please read them and proceed to execute the move and regeneration.  
Use Fable if you can. Log your work.

## Q&A

### Move

*(Questions from prompt Move/1, 2026-09-12. Each carries my recommended answer;
overwrite or correct them in place.)*

**Q1. Migration mode: fresh copy or preserved git history?**
The `moana` branch is 57 commits ahead of `main`, only 14 of which touch MOANA
paths, and the messages are `mo`, `wip`, `stuck`, `whoa` — history that carries
little information. The branch also interleaves MOANA work with unrelated
IOPtics stage-05/06 work, so a `git filter-repo`/subtree extraction would be
fiddly and would drag that along.
*Recommended:* copy the files fresh, and record provenance in the log and in the
migrated docs as "imported from IOPtics `moana` @ <sha>". No history surgery.

>A. Yes, I agree with your recommendation.

**Q2. Where does the code land — `epft_up/moana/` or flattened into `epft_up/`?**
`ioptics/moana/` is a clean, self-contained subpackage (see Q4).
*Recommended:* `epft_up/moana/`, kept intact as the NASA-reproduction baseline,
with the new global-training work going in its own sibling module. Goal 2
("improve upon the original algorithm") then reads as a diff against a fixed,
tested baseline rather than an edit that destroys it.

>A. I am likely to delete the `moana` branch altogether on `IOPtics`, i.e. without merging.

**Q3. What happens to the IOPtics copy after the move?**
Two copies that both get edited is the failure mode here.
*Recommended:* EPFT-UP becomes the sole home. You delete `ioptics/moana/`,
`ioptics/data/moana/`, `ioptics/tests/test_moana.py`, `reports/MOANA_*`,
`reports/moana_blocked.md`, `reports/scripts/moana_report_figs.py`,
`docs/design/moana_design.md` and `requests/PML_follow_up.md` from the `moana`
branch before merging it to IOPtics `main`, leaving a one-line pointer to this
repo. Alternatively the `moana` branch is simply abandoned unmerged. Which?

>A. See above.

**Q4. Does any IOPtics *core* code need to come with it?**
I checked the imports: `ioptics/moana/{io,algorithm,pipeline,train,validation}.py`
import **nothing** from IOPtics outside `ioptics.moana` itself — only numpy,
scipy, pandas, h5py and the stdlib. The one coupling is
`ioptics/tests/conftest.py::needs_amt24` (a `$OS_COLOR`/AMT24 skip marker), and
`reports/scripts/moana_report_figs.py` imports `ioptics.moana` only. So the
answer is "nothing", and the migration is a rename away from being mechanical.
*Recommended:* port just the Tier-1/Tier-2 skip-marker idiom from `conftest.py`
(~40 lines) and leave the rest of IOPtics behind.
*But* — goal 3 is "publish it all to RTD with a DOI", and IOPtics has a working
Sphinx/RTD setup (`docs/source/`, `.readthedocs.yaml`, the report → `.rst`
machinery in `ioptics/report/`). Do you want that infrastructure copied now, or
built fresh when we get to goal 3?

>A. Yes, I like your recommendation

**Q5. The 2,595-line `claude_prompts/moana_prompts.md` — verbatim or distilled?**
It holds 20 executed prompts, Q&A #1–36, eight embedded Reports, and 20 log
entries. The design doc and the report both cite "Q&A #n" as their authority for
locked decisions, so those answers are load-bearing, not archive.
*Recommended:* copy it verbatim to `claude_prompts/moana_IOPtics_history.md`
(read-only record), and keep the new `moana_prompts.md` as the live doc. The
citations then still resolve.
>A. Yes, I like your recommendation

**Q6. The papers.** `papers/*.pdf` is gitignored in IOPtics, so `lange2020.pdf`
(NTRS accepted manuscript) and `moana_atbd.pdf` exist only on this disk — and
they were already lost once in a machine migration.
*Recommended:* copy both to `EPFT-UP/papers/`, gitignored the same way, and add
their URLs to the doc so a future loss is recoverable.
>A. Yes, I like your recommendation

**Q7. Data root.** IOPtics reads `$OS_COLOR/AMT24/`, `$OS_COLOR/AMT/AMT24/`,
`$OS_COLOR/PACE/` and writes derived products to `$OS_COLOR/AMT24/derived/`.
*Recommended:* keep `$OS_COLOR` exactly as is — same tree, same paths, repo stays
data-free. EPFT-UP writes its derived products under `$OS_COLOR/EPFT-UP/`.
Confirm, or name a different convention.
>A. Confirmed

**Q8. Scope — what is EPFT-UP's algorithm, relative to MOANA?**
The repo name promises *uncertainty* and *provenance*, which are precisely
§9.3 and §9.7 of the MOANA report (its two software-side recommendations), and
goal 2 promises global training data, which is §9.11.
*Recommended reading, please confirm:* EPFT-UP = MOANA's reproduction kept as a
baseline, plus (a) a new empirical PFT algorithm with per-retrieval uncertainty
and out-of-domain detection, (b) trained on more than AMT24 alone. Still three
taxa (Pro / Syn / picoeukaryotes), or is the taxa list open too?
>A. Yes, we are going to improve upon the original MOANA.  Maybe call it X-MOANA. But all of that will come after the Move.

**Q9. Does the AMT24 underway-FCM blocker travel with us?**
`reports/moana_blocked.md` and `requests/PML_follow_up.md` §1 document a live,
unresolved data ask to PML (the underway flow-cytometry samples), last chased
2026-09-10 via Tom. It caps target (i) at n = 30 instead of Lange's n = 73–78.
*Recommended:* yes — it moves here with everything else, since this repo now
owns the reproduction. But if you would rather keep the PML correspondence
threaded through IOPtics, say so and I will leave `requests/` behind.
>A. Yes, your Recommendation

---

**One process note.** The prompt says "Use Fable if you can." This session is
Opus 5, and I cannot switch my own model mid-session — start a session with
Fable selected if you want the migration itself executed by Fable. Everything
above was produced by Opus 5.

---

*(Round 2, from prompt Move/2, 2026-09-13. Same convention: each question
carries my recommended answer; overwrite or correct in place.)*

**Q10. Q2 is still open — `epft_up/moana/` or flattened into `epft_up/`?**
The answer under Q2 ("I am likely to delete the `moana` branch altogether")
answers Q3, not Q2, so the code's landing spot is undecided. Q8's answer
("improve upon the original MOANA. Maybe call it X-MOANA. But all of that will
come after the Move") argues for keeping them separable.
*Recommended:* `epft_up/moana/`, imported unchanged apart from the import
prefix, as the frozen NASA-reproduction baseline; X-MOANA later becomes a
sibling `epft_up/xmoana/` (or similar) so that "improve upon the original" is a
diff against something that still runs and still passes its 36 tests.
>A. Use a `moana` folder, i.e. `epft_up/moana/`

**Q11. Deleting the `moana` branch unmerged discards 17,410 lines of
*non-MOANA* IOPtics work — is that intended?**
I measured `main...moana`: 24,319 added lines, of which 6,909 are MOANA/PML and
**17,410 are not**. The non-MOANA remainder is stage-05/06 IOPtics work —
`ioptics/metrics.py`, `diagnostics.py`, `evaluate.py`, `report/{rst,bokeh,
leaderboard,standard,tables,figures}.py`, `docs/source/reports/*` (two published
sweep reports), `runs/first_full_run.src`, and four `claude_prompts/` stage docs.
Abandoning the branch throws that away too, and it is the very Sphinx/RTD
machinery Q4 asked about.
Also: the migrated docs will carry "imported from IOPtics `moana` @ `3aa3b6e`"
(Q1). If the branch is deleted and was never pushed, that sha becomes
unresolvable — the provenance line points at nothing.
*Recommended:* before deleting anything, (a) land the non-MOANA work on IOPtics
`main` by whatever route you prefer, and (b) create a lightweight tag —
`git tag moana-final 3aa3b6e && git push origin moana-final` — so the branch can
be deleted while the provenance sha stays reachable forever. Then delete the
branch. If you would rather simply leave the branch sitting unmerged, that also
preserves the sha and I will cite it as-is.
>A. I will deal with deleting the branch.  But please check to see if the changes from
stage-05 and 06 are on another branch. I think they are.

**Q12. Q4's second half was not answered: RTD/Sphinx infrastructure now, or at
goal 3?**
IOPtics has a working setup (`.readthedocs.yaml`, `docs/source/conf.py`,
`docs/requirements.txt`, and the report → `.rst` machinery in `ioptics/report/`).
*Recommended:* **not now.** The Move ports the docs as Markdown exactly as they
are; RTD gets built fresh at goal 3, at which point we copy IOPtics'
`.readthedocs.yaml` + `conf.py` pattern against whatever the docs have become.
Porting `ioptics/report/*` now would import the one piece of IOPtics coupling
this migration is otherwise free of — and note it is on the branch you may
delete (Q11).
>A. Not now

**Q13. Rewriting the "Q&A #n" citations.**
`moana_design.md` and `MOANA_Claude_Report.md` carry 23 `Q&A #n` citations and
4 by-name references to `claude_prompts/moana_prompts.md`, plus 19 `ioptics/...`
path references.
*Recommended:* leave the 23 numbers untouched (they index the frozen record),
rewrite the 4 by-name references to `claude_prompts/moana_IOPtics_history.md`
so each citation resolves to the file that actually contains Q&A #1–36, and
rewrite the 19 `ioptics/` paths to `epft_up/`. Mechanical, and I will list every
changed line in the log.
>A.  Ok, let's do your Recommendation.

**Q14. What counts as "the Move succeeded"?**
Baseline measured today in IOPtics: `pytest ioptics/tests/test_moana.py -q` →
**36 passed in 4.06 s, zero skips** (the `$OS_COLOR` AMT24 Tier-2 gate is
satisfied on this machine, so the data-dependent tests really ran).
*Recommended acceptance test:* (1) `pytest -q` in EPFT-UP reproduces 36 passed,
zero skips; (2) sha256 of `pca_picophyto.h5` and `picophyt.json` match their
pinned values; (3) the six `reports/figures/moana_*.png` are copied across
byte-identical. I would **not** re-run the target (i)/(ii) reproductions or
regenerate the figures as part of the Move — the figure script needs the PACE
granule and the held-out cruises, that is real compute, and nothing about it
changes under a rename. Say the word if you want a full re-derivation as the
acceptance bar instead.
>A. Yes, I want a full re-derivation as the acceptance bar.  I want to make sure 
everything is working as expected.

**Q15. After the Move, are the report and design doc living or frozen?**
Q8 defers X-MOANA until after the Move, but the two docs are currently *living*
(report rev. 4, design rev. 2).
*Recommended:* they stay the MOANA-reproduction record and receive only
provenance-header and path corrections during the Move. X-MOANA gets its own
design doc and its own report rather than editing these, so the reproduction
baseline stays citable. The report's §13 open items and §11 NASA list stay live
and move with it — those are work queues, not history.
>A.  Let's follow your Recommendation.

**Q16. Housekeeping I will do unless you object** (no answer needed):
add `papers/*.pdf` to EPFT-UP's `.gitignore` (it is missing — IOPtics has it at
line 11, and Q6 asked for the same treatment); add `earthaccess` to
`requirements.txt` and `setup.py` (needed by `validation.py`'s target-(iii) code
paths, skip-guarded); copy `reports/moana_blocked.md` and
`requests/PML_follow_up.md` verbatim (Q9); port the Tier-1/Tier-2 skip markers
from `ioptics/tests/conftest.py` as `epft_up/tests/conftest.py`, keeping the
hang guard and renaming `$IOPTICS_TEST_TIMEOUT` → `$EPFT_UP_TEST_TIMEOUT`.
>A. Those are all fine

**Process note (unchanged from round 1).** "Use Fable if you can" — this session
is again Opus 5, which cannot switch its own model. Start a session with Fable
selected if you want the migration itself executed by Fable.

---

*(Round 3, from prompt Move/3, 2026-09-13, run as Fable. Same convention:
each question carries my recommended answer; overwrite or correct in place.
Q17 is a finding you asked for, Q18–Q21 are the last decisions before the
Move can execute, Q22 is housekeeping.)*

**Q17. Q11 follow-up — yes, the stage-05/06 work is on `develop`, and
deleting `moana` loses nothing but MOANA itself.**
Measured, not recalled. `moana` and `develop` share the merge-base `eed3b32`;
the 14 commits unique to `moana` (`32006b8`…`3aa3b6e`) are the MOANA commits,
and the 5 unique to `develop` are the PR #9/#10 merges. Every stage-05/06
file I named in Q11 — `ioptics/{metrics,diagnostics,evaluate}.py`,
`ioptics/report/*`, `docs/source/reports/*`, `runs/first_full_run.src` — is
**byte-identical** between `develop` and `moana` (`git diff --numstat` = 0
lines), and `develop` = `origin/develop` (`26a9801`, pushed). The 17,410-line
alarm in Q11 was `main...moana`; measured against `develop` the non-MOANA
delta collapses to **two lines**: `papers/*.pdf` in `.gitignore` and
`earthaccess` in `requirements.txt`, both from MOANA commit `32006b8`.
The other non-MOANA files that differ (`noise.py`, `prep.py`, `test_noise.py`,
`test_prep.py`, `coding_prompts_stage06.md`) differ because `moana` is
*behind* `develop` on them — the PR #9 `impute_frac=False` fix landed after
the branch point — so `moana` holds nothing newer there.
*Recommended:* delete `moana` whenever you like; the two lines are MOANA's and
leave with it. The one remaining point from Q11 stands: `origin/moana` is
what keeps `3aa3b6e` reachable, so if you delete the remote branch too,
`git tag moana-final 3aa3b6e && git push origin moana-final` first. The
migrated docs will cite "IOPtics `moana` @ `3aa3b6e`" either way.
>A. Thanks for checking

**Q18. The full re-derivation (Q14) — protocol, outputs and pass criterion.**
Good news first: it is cheap. I timed one mid-cruise day through
`process_day` at **0.13 s** (`scripts/probe_hsas_timing.py`, saved to this
repo), so the whole 37-day Level-2 chain is ~5 s, `validate_amt24` with its
100 bootstrap refits is well under a minute, and the prompt-18 log records
the figure script at ~10 s with the granules cached — which they are
(`$OS_COLOR/PACE/moana_validation/`, 175 MB). Everything the re-derivation
needs is on disk: 37 `.sav` days, the AMT24 FCM + Jordan SST, the Brewin 2023
Rrs, the AMT23/25/28 FCM tables, both PACE granules, and `~/.netrc` (only
needed if a granule were missing). Round 2's "that is real compute" was wrong.
*Recommended protocol:* one script, `reports/scripts/moana_rederive.py`,
that runs `validate_amt24()`, `validate_heldout_cruises()` (with the three FCM
tables) and `bitexact_pace()` in EPFT-UP, then regenerates the six figures via
`moana_report_figs.py`, and writes a Markdown table setting every number
beside its counterpart in the report — §12.2 (n = 30; the three
published/retrained/CV rows per taxon), §12.3 (per-cruise + pooled), §7.1/§12.4
(median Δlog₁₀ −0.005 / +0.084 and the verdict), and the §12.2 basis
|cos| values — to `reports/moana_rederivation.md`, with the raw results as
JSON + the matchup table as parquet under `$OS_COLOR/EPFT-UP/moana/` (the
first products under the Q7 root). Target (iii-b) is excluded — no SeaBASS
counts exist on disk, so there is nothing to re-derive.
*Pass criterion:* same conda env, fixed seeds, deterministic code → I expect
every number to match the report **exactly at the printed precision** (2 dp).
Any mismatch stops the Move and gets investigated as a migration bug; I will
not edit the report's numbers to match the run. Confirm, or name a tolerance.
>A. This is a good plan.

**Q19. Do the regenerated figures replace the imported PNGs?**
Q14 originally had the six `reports/figures/moana_*.png` copied
byte-identical; with Q18 they will also be regenerated here.
*Recommended:* commit the **regenerated** ones, so what the repo ships was
made by the repo's own script, after checking they are pixel-identical to the
IOPtics originals (a byte comparison will fail on matplotlib metadata alone;
I will compare decoded pixel arrays and report any difference).
>A. Yes, use your Recommendation.

**Q20. Sequencing — one prompt or two?**
Given Q18's runtime there is no reason to split.
*Recommended:* **Move/4 does everything**: copy + rename + provenance headers
+ path rewrites (Q13) + housekeeping (Q16, Q22) + `pytest` (36 passed, 0
skipped) + the Q18 re-derivation + Q19 figures, logged as one entry with the
full file list. If any acceptance check fails I stop there, leave the tree in
place, and report — nothing gets "fixed" silently. Alternatively Move/4 =
copy + tests, Move/5 = re-derivation, if you want to inspect the tree between.
>A. Yes, use your Recommendation.

**Q21. Where the imported code sits relative to the package's future.**
Q10 fixed `epft_up/moana/`. Two small consequences to confirm:
(a) `epft_up/__init__.py` stays as it is — it does **not** re-export MOANA
(`from epft_up.moana import run_moana` is the public path), so X-MOANA later
gets an equal footing rather than inheriting a namespace;
(b) the directory layout for the docs mirrors IOPtics exactly —
`docs/design/moana_design.md`, `reports/MOANA_Claude_Report.md`,
`reports/moana_blocked.md`, `reports/figures/`, `reports/scripts/`,
`requests/PML_follow_up.md`, `papers/` — so the 19 path rewrites are a pure
`ioptics/` → `epft_up/` prefix swap and nothing else moves.
>A.  I confirm

**Q22. Housekeeping I will do unless you object** (no answer needed; extends
Q16): `epft_up/data/moana/README.md` (116 lines, the LUT provenance and
sha256s — round 1's table omitted it) moves with the two LUTs; `setup.py`
gains `package_data={'epft_up': ['data/moana/*']}` so the LUTs ship in a
non-editable install (IOPtics never declared this — it only ever ran `-e`);
`README.md` gets a short "MOANA baseline" paragraph with the provenance line
and pointers to the report and design doc; `scripts/probe_hsas_timing.py`
(added this round) flips its import to `epft_up.moana` once the package
exists. The provenance header on each migrated doc reads
"Imported from IOPtics `moana` @ `3aa3b6e` on <date>; original path
`<path>`", above the existing header, and touches nothing else (Q15).
>A. Ok

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs
### 2026-09-12 (Move/1 — read the IOPtics MOANA work; proposed a migration plan; 9 questions raised)

Executed prompt Move/1: read the Context (the `moana` branch of
`/Users/xavier/Oceanography/python/IOPtics`, head `3aa3b6e`, 2026-09-12), and
proposed a migration plan. Questions are in Q&A/Move above; no files were moved.

**What is on the branch.** 57 commits ahead of `main`, 121 files changed. The
MOANA work itself, 20 executed prompts spanning 2026-08-01 → 2026-09-10:

| Artifact | Lines | Role |
|---|---|---|
| `claude_prompts/moana_prompts.md` | 2595 | 20 prompts, Q&A #1–36, 8 Reports, 20 log entries |
| `docs/design/moana_design.md` | 351 | the *how* — retrieval, AMT24 L2→training pipeline, code layout |
| `reports/MOANA_Claude_Report.md` | 1223 | the *what and why* — 13 sections + 2 appendices |
| `reports/moana_blocked.md` | 207 | why target (i) is capped at n = 30 |
| `reports/scripts/moana_report_figs.py` | 922 | 12 report figures (`reports/figures/*.png`) |
| `requests/PML_follow_up.md` | 211 | live data asks to PML |
| `ioptics/moana/` | 1521 | `io`, `algorithm`, `pipeline`, `train`, `validation` |
| `ioptics/data/moana/` | — | `pca_picophyto.h5` + `picophyt.json`, sha256-pinned, vendored from OCSSW |
| `ioptics/tests/test_moana.py` | 558 | 31 tests |

Status per the report §10: targets (i) reproduce Lange Tables 1–2 (done, CTD-only
configuration), (ii) held-out AMT23/25/28 (done — previously unpublished numbers),
(iii-a) PACE granule / PC-mapping question (settled: NASA ships the operational
mapping) are complete; (iii-b) PACE product vs in-situ counts awaits SeaBASS
credentials. §13 lists 10 open items.

**Proposed plan (4 steps, pending the Q&A answers):**

1. *Verbatim record.* Copy `moana_prompts.md` → `claude_prompts/moana_IOPtics_history.md`
   so that every "Q&A #n" citation in the design doc and report still resolves.
2. *Docs.* `docs/design/moana_design.md`, `reports/MOANA_Claude_Report.md`,
   `reports/moana_blocked.md`, `requests/PML_follow_up.md`, and
   `reports/figures/*.png` move across unchanged apart from a provenance header
   naming the source commit and repo-relative path fixes.
3. *Code.* `ioptics/moana/` → `epft_up/moana/`, `ioptics/data/moana/` →
   `epft_up/data/moana/`, `ioptics/tests/test_moana.py` → `epft_up/tests/`,
   `reports/scripts/moana_report_figs.py` → `reports/scripts/`. Rewrite the
   `ioptics.moana` import prefix, port the `needs_amt24` / `$OS_COLOR` skip markers
   into a new `epft_up/tests/conftest.py`, add `pyarrow`/`h5py`/`earthaccess` to
   `requirements.txt`, then verify by SHA-256 of the two vendored LUTs and
   `pytest -q`.
4. *Verify and hand back.* Confirm the test count matches IOPtics' 31, then you
   decide what happens to the IOPtics copy (Q3).

**What I learned:**

- **`ioptics/moana/` is fully self-contained.** All five modules import nothing
  from IOPtics outside `ioptics.moana` — only numpy, scipy, pandas, h5py and the
  stdlib. The entire external coupling of the MOANA track to IOPtics is *one*
  symbol: `needs_amt24` from `ioptics/tests/conftest.py`. This was deliberate —
  design doc §1 records Q&A #3 as a decision that MOANA is a "separate track"
  from the IOP work, with no `AlgorithmSpec`, no BING wrapper, no shared metrics
  table. That decision is what makes this migration cheap, and it is the single
  most useful thing I found.
- The branch's other 100-odd changed files (`metrics.py` +865, `diagnostics.py`,
  `report/*`, `docs/source/`) are IOPtics stage-05/06 work that happens to share
  the branch. None of it is MOANA. It should stay in IOPtics.
- Commit messages on the branch are `mo`, `wip`, `stuck`, `whoa` — history worth
  copying files past, not surgically extracting (Q1).
- `papers/*.pdf` is gitignored; `lange2020.pdf` and `moana_atbd.pdf` exist on this
  disk only, and the design doc footnote records that they were already lost once
  in a machine migration and had to be re-fetched. Worth copying deliberately (Q6).
- `$OS_COLOR` resolves on this machine to
  `/Users/xavier/Projects/Oceanography/data/Color/` and already holds `AMT/`,
  `PACE/`, `ocssw/` — so the Tier-2 data-dependent paths should work here unchanged.
- The three goals at the top of this doc map cleanly onto the MOANA report's own
  §9 recommendations: "uncertainty and provenance" is §9.3 + §9.7, and "additional
  data from across the globe" is §9.11 (a second basin / Koji Suzuki's
  ship-of-opportunity route). The report already contains the case for this
  repository's existence.

The prompt asks for Fable; this session is Opus 5 and cannot change its own model.
Flagged in Q&A/Move rather than silently ignored.

No git commands that change repository state were run.

### 2026-09-13 (Move/2 — re-read Context in a fresh session; answers folded in; round-2 questions Q10–Q16)

Executed prompt Move/2. New session, so I re-read the Context from scratch
(`/Users/xavier/Oceanography/python/IOPtics`, branch `moana`, head `3aa3b6e`)
rather than trusting the round-1 log, then read the answers to Q1–Q9. Seven new
questions are in Q&A/Move above. **No files were moved and no git state was
changed** — Move/2 asks only for another round of questions.

**What the answers settle.** Fresh copy with a provenance line, no history
surgery (Q1); the IOPtics `moana` branch is likely abandoned unmerged, making
EPFT-UP the sole home (Q2/Q3); port only the skip-marker idiom from
`conftest.py` (Q4); `moana_prompts.md` copied verbatim to
`claude_prompts/moana_IOPtics_history.md` so the "Q&A #n" citations resolve
(Q5); both PDFs copied and gitignored (Q6); `$OS_COLOR` unchanged, EPFT-UP
derived products under `$OS_COLOR/EPFT-UP/` (Q7); improvements ("X-MOANA") come
strictly after the Move (Q8); the PML correspondence moves here (Q9).

**What I verified this session, rather than assuming:**

- **The test baseline is 36, not 31.** `conda run -n ocean14 python -m pytest
  ioptics/tests/test_moana.py -q` → **36 passed in 4.06 s, zero skips**. Zero
  skips matters: `$OS_COLOR` resolves and `$OS_COLOR/AMT/AMT24/` exists, so the
  Tier-2 `needs_amt24` tests genuinely ran. That is the acceptance number for the
  Move (Q14); round 1's "31 tests" was wrong.
- **The figure count is 6, not 12.** `moana_report_figs.py` writes exactly six
  PNGs (`moana_{clipping,heldout_skill,mapping_consequence,mapping_verdict,
  pc_loadings,product_masks}.png`); the other 11 files in `reports/figures/` are
  IOPtics GLORIA/L23 figures that stay behind. Round 1's "12" was wrong.
- **Self-containment re-confirmed by import audit.** The five `ioptics/moana/`
  modules import only numpy, scipy, pandas, h5py, pathlib and each other. The
  whole external coupling is one symbol, `needs_amt24`. `moana_report_figs.py`
  does import `ioptics.moana` — inside functions, at lines 500, 580, 755–756 —
  so the rename touches it too.
- **Nothing is cached under `$OS_COLOR`.** No `AMT24/derived/` exists (neither at
  `$OS_COLOR/AMT24/` nor at the real path `$OS_COLOR/AMT/AMT24/`). So Q7's switch
  to `$OS_COLOR/EPFT-UP/` costs nothing — there are no products to relocate or
  invalidate. Note the data tree's real layout is `$OS_COLOR/AMT/AMT24/`, with
  siblings AMT23/25/28, which is what the held-out target (ii) uses.
- **EPFT-UP's scaffolding is closer to ready than expected**, with two gaps:
  `.gitignore` has **no** `papers/*.pdf` line (IOPtics has it at line 11), and
  `requirements.txt`/`setup.py` lack `earthaccess`. `pyarrow`, `h5py`, `netcdf4`,
  `xarray` are all already there, so round 1's dependency list was overstated —
  `earthaccess` is the only real addition. `pytest.ini` (`testpaths =
  epft_up/tests`) and the CI workflow already suit the imported suite.
- **Cross-references to rewrite:** 23 `Q&A #n` citations (leave alone), 4 by-name
  references to `claude_prompts/moana_prompts.md` (repoint at the history file),
  19 `ioptics/...` path references in the four migrating docs (repoint at
  `epft_up/`).

**The one thing worth pushing back on (Q11).** Abandoning the `moana` branch is
not free. `git diff --numstat main...moana`: 24,319 added lines, **6,909
MOANA/PML and 17,410 not** — `metrics.py`, `diagnostics.py`, `evaluate.py`, the
whole `ioptics/report/` package, `docs/source/reports/*` (two published sweeps),
`runs/first_full_run.src`, four stage-05/06 prompt docs. Deleting the branch
discards all of it, and it is precisely the Sphinx/RTD machinery Q4's second
half asked about. Separately, the provenance line the migrated docs will carry
("imported from IOPtics `moana` @ `3aa3b6e`") becomes unresolvable the moment an
unpushed branch is deleted. Recommended in Q11: land the non-MOANA work first,
and tag `3aa3b6e` before deleting so the citation survives.

**What I learned about the repository.** Two round-1 numbers (31 tests, 12
figures) did not survive checking, which is the argument for the fresh re-read
this prompt asked for. The deeper point is unchanged and now measured: MOANA's
"separate track" decision (Q&A #3, design §1) means the migration is a rename
plus one skip-marker port — the risk in this move is not the code, it is the
17,410 lines of *other* work sharing the branch it is leaving.

The prompt asks for Fable; this session is Opus 5 and cannot change its own
model. Flagged in Q&A/Move rather than silently ignored.

No git commands that change repository state were run.

### 2026-09-13 (Move/3 — run as Fable; Q11 measured against `develop`; re-derivation costed at seconds; round-3 questions Q17–Q22)

Executed prompt Move/3 as **Fable 5.1** (the "Use Fable if you can" request is
met for the first time). New session, so the Context was re-read from the
IOPtics `moana` branch (head `3aa3b6e`, working tree clean) rather than from
the earlier logs: the design doc in full, the report's §10/§12/§13 and
appendices, `moana_blocked.md`, the five `ioptics/moana/` modules' entry
points, `validation.py`'s three target functions, the figure script's
`main()`, `conftest.py`, and the prompt-14/15/18 run records in the 2,595-line
history. Then the answers to Q10–Q16. Six new questions are in Q&A/Move above.
**No files were moved and no git state was changed.** One file was added to
this repo: `scripts/probe_hsas_timing.py` (the calculation below, kept on disk
per the project convention).

**What the round-2 answers settle.** `epft_up/moana/` (Q10); no RTD/Sphinx
until goal 3 (Q12); leave the 23 `Q&A #n` numbers, repoint 4 by-name
references and 19 `ioptics/` paths (Q13); **full re-derivation is the
acceptance bar** (Q14); report and design doc are frozen apart from
provenance/path edits (Q15); the Q16 housekeeping list is approved.

**Q11 — checked as asked, and the alarm shrinks to two lines.** The user's
recollection was right: the stage-05/06 work is on `develop`
(= `origin/develop` `26a9801`, pushed; also on `first-go`). Method:
`git diff --numstat moana <branch>` restricted to the stage-05/06 paths, for
every local branch, plus merge-bases. `moana` and `develop` fork at
`eed3b32`; the 14 commits unique to `moana` are exactly the MOANA commits;
`ioptics/{metrics,diagnostics,evaluate}.py`, `ioptics/report/*`,
`docs/source/reports/*` and `runs/` are byte-identical on both. Round 2's
17,410-line figure was `main...moana` — real, but already preserved
elsewhere. Against `develop`, the only non-MOANA-path content unique to
`moana` is `papers/*.pdf` in `.gitignore` and `earthaccess` in
`requirements.txt` (both from MOANA commit `32006b8`). The other differing
files (`noise.py`, `prep.py`, two tests, `coding_prompts_stage06.md`) are
cases where `moana` is *behind* `develop` — the PR #9 review fixes landed
after the branch point — so deleting `moana` cannot lose them. The
reachability point survives: `3aa3b6e` lives on `origin/moana`, so a tag is
needed only if the remote branch is deleted too.

**Q14 — the re-derivation is seconds, not "real compute".** Measured with
`scripts/probe_hsas_timing.py` (run against the IOPtics copy via
`PYTHONPATH`, since the code has not moved yet): one mid-cruise day,
24,533 raw spectra → 18,994 after geometry → 522 one-per-minute → 178 after
QC, in **0.13 s**; 37 days ≈ 5 s. The `.sav` reads are fast because the
15 GB Level-2 tree is mostly the (unused, broken-ES) CSVs. The prompt-18 log
records the six figures at ~10 s with granules cached, and both granules are
in `$OS_COLOR/PACE/moana_validation/` (173 MB AOP + 2 MB MOANA). Every input
for targets (i), (ii), (iii-a) is on disk and `~/.netrc` carries an Earthdata
entry; only (iii-b) has no data and is excluded from the bar. Round 2 argued
against re-deriving on cost grounds; that argument was wrong and is
withdrawn — Q18 proposes the protocol, outputs and an exact-match criterion.

**Smaller facts verified this session, for the executing prompt:**

- `ioptics/data/moana/` has **three** files, not two — `README.md` (116
  lines, LUT provenance) travels too. Working-tree sha256s match the pinned
  `478ace87…deb57` / `cd6b5c1c…756d4`.
- Nothing in the codebase writes to `$OS_COLOR/AMT24/derived/` despite design
  §5.7 — `process_cruise` returns in memory. So Q7's `$OS_COLOR/EPFT-UP/`
  root is currently empty of obligations; the re-derivation outputs (Q18)
  would be its first contents.
- The figure script's `DATA_DIR` is built from `ioptics/data/moana` relative
  to the repo root and its `ioptics.moana` imports sit inside functions
  (lines 500, 580, 755–756): five edit sites in that file.
- IOPtics' `setup.py` never declared `package_data`; the LUTs only worked
  because everything ran `pip install -e .`. EPFT-UP's CI does the same, but
  Q22 adds the declaration anyway.
- The `moana` branch is checked out and clean in IOPtics, so the copy source
  is unambiguous.
- Machine: 18 cores, 128 GB — irrelevant now that the run is seconds, noted
  for X-MOANA.

**What I learned about the repository.** The migration risk I flagged in
round 2 — the other work sharing the branch — was a measurement against the
wrong baseline. Measured against `develop`, the `moana` branch is exactly
what its name says, plus two lines. That, and the pipeline running in
seconds, means the Move has no remaining unknowns of scale: it is a copy, a
prefix rename, one skip-marker port, and a re-run that must reproduce every
printed number. Round-2 answers plus Q17–Q22 are everything Move/4 needs.

No git commands that change repository state were run.

### 2026-09-13 (Move/4 — the Move executed: 23 files imported from IOPtics `moana` @ `3aa3b6e`, renamed, tested (37 passed, 0 skipped) and fully re-derived — all 35 report-printed numbers and all 6 figures reproduce)

Executed prompt Move/4 as **Fable 5.1**. Every decision in Q1–Q22 was applied
as answered; nothing beyond them was done. **No git commands that change
repository state were run** — the tree below is uncommitted and yours to
review (`git status` lists it).

**Source.** IOPtics branch `moana`, head `3aa3b6e`, working tree clean
(asserted before copying). Files were exported with `git archive` from that
commit, not copied from the working tree, except the two gitignored PDFs.

**What landed (23 files imported, 5 created, 5 edited):**

| EPFT-UP path | from IOPtics | notes |
|---|---|---|
| `epft_up/moana/{__init__,algorithm,io,pipeline,train,validation}.py` | `ioptics/moana/` | import prefix only |
| `epft_up/data/moana/{pca_picophyto.h5,picophyt.json,README.md}` | `ioptics/data/moana/` | sha256 verified `478ace87…deb57`, `cd6b5c1c…756d4`; README's one code path rewritten |
| `epft_up/tests/test_moana.py` | `ioptics/tests/test_moana.py` | import prefix + conftest import |
| `epft_up/tests/conftest.py` | ported from `ioptics/tests/conftest.py` | **new** (145 lines): hang guard + `needs_data`/`needs_amt24`/`needs_netrc`; `$IOPTICS_TEST_TIMEOUT` → `$EPFT_UP_TEST_TIMEOUT`; the six ocpy-specific markers left behind |
| `docs/design/moana_design.md` | same path | provenance header + 9 path lines |
| `reports/MOANA_Claude_Report.md` | same path | provenance header + 6 path lines + 2 by-name refs |
| `reports/moana_blocked.md` | same path | provenance header + 2 path lines |
| `requests/PML_follow_up.md` | same path | provenance header + 2 path lines + 1 by-name ref |
| `reports/figures/moana_{clipping,heldout_skill,mapping_consequence,mapping_verdict,pc_loadings,product_masks}.png` | same paths | imported byte-identical, then **regenerated here** (Q19) — the regenerated files are byte-identical to the originals too |
| `reports/scripts/moana_report_figs.py` | same path | `DATA_DIR` + 4 in-function imports + 2 docstring lines |
| `claude_prompts/moana_IOPtics_history.md` | `claude_prompts/moana_prompts.md` | verbatim (2,595 lines) + provenance header; Q&A #1–36 citations now resolve here |
| `papers/{lange2020.pdf,moana_atbd.pdf}` | `papers/` (gitignored) | copied from disk; gitignored here too |
| `reports/scripts/moana_rederive.py` | — | **new** (Q18): the acceptance script |
| `reports/moana_rederivation.md` | — | **new**: its output, the side-by-side table |
| `scripts/probe_pace_subsample.py` | — | **new**: diagnostic written while resolving the one discrepancy (below) |
| `scripts/probe_hsas_timing.py` | — | import flipped to `epft_up.moana` (Q22) |
| `.gitignore`, `requirements.txt`, `setup.py`, `README.md` | — | Q16/Q22: `papers/*.pdf`; `earthaccess` (both files); `package_data` for the LUTs; a "MOANA baseline" section + layout list |

**Q13 rewrites, every changed line.** Path prefix `ioptics/` → `epft_up/`
(19 lines): design doc 6, 27, 28, 45, 46, 93, 304, 311, 312; report 10, 36,
360, 898, 1139, 1150; `moana_blocked.md` 84, 179; `PML_follow_up.md` 86, 129
(line numbers as in the source files, before the provenance header shifted
them by two). By-name `claude_prompts/moana_prompts.md` →
`claude_prompts/moana_IOPtics_history.md` (4 lines): design 9; report 780,
901; `PML_follow_up.md` 4. The 23 `Q&A #n` numbers were not touched. In code:
`algorithm.py` 23, 42; `__init__.py` 10–13, 15, 18–20, 22; `train.py` 29,
30, 39; `io.py` 6, 125; `validation.py` 35, 36, 38, 39, 111, 230;
`pipeline.py` 31, 88, 361, 382; `test_moana.py` 1, 20–24, 30, 425, 437, 522,
537; `moana_report_figs.py` 6, 51, 54, 500, 580, 755, 756; `data/moana/
README.md` 90. Two prose mentions of "IOPtics" as a repository (design 321,
report 587) were deliberately left. The provenance header on each doc is
`> **Provenance.** Imported from IOPtics `moana` @ `3aa3b6e` on 2026-09-13;
original path `<path>`.`, inserted directly under the H1; nothing else in the
frozen docs changed (Q15).

**Acceptance — tests.** `conda run -n ocean14 python -m pytest -q` →
**37 passed in 3.6 s, 0 skipped** = the 36 MOANA tests (Tier-2 included,
since `$OS_COLOR/AMT/AMT24/Radiometry/level2` is present) + `test_version`.
No `pip install -e .` was run; pytest's rootdir import resolves `epft_up`,
and the scripts put the repo root on `sys.path` themselves.

**Acceptance — full re-derivation (Q14/Q18).**
`reports/scripts/moana_rederive.py` re-ran targets (i), (ii) and (iii-a) from
the raw inputs and regenerated the six figures, ~25 s total (target (i) 5.7 s
for the 37-day chain + matchup + retraining + 100 bootstrap refits; (iii-a)
6.5 s on 100k pixels; figures 8.2 s). Result, in `reports/moana_rederivation.md`:
**all 35 numbers the report prints reproduce exactly at their printed
precision** — §12.2 (n = 30; nine bias/MAE/R² triples; the +96 %/−39 %
transfer biases; basis |cos| 0.999/0.998 and PC3–5 within 0.96–0.98, in
order), §12.3 (twelve triples incl. pooled; 66–71 matchups; 90 stations),
§7.1 (verdict `operational`; ~12 % exact; the box's four Δlog₁₀ statistics).
**All six figures are pixel-identical** (and, it turns out, byte-identical) to
the IOPtics originals. Raw results and the n = 30 AMT24 matchup table are the
first products under the Q7 root: `$OS_COLOR/EPFT-UP/moana/
{rederivation,amt24_matchups}_2026-09-13.{json,parquet}`. The script exits 0
on the report-printed criterion.

**The one discrepancy, and what it is.** The first pass showed 3 mismatches
of 36; two were my checks, one is real but outside the report:

1. *PC3 |cos| 0.983 vs "0.96–0.98"* — my range check compared at full
   precision; at the report's 2 dp, 0.98 is inside the range. Check fixed.
2. *ATBD median Δlog₁₀ +0.084 vs the report's +0.082* — provenance, not
   arithmetic. The §7.1 box's four numbers are the **annotation of
   `moana_mapping_verdict.png`** (prompt 16; `plot_mapping_verdict` draws its
   own 60k-pixel subsample with a `nasa > 0` mask and no `nasa_compat`), not
   the `bitexact_pace` call the box names (100k pixels, seed 0, `nasa_compat`).
   Replicating the figure's method gives −0.0049/0.0066 and +0.0824/0.0905 →
   exactly the box's −0.005/0.007/+0.082/0.090; `bitexact_pace` at 100k gives
   −0.0048/0.0064 and +0.0838/0.0902, matching the prompt-14 log's +0.084.
   The script now checks each against its own source. The report is right
   as printed; its parenthetical attribution to `bitexact_pace` is loose —
   noted here, not edited (Q15).
3. *`bitexact_pace` operational MAD: 0.0064 vs the prompt-14 log's "0.007"*
   — **unresolved.** The number exists only in the IOPtics log (the report
   prints the figure's 0.007, which reproduces). `validation.py` did not
   change after prompt 14 except `validate_heldout_cruises`, and the
   statistic is flat across subsample sizes (0.0063–0.0064 for 60k–183k
   pixels; `scripts/probe_pace_subsample.py`), so the log's value did not
   come from the pinned code with these inputs — most likely a pre-commit
   code state or an unrecorded call. The medians, the ATBD MAD and the
   verdict all match. I am flagging it rather than explaining it away; it
   changes nothing in the report.

**What I learned about the repository.** The Move was as mechanical as the
"separate track" decision (Q&A #3) promised: 64 rewritten lines, one ported
conftest, zero logic changes, and every scientific number came back. The
useful surprise is that matplotlib's PNGs are byte-reproducible here, so
figure provenance can be checked with `cmp`, not just by eye. The less
comfortable lesson is item 3: a number that lives only in a prompt log, from
a run whose exact code state was never pinned, cannot be reproduced later —
which is an argument for the re-derivation script itself being the record
from now on (it writes dated JSON to `$OS_COLOR/EPFT-UP/moana/`).

**Not done, by design:** target (iii-b) (no SeaBASS counts on disk); the
IOPtics side of the Move (branch deletion / tag — yours, Q17); RTD (Q12);
X-MOANA (Q8).
