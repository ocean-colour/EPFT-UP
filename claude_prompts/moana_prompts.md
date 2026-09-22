# MOANA

## Goals

First, move the work that I've done in the IOPtics repo into this one.
Second, improve upon the original algorithm by using additional data from across the globe.
Third, publish it all to RTD with a DOI and all

Actually, I've decided to defer the 2nd Goal for now.

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

### Website and RTD

1. I wish to: (i) Improve the main README.md and (ii) build out an RTD website.  I will do the work necessary on RTD itself.  Before we do so, ask discuss with me the work needed.  Ask a set of questions in the Q&A/RTD section.  Use Fable if you can. Log your work.

2. I have answered your questions in Q&A/RTD.  Please read them and ask another round.  Use Fable if you can. Log your work.

3. I have answered your questions in Q&A/RTD.  Please read them and proceed to execute the plan.  Use Fable if you can. Log your work.

### Dove+2026

1. I have downloaded the Dove+2026 paper that uses MOANA into `papers/dove2026.pdf`.  Please read it and then comment on how compomised (if at all) the results are given what we have learned about the MOANA algorithm.  Write your findings in the Report/Dove+2026 section.  
Use Fable if you can. Log your work.

## Report

### Dove+2026

*Dove, L. A. & Freilich, M. A. (2026), "Biophysical dynamics of mesoscale
eddies: coincident observations from SWOT and PACE", Oceanography 39(1),
14–21. Read 2026-09-22 (`papers/dove2026.pdf`, 10 pp). Assessment against
`reports/MOANA_Claude_Report.md` (rev. 4); numbers below that are ours come
from `reports/scripts/dove2026_check.py`, written for this prompt.*

#### What the paper does with MOANA

Dove & Freilich combine SWOT altimetry (the provisional 0.1° MIOST L4
product, Rossby number `Ro = ζ/f` from geostrophic velocities) with PACE OCI
products over **63–75°W, 30–40°N, April–November 2024** — the Gulf Stream
and the subtropical gyre south of it, with waters north of the Stream masked
(ADT < 0.3 m). PACE fields are the OB.DAAC L3/L4 browser products; MOANA
gives *Prochlorococcus*, *Synechococcus* and picoeukaryotes in cells mL⁻¹.
Processing: ±2-day windows → 5-day composites; days with > 50 % cloud in the
region dropped; the composites are **linearly averaged onto the 0.1°
altimetry grid**; pixels are then binned by Ro (Figs 2c, 3) or split into
"cyclonic" (Ro > 0.4) vs "background" (|Ro| < 0.1) and tracked in time
(Fig. 4), with standard errors of the mean across pixels.

The MOANA-dependent claims:

1. **Picoeukaryotes and *Synechococcus* are enhanced in cyclonic eddies** in
   spring/early summer, up to 2.5× background (Fig. 4d), with their peaks
   co-timed; concentrations rise monotonically with Ro (Fig. 3c–f), and the
   strongest signals sit at the Ro extremes only SWOT resolves.
2. **Both decline from spring to fall**, converging on cyclonic/background
   ratios of ~1.25 by August–November.
3. ***Prochlorococcus* shows little differential enhancement** (ratio
   1.0–1.2), is elevated in *both* cyclonic and anticyclonic features
   (Fig. 3a), and "remains somewhat constant" through the season (Fig. 4a).
   They interpret this as *Prochlorococcus* being weakly responsive to eddy
   nutrient enrichment, with a nod to unresolved ecotypes.

They add one sentence of caveat: algorithms "may be subject to geographic
biases in their training data as well as additional biases introduced by
reliance on sea surface temperature or salinity." That sentence is correct;
what follows puts sizes on it.

#### Verdict, taxon by taxon

| claim | rests on | status | why |
|---|---|---|---|
| picoeukaryote enhancement in cyclones (2.5×), Ro-monotonic, seasonal decline | picoeuk product | **robust** | the one taxon that transfers (held-out MAE 1.47, R² 0.78, §12.3); no SST term; its coefficient assignment is undisputed, so the §7.1 mapping error does not touch it (bit-identical under both mappings); no clipping in their domain |
| *Synechococcus* enhancement, co-timed with picoeuk, seasonal decline | Syn product | **qualitatively supported; absolute values ~20 % low** | the operational mapping makes Syn a *uniform* 0.80× the as-published value (§7.1) — cancels in the ratios and time series they show, but every cells mL⁻¹ in Figs 3–4 is low by that factor; Syn's held-out transfer is the fragile one (pooled R² −3.0, collapses on AMT28, §12.3), so the *magnitude* of the 2.5× is less certain than for picoeuk |
| *Prochlorococcus*: ratio 1.0–1.2, elevated in anticyclones too, seasonally constant | Pro product | **not supported by the product** | three independent effects, each the size of the signal claimed, detailed below |

The paper's headline — eddies as transient niches that restructure the
picophytoplankton community, seen from space alongside the physics — stands
on picoeukaryotes and *Synechococcus* and is **not compromised in
substance**. The *Prochlorococcus* paragraphs are.

#### Why the *Prochlorococcus* conclusions do not hold

**(a) The SST term is anti-correlated with their eddy classes, at the size of
the effect they report.** MOANA's Pro model is
`Pro = p₀ + 770448·log₁₀(SST) + Σ pᵢUᵢ` (§4.2, §9.10). Computed from the
vendored coefficient: **dPro/dSST = 13,400 cells mL⁻¹ per °C at 25 °C**. A
cyclonic Gulf Stream cold-core eddy runs 1–2 °C below its surroundings, so
the SST term alone *removes* **14,000–28,000 cells mL⁻¹ (5–9 % of a
3 × 10⁵ retrieval) inside cyclones** and *adds* a similar amount in warm-core
anticyclones — before any optical information enters. Their finding that
Pro is elevated in anticyclones as well as cyclones, and that the cyclonic
ratio is only 1.0–1.2, is exactly the pattern the SST crutch produces on its
own; the optical enhancement in cyclones could be larger than they report,
or absent — the product cannot say. Seasonally the effect is larger still:
regional SST rises from ~21 °C in April to ~28 °C in August, and
`770448·log₁₀(28/21)` = **+96,000 cells mL⁻¹** — a third to a half of the
entire 200–400 × 10³ range of their Fig. 4a — supplied purely by the
temperature record. "Somewhat constant" Pro through the season therefore
means either that the optical component *fell* by ~10⁵ cells mL⁻¹ or that the
retrieval is climatology-dominated; Fig. 4a cannot be read as biology until
the SST term is subtracted. This is the §9.6 circularity — anyone studying
*Prochlorococcus*–temperature coupling reads back an assumption — with a
magnitude attached, and their region (cold-core vs warm-core eddies) is the
worst case for it.

**(b) The operational coefficient mapping perturbs Pro by ~21 % with a
spatially varying sign.** The shipping product places Pro's last coefficient
on U7 instead of the published U17 (§7.1, settled empirically §12.4). Under
the as-published mapping Pro differs by a **median 54 × 10³ cells mL⁻¹ (21 %
of NASA's value)** and the difference **changes sign regionally** — it is a
distortion of spatial gradients, not a scale factor. Inside a 12° × 10° box
the sign is probably uniform (negative through the subtropical gyre), but the
term is `p₇·U₇ − p₁₇·U₁₇`, a function of each pixel's spectrum, so it varies
across eddy cores, edges and background at the ~20 % level — comparable to,
and independent of, the 0–20 % cyclonic enhancement claimed.

**(c) The Pro retrieval has not been shown to track real *Prochlorococcus*
anywhere off the training cruise.** On Lange's own held-out cruises, with
*in-situ hyperspectral* reflectance (no atmospheric correction to blame), the
published model returns 1–4 × 10⁵ cells mL⁻¹ almost regardless of the counted
abundance — pooled bias +119 %, R² −0.14, the "flat" panel of
`figures/moana_heldout_skill.png` (§12.3); applied to our own AMT24 matchups
it is +96 % biased (§12.2). Dove & Freilich's Pro values span 2–4 × 10⁵
(Fig. 3a): the whole dynamic range they interpret lies inside the band the
model emits when it is not tracking anything. Their 1.0–1.2 ratio is
therefore a statement about the product's response to eddies, not
demonstrably about the organism's.

Any one of (a)–(c) would make the Pro paragraphs unsafe; together they mean
the product's *Prochlorococcus* field in this region is substantially
SST-plus-offset, and the eddy interpretation offered for it should be
withdrawn or restated as "the MOANA product shows…".

#### What we checked in their box, and what turned out *not* to be a problem

Using the 2025-07-01 daily 0.1° L4M granule we hold (different year, same
coefficients — unchanged since OCSSW T2023.31), inside 63–75°W, 30–40°N:

| | *Pro* | *Syn* | *peuk* |
|---|---|---|---|
| pixels in box | 12,000 | | |
| fill / cloud | 5,676 (47 %) | | |
| **land sentinel 254** | **64** | 64 | 64 |
| real retrievals | 6,260 | 6,260 | 6,260 |
| clipped to 0 | **1** | 1 | 0 |
| median (p10–p90), cells mL⁻¹ | 217 k (164–413 k) | 4.2 k (2.6–8.2 k) | 953 (753–2,556) |

- **Clipping (§7.3) is not their problem.** The 18 % of zero-clipped Pro
  pixels live poleward of 40° and in cold, picoeuk-rich water; their box
  stops at 40°N and masks the Stream's north side. One zero in 6,260. The
  same goes for the `int32` truncation (§7.3): irrelevant at 10³–10⁵ cells
  mL⁻¹. This is the strongest thing to say in the paper's favour: they are
  working in MOANA's home regime — the warm, oligotrophic Atlantic, Chl
  0.1–0.4 mg m⁻³, the training basin.
- **The land-254 sentinel (App. A.1) is a small, real exposure.** 64 pixels
  per day at `254 cells mL⁻¹` for all three taxa — Bermuda (32.3°N, 64.8°W,
  inside their box and south of the Stream) and the Carolina coast in the
  north-west corner — pass a `valid_min ≤ x ≤ valid_max` screen and enter a
  linear mean as absurdly low values. Regionally it is ~1 % of the Pro mean
  and < 1 % of peuk (64 × (254 − 953)/6,260 ≈ −7 cells mL⁻¹) — negligible —
  but in a sparsely populated Ro bin a handful of 254s is an outlier tail.
  Worth a one-line check on their side; not a threat to the conclusions.
- **Their absolute levels match the product.** Pro 200–400 k, peuk ~1–2 k in
  summer background agree with our medians; their summer *Synechococcus*
  (~10 k) is above our 4.2 k median for the July granule — plausibly year-
  to-year (2024 vs 2025) or version (provisional vs V3.2) — not diagnostic.

#### Two further caveats they should carry

- **Scale and compositing (§9.9, §12.4).** Syn and peuk are `10^x` models:
  the retrieval of a mean spectrum is not the mean of retrievals. NASA already
  composites L2 retrievals to L3/L4 (which is why bit-reproduction from L3M
  Rrs is impossible, §12.4); Dove & Freilich then average 5 days and then
  average again onto 0.1° altimetry pixels. The Jensen bias this introduces
  grows with sub-pixel variance, which is largest at eddy edges and high
  |Ro| — precisely the bins (Ro ±0.7–0.9, SWOT-only) where they report the
  tightest coupling and the largest concentrations. Direction: inflates the
  patchier class relative to the quieter one, i.e. cyclonic edges relative
  to background. Unquantified here; it argues for treating the *magnitude*
  of the 2.5× as an upper-ish estimate while leaving the sign intact.
- **"Validated in the Atlantic sector" is an overstatement.** The satellite
  product has never been compared with in-situ cell counts (§10 item 4, §13
  item 2 — the validation NASA itself says has not been done). What exists is
  Lange et al.'s in-situ-radiometry training and MODIS-based held-out test,
  and now our §12.3 in-situ hyperspectral held-out — which is where the Pro
  transfer failure comes from.
- (General, not MOANA-specific:) standard errors of the mean across
  spatially autocorrelated 0.1° pixels understate the uncertainty of every
  "statistically significant" difference in Fig. 4; the ratios are the
  robust quantities.

#### Bottom line

The SWOT–PACE demonstration and the picoeukaryote / *Synechococcus* eddy
results survive what we know about MOANA: relative, same-region, same-version
comparisons of the taxa whose retrieval works, in the algorithm's home
regime, with the systematic Syn bias cancelling in ratios. The
***Prochlorococcus* results do not survive**: the SST term moves Pro by
5–9 % per eddy and ~10⁵ cells mL⁻¹ per season in the direction that mimics
their findings, the shipping coefficient mapping adds a ~21 % spectrally
varying perturbation, and the model has not been shown to track Pro off its
training cruise. The paper is a clean, real-world instance of §9.6/§9.10 —
the product's temperature crutch being read back as ecology — published
within two years of launch, and the best single argument yet for EPFT-UP's
Rrs-only Pro variant and per-retrieval uncertainty (§9.3, §9.6).

**What we could do about it (an afternoon each, data in hand):** run
`run_moana` on the L3M AOP Rrs over their box under both mappings and with
the SST term zeroed, to give them (i) the Pro field with the climatology
removed, (ii) the mapping-induced spread, and (iii) picoeuk/Syn as computed
from Rrs directly with no compositing beyond L3M. Candidate Open item for the
report's §13; also worth a courteous note to the authors, since their own
caveat sentence shows they anticipated the problem without a way to size it.

## Q&A

### RTD

*(Questions from prompt "Website and RTD"/1, 2026-09-13, run as Fable. Same
convention as Q&A/Move: each carries my recommended answer; overwrite or
correct in place. R1–R2 are the README, R3–R8 the site, R9 housekeeping.)*

**R1. README — who is it for, and what leads?**
Today's README (my Move/4 edit) is a repo map: install, a MOANA paragraph,
layout, license. IOPtics' README leads with a one-line pitch, a "What it
does" list, install, status, authors. For a public repo whose first visitors
will be ocean-colour scientists arriving from a paper, a talk or the RTD
site, the README should answer "what is this and is it any good?" before
"where are the files?".
*Recommended shape:* badges (CI · RTD · license · DOI once minted) → one-
paragraph pitch → **"What it does"** (three bullets: the MOANA baseline with
its headline numbers — basis recovered at |cos| 0.999/0.998, Lange Table 1
reproduced to within the CTD-only handicap, first held-out hyperspectral
skill numbers — then uncertainty/provenance and global training as the
roadmap) → a **5-line quickstart** (`from epft_up.moana import run_moana`
on a spectrum) → install + `$OS_COLOR` note → documentation link → status →
citing (points at `CITATION.cff`, R8) → authors → license. One figure
(`moana_heldout_skill.png`) inline, as IOPtics does not do but GitHub
renders well. Layout section dropped to the docs site. Authors line as in
IOPtics: "J. Xavier Prochaska (UC Santa Cruz) and Claude"?
>A. The README should be for the ocean optics community, i.e. semi-experts.  
I like your ideas overall. But MOANA is only the starting point, and not 
anything original (yet).  So, be sure to indicate that the future will bring
a series of empirical models for PFT analysis.

**R2. How much of the roadmap goes public?**
The three Goals at the top of this doc (baseline → X-MOANA on global data →
RTD + DOI) and the report's §13 open items are candid about what is not done
(underway FCM blocked at PML, SeaBASS credentials pending).
*Recommended:* a short "Status and roadmap" section in the README that names
the three goals and links the report's §13 for detail, without the PML
correspondence itself (R4). Say if you would rather the README stay silent
on X-MOANA until there is code.
>A. No, don't bother.  Do link to the MOANA Report but no need to provide details.

**R3. Tooling — copy the IOPtics pattern, plus MyST for the Markdown.**
IOPtics' site (live at ioptics.readthedocs.io, Sphinx + **furo** + autodoc/
napoleon, `.readthedocs.yaml` installing the package with `pip install .`,
`docs/requirements.txt` = sphinx + furo) is a working template and I would
copy it near-verbatim. The one thing IOPtics does *not* do is publish its
Markdown design docs — its site is all `.rst`. EPFT-UP's substance is
Markdown: the 1,223-line report, the design doc, `moana_blocked.md`,
`moana_rederivation.md`. I checked what those files use: tables, six images
(`figures/moana_*.png`, relative), blockquotes, Unicode footnote marks — **no
HTML, no math**. MyST (`myst_parser` 5.1.0, already in `ocean14`) renders
all of that natively.
*Recommended:* Sphinx + furo + `myst_parser` (+ `sphinx_copybutton`,
`sphinx_design`, both installed), so the Markdown is published **as-is**, no
conversion and no second copy to keep in sync. Theme: furo (matches IOPtics)
unless you want the classic RTD look.
>A.  I like your Recommendation.

**R4. What goes on the site — and what stays repo-only.**
*Recommended toctree:*
1. **Overview** (`index`) — the pitch + the held-out-skill figure;
2. **Installation & data** — install, `$OS_COLOR` layout (AMT/, PACE/,
   EPFT-UP/ for products), Earthdata `~/.netrc`;
3. **MOANA baseline** — the report, in full;
4. **Design** — the design doc;
5. **Reproduction record** — `moana_rederivation.md` (the acceptance run:
   every report number re-derived, figures byte-identical);
6. **Open items** — `moana_blocked.md`;
7. **API reference** — autodoc of the five `epft_up.moana` modules;
8. **Provenance** — where the code came from (IOPtics `moana` @ `3aa3b6e`),
   the vendored-LUT README with its sha256s, how to cite.
*Excluded, deliberately:* `requests/PML_follow_up.md` (live correspondence
naming five colleagues and an unresolved ask — repo yes, website no) and
`claude_prompts/` (the working record, including the 2,595-line history).
Both stay in git and can be linked from the site as GitHub URLs. Confirm, or
move items across the line.
>A.  This is too MOANA focused.  Again, MOANA is just one of many models we
will be developing.  Do include it, but have it be a sub-page.  Also, the
community will not have access to $OS_COLOR.  If you think it is important
that I provide the data (it probably is), then I will need to place it in Dryad
or something similar.

**R5. Where the docs source lives, and how the Markdown reaches it.**
The report's figures are regenerated into `reports/figures/` by
`moana_report_figs.py`, and the acceptance script writes
`reports/moana_rederivation.md` — those files should stay where the scripts
put them, single-sourced.
*Recommended:* Sphinx root at `docs/source/` (as IOPtics); each Markdown
document gets a two-line **shim page** there using MyST's
`` ```{include} ../../reports/MOANA_Claude_Report.md `` `` with
`:relative-images:` / `:relative-docs:` so the `figures/` paths resolve
without copying. The shim is also where a site-only banner can go — e.g. the
design doc's header still says "Status: design; no code exists yet. Prompt
12 implements this." (true when written, frozen per Q15): a one-line
admonition above the include ("historical design document, written before
implementation; see the report's §12 for what was built") fixes the reader's
expectation without touching the document. Alternative: symlinks under
`docs/source/` (git-tracked, work on RTD's Linux builders, break for some
Windows checkouts). Include-shims or symlinks?
>A.  I like your Recommendation.

**R6. API reference scope and the RTD build environment.**
`.readthedocs.yaml` installs the package (`pip install .`), so autodoc
imports the real modules; EPFT-UP's runtime deps are light (numpy, scipy,
pandas, h5py, xarray, netcdf4, earthaccess — all pip-installable), so unlike
IOPtics there is nothing unresolvable to mock. The docstrings are NumPy-style
throughout (napoleon handles them; note the `name : type — description`
dash idiom renders as the type string, fine).
*Recommended:* one page per module (`io`, `algorithm`, `pipeline`, `train`,
`validation`) with `automodule :members:`, private `_helpers` excluded;
`autodoc_mock_imports = ['earthaccess']` only (network-touching, not needed
to import); `viewcode` + `intersphinx` to numpy/scipy/pandas/xarray as in
IOPtics. `fail_on_warning: false` to start, tightened later.
>A. I like your Recommendation.

**R7. The RTD project itself — your side, but two facts you need first.**
(a) `epft-up.readthedocs.io` returns 404 today, so the slug is free.
(b) **`main` is still the initial commit** — `main cd07fd1 [behind 4]`; all
the work is on `moana` (and `start-up`). RTD builds the repository's
default branch (`main` on GitHub) unless told otherwise, so enabling it
today would publish an empty site.
*Recommended:* merge `moana` → `main` (fast-forward; nothing on `main` to
conflict) before importing the project on RTD, and keep `main` as the doc
source with `latest` tracking it; add a `stable` version when the first
release is tagged (R8). Enable pull-request previews (RTD's "Build pull
requests" toggle) so doc changes are checked before merge. Or tell me you
want RTD pointed at `moana` and I will write the config for that.
>A.  That all sounds good.

**R8. The DOI (Goal 3) — Zenodo, and what it needs from the repo.**
The standard route is Zenodo's GitHub integration: a **concept DOI** for the
software plus a DOI per tagged release, minted automatically on each GitHub
release. From the repo side it needs (a) a `CITATION.cff` (GitHub renders a
"Cite this repository" button from it; Zenodo reads it for metadata) and
optionally `.zenodo.json`; (b) a first release/tag (v0.1.0?). Neither exists
yet. The report itself could get a separate Zenodo DOI as a document, but it
is still living (§13 work queue, rev. 4).
*Recommended:* `CITATION.cff` now (authors, ORCIDs, license, repo URL,
abstract), software concept DOI at the first release — which I would time
for when the RTD site is up so the DOI badge lands on a populated README —
and a report DOI only when the report is declared frozen. I need: author
list and order, ORCIDs, and whether Claude is listed (IOPtics' README
credits "and Claude"; CFF has no good slot for a non-human author — a
`preferred-citation` note or the README line is the usual compromise).
>A.  Let's follow your Recommendation.  Authors = Prochaska and Claude.
My ORCID is 0000-0002-7738-6875.  Generate an abstract

**R9. Housekeeping I will do unless you object** (no answer needed): a
`docs` job in CI running `sphinx-build -b html` (advisory, like RTD's
`fail_on_warning: false`) so the site cannot silently break; `make html`
verified locally before handing over; a `linkcheck` pass on the migrated
docs; README badges (CI, RTD, license; DOI placeholder until R8); RTD's
`sphinx_rtd_theme`/`furo` pinned in `docs/requirements.txt`. **Not** in
scope unless you say so: dark-mode figure variants (report §13 item 10 —
furo has a dark mode and the figures have a fixed light `SURFACE`
background, so they will sit on a light card; acceptable for now), and any
edit to the frozen docs' content.
>A.  Sounds good

---

*(Round 2, from prompt "Website and RTD"/2, 2026-09-14, run as Fable. Same
convention. Your R1/R4 answers reframe the whole thing — EPFT-UP is a
*family* of empirical PFT models and MOANA is member #1 — and R4 opened the
data-access question, which is where most of this round goes. S5 is the
abstract you asked for; S6 is a verified quickstart.)*

**S1. The reframed site — does this structure read as "a framework, MOANA
inside it"?**
*Recommended toctree:*
1. **Overview** — EPFT-UP as a framework for empirical PFT models with
   uncertainty and provenance; MOANA named as the first model and the
   reproduction baseline, with a link to the report — no roadmap (R2).
2. **Getting started** — install; `Data access` (S2/S3); Earthdata login.
3. **Models** — an index page with one entry today, **MOANA**, as a
   sub-tree: *Report* · *Design* · *Reproduction record* · *Open items*
   (the four Markdown documents via include-shims, R5). Future models
   become sibling sub-trees.
4. **Data** — the datasets page: every input, its source DOI/URL, licence,
   size, and how to obtain it (S2).
5. **API reference** — `epft_up.moana` module pages (R6), organised so
   `epft_up.<next model>` slots in beside it.
6. **Provenance & citing** — code provenance, vendored-LUT sha256s,
   `CITATION.cff`, DOI badge.
The MOANA report keeps its own (first-person) title inside its page; the
shim/sidebar entry reads "MOANA report".
>A.  Yes, that looks great.

**S2. Data access for the community — what can actually be redistributed.**
I inventoried `$OS_COLOR` by source, size and licence:

| input | used by | size | source & licence | redistributable? |
|---|---|---|---|---|
| AMT23/24/25/28 flow cytometry (BODC) | targets (i), (ii) | 0.1–0.3 MB each | BODC DOIs; **NERC Open Data Licence** (attribution: "Contains data supplied by Natural Environment Research Council") | yes, with attribution — but a DOI already exists, so *cite, don't re-deposit* |
| Brewin et al. 2023 in-situ Rrs (BODC) | target (ii) | 0.8 MB | BODC DOI 10.5285/f3198e10-…; same licence | same |
| Jordan et al. 2025 underway SST netCDF | target (i) SST | 545 MB | Zenodo 10.5281/zenodo.12527954, **CC-BY-4.0** | yes — but re-downloadable by DOI, so cite |
| PACE L3M AOP + L4M MOANA granules | target (iii-a), 5 of 6 figures | 175 MB/day | NASA OB.DAAC via `earthaccess`, public (Earthdata login) | re-downloadable; code already does it |
| **AMT24 HyperSAS Level-2 (PML)** | **target (i) — the training radiometry** | **15 GB** (4.8 GB `.sav` actually read + 10 GB unused CSVs) | Brewin+ **private communication**, 2026-08; unpublished; no licence | **no — not without PML's permission** |
| vendored MOANA LUTs | everything | 25 KB | NASA OCSSW, in the repo | already shipped |

So the community can reproduce targets (ii) and (iii-a) and five figures
from public sources today; **only target (i) is gated**, by the one dataset
that is not ours to publish.
*Recommended, three parts:* (a) **don't deposit anything that has a DOI** —
instead ship `scripts/fetch_data.py` that pulls the BODC/Zenodo/PACE inputs
into the data root with checksums, which is better provenance than a copy;
(b) for the PML radiometry, **ask Brewin/Tilstone for permission** to
deposit either the 4.8 GB `.sav` subset or — far more useful and far smaller
— our **derived, screened 1-minute Rrs stream + the n = 30 matchup table**
(a few MB; still a derivative of their data, so still their call). This ask
belongs on `requests/PML_follow_up.md` as a new item; I can draft it.
(c) Repository: **Zenodo** rather than Dryad — Zenodo takes 50 GB/record for
free, mints a DOI that Zenodo↔GitHub already integrates with (R8), and does
not charge; Dryad is curated, data-only, has a fee unless the institution is
a member (UC is, via CDL), and is the better fit if you want a *curated* data
paper later. Say which, or tell me to skip deposition entirely and document
target (i) as "reproducible with the PML delivery, available on request".
>A.  I like your Recommendation, and I will ask PML.  I am confident they will give permission.

**S3. The data root — `$OS_COLOR` is ours, not the community's.**
The code resolves the root in six places (`io._os_color`,
`validation._pace_dir`, two conftest probes, two scripts), all through
`$OS_COLOR`, with the fixed sub-layout `AMT/AMT<nn>/`, `PACE/`, `EPFT-UP/`.
*Recommended:* keep the variable name `$OS_COLOR` (renaming touches the
frozen MOANA code for no functional gain) and **document it** on the Data
page as "set `$OS_COLOR` to any directory; `fetch_data.py` (S2a) creates the
layout under it". If you would rather the public-facing name be
`$EPFT_UP_DATA`, I would add it as an *alternative* read in `_os_color()`
(two lines, one test) with `$OS_COLOR` still honoured.
>A. Ok, use your Recommendation.

**S4. `CITATION.cff` — the fields I will fill.**
`title`: "EPFT-UP: Empirical Phytoplankton Functional Types with Uncertainty
and Provenance"; `authors`: Prochaska, J. Xavier (UC Santa Cruz, ORCID
0000-0002-7738-6875) and — since CFF requires a person or an *entity* —
`entity: name: "Claude (Anthropic)"`; `license: BSD-3-Clause`;
`repository-code`; `abstract` (S5); `version` and `date-released` set at
the first tag. Open: **version number for the first release** — I would
tag `v0.1.0` and bump `__version__` from `0.0.dev0` accordingly when the RTD
site is up, per R8's timing. Also whether to add `.zenodo.json` (lets
Zenodo carry keywords/communities the CFF cannot) — recommended, small.
>A. Yes, use `v0.1.0`

**S5. Abstract — draft for `CITATION.cff` / Zenodo / the Overview page
(~170 words). Edit freely.**
> EPFT-UP is an open Python framework for building, retraining and
> validating *empirical* phytoplankton functional type (PFT) algorithms —
> models that map hyperspectral remote-sensing reflectance (and ancillary
> fields such as sea-surface temperature) to the abundance or composition of
> phytoplankton groups — with per-retrieval uncertainty and end-to-end
> provenance as first-class outputs. Every product carries the data, code
> version and processing decisions that produced it, and every reported
> number is regenerated by a script in the repository. The first model in
> the family is a from-scratch reimplementation of NASA's PACE MOANA
> algorithm for *Prochlorococcus*, *Synechococcus* and picoeukaryotes: it
> recovers NASA's operational PCA basis from the AMT24 training cruise,
> reproduces the published skill within the limits of the archived data,
> reports the first held-out skill numbers for hyperspectral in-situ
> reflectance (AMT23/25/28), and documents where the shipping product
> departs from its publication. EPFT-UP is intended as a shared baseline
> against which the ocean-colour community can develop and compare the
> next generation of empirical PFT models trained on data from across the
> global ocean.
>A.  that's great, thanks

**S6. README quickstart — verified, one real PACE pixel.**
R1 asked for a five-line quickstart. A synthetic spectrum gives meaningless
abundances, so I wrote it against a real pixel (`scripts/quickstart_moana.py`,
in the repo and run this session): fetch the 2025-07-01 daily PACE Rrs
granule with `earthaccess` (cached after the first call), select the
**cloud-free pixel nearest** 30°N 60°W in the Sargasso Sea (a daily
composite is ~89 % cloud/fill — the point itself was empty and my first
version returned NaN with `FLAG_TOO_FEW_BANDS`, which is the QC doing its
job), call `run_moana` with an SST, print the result. Output this session:
`pixel 31.55°N 60.05°W · pro 150219 · syn 3337 · apeuk 893 cells/mL ·
flags 0`. It needs an Earthdata login (free) and nothing from `$OS_COLOR`
beyond a cache directory.
*Recommended:* this is the README example, with that output shown and the
sentence "Prochlorococcus needs SST; here a nominal 27 °C — the pipeline
pulls it from GHRSST in real use". Alternative: a data-free example on the
vendored LUTs only (loadings, coefficients, the standardise → project →
evaluate steps), which shows the mechanics but no science.
>A.  Use your recommendation

**S7. Where the "series of models" language goes without becoming a
roadmap.** R1 wants the README to say the future brings a series of
empirical PFT models; R2 says no roadmap details.
*Recommended wording, README and Overview alike:* "MOANA is the first of a
planned series of empirical PFT models in EPFT-UP; each will be trained,
validated and documented with the same uncertainty and provenance
machinery." One sentence, no names, no dates, no §13. The Models index page
carries the same sentence and nothing more.
>A.  Your recommendation is good

**S8. Housekeeping I will do unless you object** (no answer needed, extends
R9): the datasets page's attribution statements exactly as the licences
require (NERC's fixed sentence; CC-BY citation for Jordan; NASA OB.DAAC
acknowledgement); `scripts/quickstart_moana.py` kept in sync with the README
by a Tier-2 test that imports and runs it when the granule is cached; the
Overview figure stays `moana_heldout_skill.png` (R1) but is captioned as
"the MOANA baseline's held-out skill" so the framework, not the model, is the
subject of the page.

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

### 2026-09-13 (Website and RTD/1 — scoped the README rewrite and the RTD site; nine questions R1–R9 in Q&A/RTD)

Executed prompt "Website and RTD"/1 as **Fable 5.1**. Discussion only: nine
questions with recommended answers are in Q&A/RTD above. **No files other
than this one were changed; no git state was changed.** The Move commits
(`c1eaf8e`, `c2ac705`) are yours.

**What I checked, so the questions rest on facts rather than assumptions:**

- **IOPtics' RTD setup is a working template.** `ioptics.readthedocs.io`
  is live (HTTP 200, `/en/latest/`), built from `.readthedocs.yaml`
  (ubuntu-24.04, py3.12, `pip install .` + `docs/requirements.txt`),
  `docs/source/conf.py` (Sphinx + furo, autodoc/autosummary/napoleon/
  viewcode/intersphinx/mathjax, heavy imports mocked) and an `.rst`-only
  tree (`index`, `installation`, `datasets`, `models`, `reports/`, `api/`).
  Its Markdown design docs are **not** on the site — there is no MyST.
- **EPFT-UP's substance is Markdown**, and it is MyST-clean: the report
  (44 headings, 98 table rows, 6 images at `figures/moana_*.png`, 35
  Unicode footnote marks), design doc, `moana_blocked.md`,
  `moana_rederivation.md` — **zero HTML tags, zero `$` math** across all of
  them. `myst_parser` 5.1.0, `furo`, `sphinx_rtd_theme`, `pydata_sphinx_theme`,
  `sphinx_copybutton`, `sphinx_design` and Sphinx 9.1.0 are already in
  `ocean14`; `nbsphinx`, `numpydoc`, `sphinx_autodoc_typehints` are not (and
  are not needed).
- **No dangling cross-repo links.** After the Move's Q13 rewrites, the four
  migrated docs contain no relative links or paths outside this repo; the
  only `ioptics/` mention left is `conftest.py`'s own provenance line.
- **GitHub:** `ocean-colour/EPFT-UP` is public, BSD-3, default branch
  `main`, homepage empty. **`main` is still the initial commit**
  (`cd07fd1`, behind `moana` by 4) — RTD would build an empty site if
  enabled today (R7). `epft-up.readthedocs.io` is 404, so the slug is free.
- **No citation metadata exists** — no `CITATION.cff`, `.zenodo.json` or
  `codemeta.json`; the README has no badges (R8/R9).
- Docstrings are NumPy-style with the `name : type — description` idiom
  throughout `epft_up/moana/`, which napoleon renders (R6).

**What I learned about the repository.** The two halves of this prompt pull
in different directions and the questions are drawn along that seam: the
site wants the report and design doc published *unchanged* (they are the
frozen baseline record, Q15) while the README wants a fresh, outward-facing
voice — so the plan keeps the Markdown single-sourced under `reports/` and
`docs/design/` (where the scripts write it) and reaches it from `docs/source/`
through MyST include-shims, with any reader-facing caveats (e.g. the design
doc's pre-implementation "Status" line) living in the shim, not the document.
The other seam is public vs. repo-only: `requests/PML_follow_up.md` names
colleagues in a live data ask and `claude_prompts/` is the working record;
both belong in git, neither on the website (R4).

No git commands that change repository state were run.

### 2026-09-14 (Website and RTD/2 — answers folded in; the site reframed around a model *family*; data-access inventory; round-2 questions S1–S8, abstract drafted, quickstart written and run)

Executed prompt "Website and RTD"/2 as **Fable 5.1**. Read the R1–R9
answers, asked a second round (S1–S8 in Q&A/RTD). Two files added:
`scripts/quickstart_moana.py` (S6 — the proposed README example, executed
against the cached PACE granule) and this log. **No other files changed; no
git state changed.**

**What the answers settle.** README for semi-expert ocean-colour readers,
leading with results and a quickstart, but framing MOANA as the *starting
point* of a series of empirical PFT models (R1); no roadmap, link the report
(R2); Sphinx + furo + MyST, include-shims, autodoc per module, `main` as the
doc source after a fast-forward merge, Zenodo DOI at first release with
`CITATION.cff` now (R3, R5–R9); authors Prochaska + Claude, ORCID
0000-0002-7738-6875, abstract requested (R8). **R4 redirected the site**: not
MOANA-centric — MOANA is a sub-page under Models — and raised that the
community has no `$OS_COLOR`, with Dryad offered if the data matter.

**What I measured for the data question (S2), rather than assumed:**

- `$OS_COLOR` by source: BODC flow cytometry 0.1–0.3 MB per cruise; Brewin
  2023 Rrs 0.8 MB; Jordan 2025 netCDF 545 MB; PACE granule pair 175 MB;
  **AMT24 HyperSAS Level-2 15 GB, of which the `.sav` files the code reads
  are 4.8 GB and the never-read CSVs 10 GB**; the derived products written so
  far, 32 KB.
- Licences: the BODC deposits carry the NERC Open Data Licence (mandatory
  attribution sentence extracted from the licence document on disk); the
  Jordan record is CC-BY-4.0 (checked via the Zenodo API, record 12527954,
  eight AMT netCDFs); PACE is public via Earthdata; the PML Level-2 delivery
  is a private communication (report §12.1 and its references) with no
  licence — the single non-redistributable input, and it gates only target
  (i).
- The data root is resolved through `$OS_COLOR` at six code sites, all via
  one helper in `io.py` plus `validation._pace_dir`, the conftest probes and
  two scripts — so an alternative public-facing variable is a two-line
  change if wanted (S3).

**The quickstart runs — after one correction worth recording.**
`scripts/quickstart_moana.py`: `fetch_pace_pair` (cache hit, no download) →
`run_moana` with a nominal 27 °C SST → three abundances and a zero flag
word. The first version selected the pixel *at* 30°N 60°W and returned NaN
with flag 4 (`FLAG_TOO_FEW_BANDS`): the 2025-07-01 daily 0.1° composite has
complete spectra in only 11 % of pixels (712,412 of 6.48 M), and a 2°×2°
box around that point has none. The script now takes the nearest cloud-free
pixel (31.55°N 60.05°W → Pro 150,219, Syn 3,337, picoeuk 893 cells mL⁻¹,
flags 0). It is the executable form of the README example so the README
cannot drift from working code (S8 proposes a Tier-2 test around it).

**What I learned about the repository.** R4's correction is the important
one: the Move left the repo *shaped* like a MOANA repo (one model, one
report), and I had carried that shape into the site plan. The fix is
structural, not cosmetic — a `Models/` sub-tree and `epft_up.<model>`
namespaces that make the second model a sibling, not an appendix. On data,
the useful finding is how little is actually gated: everything except the
PML radiometry is public with a DOI, so the honest community story is
"targets (ii) and (iii-a) and five of six figures reproduce from public
sources today; target (i) needs the PML delivery" — and the smallest thing
worth asking PML to release is not their 15 GB but our few-MB screened
stream and matchup table.

No git commands that change repository state were run.

### 2026-09-14 (Website and RTD/3 — plan executed: README rewritten, Sphinx/RTD site built (16 pages, local build green), CITATION.cff + .zenodo.json, fetch_data.py, quickstart test, CI docs job; version → 0.1.0)

Executed prompt "Website and RTD"/3 as **Fable 5.1**, applying R1–R9 and
S1–S8 as answered. **No git state was changed** — the tree is uncommitted
and yours (`git status` lists 13 entries). The RTD project itself, the
`moana` → `main` fast-forward and the `v0.1.0` tag/release are your side (R7,
R8); everything the repo needs for them is in place.

**What landed (13 files new, 4 edited):**

| file | what |
|---|---|
| `README.md` | rewritten (R1/R2/S6/S7): badge row (CI · RTD · licence; DOI placeholder comment) → pitch → "What it does" with the headline numbers → the S7 sentence → held-out-skill figure → verified quickstart with its output → install + `$OS_COLOR`/`fetch_data.py` → docs link → citing → authors → licence. The layout section moved to the site. |
| `.readthedocs.yaml`, `docs/requirements.txt`, `docs/Makefile`, `docs/make.bat` | copied from IOPtics @ `develop` and renamed; `myst-parser`, `sphinx-copybutton`, `sphinx-design` added to the toolchain (R3) |
| `docs/source/conf.py` | IOPtics' conf + MyST (`colon_fence`, `deflist`, `attrs_inline`, `heading_anchors = 3`), furo with the ocean palette and GitHub source links, `autodoc_mock_imports = ['earthaccess']` only (R6), `suppress_warnings` for the GitHub-written Markdown's heading levels |
| `docs/source/_static/custom.css` | IOPtics' ocean-colour accents, renamed |
| `docs/source/index.md` | Overview: the S5 abstract as prose, the S7 sentence, the held-out-skill figure captioned as "the MOANA baseline's" (S8), three toctrees (Getting started / Models / Reference) — the S1 structure |
| `docs/source/getting_started/installation.md` | install, Earthdata `~/.netrc`, `$OS_COLOR`, the quickstart with output, tests |
| `docs/source/data.md` | the S2 inventory as a table (source DOI, licence, size, how to obtain), the three attribution statements verbatim (S8), the `$OS_COLOR` layout tree, `fetch_data.py` usage; the PML radiometry marked "available on request; permission to deposit a derived subset being sought" |
| `docs/source/models/index.md`, `models/moana/index.md` | the Models index (S7 sentence, one entry) and the MOANA landing page with its four sub-pages |
| `models/moana/{report,design,reproduction,open_items}.md` | MyST include-shims (R5) onto `reports/MOANA_Claude_Report.md`, `docs/design/moana_design.md`, `reports/moana_rederivation.md`, `reports/moana_blocked.md` with `:relative-images:`; the design shim carries the "historical design document" note and the reproduction shim a one-paragraph explanation — the documents themselves untouched (Q15) |
| `docs/source/api/index.rst`, `api/moana/{io,pipeline,algorithm,train,validation}.rst` | one autodoc page per module (R6) |
| `docs/source/provenance.md` | code provenance (IOPtics `moana` @ `3aa3b6e`), the vendored-LUT README included with `:heading-offset: 2`, how to cite, pointer to the data attributions |
| `CITATION.cff`, `.zenodo.json` | S4/S5: title, v0.1.0, 2026-09-14, BSD-3-Clause, Prochaska (UCSC, ORCID 0000-0002-7738-6875) + `Claude (Anthropic)` as a CFF entity, keywords, the S5 abstract (HTML-wrapped for Zenodo) |
| `epft_up/__init__.py`, `setup.py` | `0.0.dev0` → **`0.1.0`** (S4) |
| `scripts/fetch_data.py` | S2a: downloads the Jordan netCDF (Zenodo API link) and the PACE pair (via `fetch_pace_pair`), verifies six public files against SHA-256s measured from the copies behind every report number, prints BODC DOIs + expected filenames for anything missing (BODC is browser-only); `--verify`, `--no-pace`, `--root` |
| `epft_up/tests/test_quickstart.py` | S8: Tier-2 test that imports `scripts/quickstart_moana.py`, swaps `fetch_pace_pair` for the cached paths (no network), runs it, and asserts flags 0, positive abundances and the README's pixel (31.55°N 60.05°W) |
| `.github/workflows/ci.yml` | R9: a `docs` job that builds the site exactly as RTD does (`sphinx-build -b html`, advisory) |

**Verified:**

- `sphinx-build -b html` → **build succeeded, 16 pages**, all six report
  figures copied into `_images/`, the design-page admonition present, the
  LUT README's headings offset under Provenance, the Overview figure
  resolved from `../../reports/figures/` — so the include-shim approach
  (R5) works without copying or symlinking anything. **29 warnings, all
  from docstrings** in the frozen `epft_up/moana/` code: the
  `name : type — description` idiom with wrapped continuation lines
  trips docutils ("unexpected indentation"), and `|cosine|` / `|Δlog₁₀|`
  parse as RST substitution references. The text still renders (checked
  `compare_loadings`: the `|cosine|` appears in the page), so this is
  cosmetic. I did **not** edit the docstrings — Q10 froze the code — but a
  docstring-only touch-up (blank lines, escaping the pipes) would clear all
  29 and is worth a decision; RTD's `fail_on_warning: false` means the site
  builds regardless.
- `pytest -q` → **38 passed in 9.0 s, 0 skipped** (37 + the new quickstart
  test, which ran against the cached granule).
- `fetch_data.py --verify` → 6 verified, 0 mismatches, 0 missing.
- `linkcheck` → 10 working, 8 redirected, **3 broken — all
  `github.com/…/blob/main/…` links** (`CITATION.cff`, the history file, the
  quickstart script), broken only because `main` is still the initial
  commit; they resolve when `moana` is fast-forwarded to `main` (R7). The RTD
  badge will likewise read "unknown" until the project exists.
- `CITATION.cff`, `.zenodo.json`, `.readthedocs.yaml`, `ci.yml` parse.

**Not done, by design:** the RTD project import and `main` merge (yours,
R7); the `v0.1.0` tag → Zenodo DOI → README badge (yours, R8 — the badge
slot is a comment in the README); the PML permission ask (yours, S2 — say
the word and I will draft it into `requests/PML_follow_up.md`); the
docstring warnings above.

**What I learned about the repository.** MyST include-shims turn out to be
the whole trick: five two-line files publish 1,900 lines of frozen Markdown
with their figures, and the site-only framing (the "historical design
document" note) lives in the shim, so the record and the presentation are
separable by construction. The docstring warnings are the one place the
frozen-code decision has a cost on the site, and it is a small one. The
data page ended up being the most useful new document — it is the first
place the project states, in one table, exactly what a stranger can and
cannot reproduce from public archives.

No git commands that change repository state were run.

### 2026-09-22 (Dove+2026/1 — read Dove & Freilich (2026); assessed their MOANA-based results against the report; picoeuk/Syn findings robust, Prochlorococcus findings not supported)

Executed prompt Dove+2026/1 as **Fable 5.1**. Read `papers/dove2026.pdf`
(Oceanography 39(1), 14–21; text via `pdftotext`, figures by rendering pp.
3–8), re-read report §7.1 (mapping and its cost), §7.3–7.5, §8, §9.6,
§9.9–9.11 and §12.2–12.4, and wrote the assessment into Report/Dove+2026.
One script added, `reports/scripts/dove2026_check.py`; **no other files
changed, no git state changed.**

**What the paper is.** SWOT (0.1° MIOST, Rossby number) × PACE over
63–75°W, 30–40°N, April–November 2024; MOANA's three taxa 5-day-composited,
linearly averaged onto the altimetry grid, binned by Ro and split cyclonic
(Ro > 0.4) vs background (|Ro| < 0.1). Claims: picoeuk and Syn enhanced up
to 2.5× in cyclones in spring, declining to ~1.25× by autumn; Pro weakly
enhanced (1.0–1.2×), elevated in anticyclones too, seasonally flat.

**The assessment, in one line per taxon.** Picoeukaryotes: robust — the taxon
that transfers (§12.3), no SST term, coefficient assignment undisputed so
the §7.1 mapping error is bit-identical for it, no clipping in their domain.
*Synechococcus*: qualitatively supported; every absolute value is a uniform
0.80× (operational mapping), which cancels in their ratios; the taxon's
fragile transfer makes the 2.5× magnitude less certain. *Prochlorococcus*:
not supported — (a) the SST term is anti-correlated with their eddy classes
at the size of the claimed signal, (b) the mapping perturbation is ~21 % with
spatially varying sign, (c) the model's Pro is flat on every held-out test.
The paper's headline stands on the first two taxa and is not compromised in
substance; its Pro paragraphs should be withdrawn or restated as properties
of the product.

**Calculations (script, not chat):** `dove2026_check.py` (1) counts the
2025-07-01 L4M granule inside their box — 12,000 pixels, 5,676 fill,
**64 land-254 sentinels** (Bermuda + the Carolina coast), 6,260 real,
**1 clipped Pro zero**; medians Pro 217 k, Syn 4.2 k, peuk 953 cells mL⁻¹ —
and (2) evaluates the SST term from the vendored coefficient:
dPro/dSST = 13,400 cells mL⁻¹ °C⁻¹ at 25 °C; a −1/−2 °C cold core removes
14–28 k (5–9 % of 3 × 10⁵); the April→August 21→28 °C warming adds
**+96 k cells mL⁻¹**, a third to a half of their Fig. 4a Pro range.

**What I learned.** Two things I did not expect. First, the paper works in
MOANA's best regime — the 18 % clipping (§7.3) is a > 40°N phenomenon and
they stop at 40°N, so the product defects I would have led with are absent
from their box; the damage is entirely the *model* problems (SST crutch,
mapping, Pro non-transfer), not the *product-hygiene* ones. Second, the
eddy setting is the worst case for the SST term specifically: cyclonic =
cold-core, anticyclonic = warm-core, so the climatology term produces their
Pro pattern (weak cyclonic enhancement, anticyclonic elevation, seasonal
constancy) with no biology at all. §9.6's "circularity for exactly the
science the product invites" now has a published example, two years after
launch. For EPFT-UP this is the clearest motivation yet for the Rrs-only Pro
variant and per-retrieval uncertainty — and a concrete Open item: rerun their
box from L3M Rrs under both mappings with the SST term zeroed.

No git commands that change repository state were run.
