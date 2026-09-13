# IOPtics — Explore MOANA

> **Provenance.** Imported from IOPtics `moana` @ `3aa3b6e` on 2026-09-13; original path `claude_prompts/moana_prompts.md`. Verbatim, read-only record of the IOPtics MOANA track (prompts 1–20, Q&A #1–36, logs); the live prompt doc for this repository is `claude_prompts/moana_prompts.md`.

## Goal

Explore the NASA MOANA algorithm.  We will attempt to reproduce the model and tests its accuracy.

## Conventions

- `ocean14`; JXP runs git; after each step run `pytest -q` where relevant, Q&A, Log.
- Run via the env interpreter directly
  (`/home/xavier/miniforge3/envs/ocean14/bin/python …`); `conda activate` fails
  non-interactively. This is a **Tier-2** operation — it needs the `$OS_COLOR`
  data tree (L23) mounted; do **not** unset `$OS_COLOR` here.
- `report`/`plotting` need matplotlib (`Agg` is fine headless) + Bokeh + Sphinx.

## Context

- Design docs in `docs/design`
- The paper that developed the algorithm initially: `papers/lange2020.pdf`
- Papers and websites on the Atlantic Meridional Transect 24 (AMT24).  

*(Added 2026-08-01 by prompt 1. All DOIs below were resolved and verified — via
`doi.org`, or via the Crossref/DataCite APIs where the publisher blocks
automated requests. Items marked **key** are load-bearing for reproducing
MOANA.)*

### MOANA — the operational NASA PACE algorithm

**MOANA = "Multiple Ordination ANAlysis"** (ATBD title; the Earthdata product
titles render it "Multi-Ordination ANAlysis"). *Moana* is also "ocean" in
Hawaiian. It is the first operational PACE Phytoplankton Community Composition
(PCC) product to exploit OCI's full spectral resolution, and it is the Lange+2020
algorithm essentially unchanged.

- **key** — **MOANA ATBD v1.2** (16 Dec 2024), P. K. Lange, P. J. Werdell,
  I. Cetinić; eds. I. Cetinić & G. Wang. NASA Algorithm Publication Tool.
  DOI [10.5067/0AV5267G5C77](https://doi.org/10.5067/0AV5267G5C77) ·
  landing page <https://www.earthdata.nasa.gov/apt/documents/moana/v1.0> ·
  PDF <https://oceancolor.gsfc.nasa.gov/files/atbd/atbd-obdaac-phytoplankton-groups-MOANA.pdf>.
  **Downloaded to `papers/moana_atbd.pdf`.** This is the single most useful new
  document: it gives the operational recipe *and tabulates the 25 regression
  coefficients* (its p. 5), which the paper does not.
- **key** — OCSSW reference implementation `src/l2gen/get_Cpicophyt.c`
  (author Minwei Zhang, created 2023-09-29). The Ocean Color doxygen pages have
  since been retired (they 301 to the OB.DAAC landing page); a working snapshot is
  <http://web.archive.org/web/20250712185034/https://oceancolor.gsfc.nasa.gov/docs/ocssw/get__Cpicophyt_8c_source.html>.
  It reads the PCA loadings from `$OCDATAROOT/common/pca_picophyto.h5`
  (datasets `component` [nwave × npc] and `wavelength` [int, nm]) and the
  regression coefficients from `$OCDATAROOT/common/picophyt.json`
  (keys `npc`, `pro_coef`, `syn_coef`, `apeuk_coef`). *Neither of those two data
  files is published in the literature* — see Q&A #1.
- PACE data-products overview (MOANA section, POC Ivona Cetinić):
  <https://pace.oceansciences.org/data.htm>
- Products — PACE OCI Regional Mapped MOANA, cell abundances [cells mL⁻¹],
  2024-03-05 to present:
  - L3M v3.1 DOI [10.5067/PACE/OCI/L3M/MOANA/3.1](https://doi.org/10.5067/PACE/OCI/L3M/MOANA/3.1)
  - L4M v3.1 DOI [10.5067/PACE/OCI/L4M/MOANA/3.1](https://doi.org/10.5067/PACE/OCI/L4M/MOANA/3.1)
  - v3.2 and a v3.2 NRT stream are listed in the Earthdata catalog
    (`ob-cloud-pace-oci-l4m-moana-3.2`, `…-nrt-3.2`) but their DOIs were **not yet
    registered** as of 2026-07-31, so cite v3.1 for now.
  - **Corrected 2026-08-01 (prompt 3), after opening a real granule:** the CMR
    collections that actually exist are **`PACE_OCI_L4M_MOANA` v3.2** and
    `PACE_OCI_L4M_MOANA_NRT` v3.2 — there is **no** L3M MOANA collection in CMR
    despite the catalog page above. And the in-file variable names are the
    misspelled L2-style **`prococcus_moana`, `syncoccus_moana`, `picoeuk_moana`**,
    *not* the `prochlorococcus_moana` / `synechococcus_moana` the Earthdata
    catalog page advertises. Use the in-file names. The Rrs input lives in
    **`PACE_OCI_L3M_AOP` v3.2** (variable `Rrs`, 172 OCI wavelengths spanning
    346–719 nm, so 414–660 is fully covered).
- **Validation status (important):** the ATBD states plainly that *"PACE OCI MOANA
  products have not been validated yet."* Operational validation is against
  SeaBASS cell counts following Bailey & Werdell (2006). The only published skill
  numbers are Lange+2020 Tables 1–3 (AMT in-situ and Aqua-MODIS).
- Cetinić, I., et al. (2024). "Phytoplankton composition from sPACE:
  Requirements, opportunities, and challenges." *Remote Sensing of Environment*
  302, 113964. DOI [10.1016/j.rse.2023.113964](https://doi.org/10.1016/j.rse.2023.113964)
  — the PACE PCC roadmap; the ATBD points here for background.
- Mouw, C. B., et al. (2017). "A Consumer's Guide to Satellite Remote Sensing of
  Multiple Phytoplankton Groups in the Global Ocean." *Frontiers in Marine
  Science* 4, 41. DOI [10.3389/fmars.2017.00041](https://doi.org/10.3389/fmars.2017.00041)

Method lineage the ATBD credits for the principal-component-regression approach:

- Craig, S. E., et al. (2012). *Remote Sensing of Environment* 119, 72–83.
  DOI [10.1016/j.rse.2011.12.007](https://doi.org/10.1016/j.rse.2011.12.007)
- Bracher, A., et al. (2015). *Ocean Science* 11, 139–158.
  DOI [10.5194/os-11-139-2015](https://doi.org/10.5194/os-11-139-2015)
- Seegers, B. N., et al. (2018). *Optics Express* 26, 7404–7422.
  DOI [10.1364/OE.26.007404](https://doi.org/10.1364/OE.26.007404) — the
  log-space bias/MAE metrics Lange+2020 uses (already in the IOPtics design doc).
- Bailey, S. W. & Werdell, P. J. (2006). *Remote Sensing of Environment* 102,
  12–23. DOI [10.1016/j.rse.2006.01.015](https://doi.org/10.1016/j.rse.2006.01.015)
  — the satellite/in-situ matchup protocol.
- Werdell, P. J., et al. (2003). "Unique data repository facilitates ocean color
  satellite validation." *Eos* 84(38), 377–387.
  DOI [10.1029/2003EO380001](https://doi.org/10.1029/2003EO380001) — SeaBASS.

### AMT24 — the MOANA training cruise

AMT24 = RRS *James Clark Ross*, cruise **JR20140922 / JR303**, Immingham (UK) →
Stanley (Falkland Islands), 22 Sep – 2 Nov 2014 (Lange+2020 uses the
30 Sep – 1 Nov window); principal scientist Tim Smyth (PML); 70 CTD stations plus
30-minute underway sampling across the South Atlantic Gyre front (25–45°S).

- AMT programme site <https://amt-uk.org/>, cruise page
  <https://www.amt-uk.org/Cruises/AMT24>, data page <https://amt-uk.org/data/>
- BODC cruise summary report (JR20140922)
  <https://www.bodc.ac.uk/resources/inventories/cruise_inventory/report/15038/>
- **key** — **AMT24 flow-cytometric cell abundances** (the MOANA training truth):
  Tarran, G. A. & Zubkov, M. V. (2020), BODC/NOC/NERC.
  DOI [10.5285/a2104adc-e98f-6789-e053-6c86abc0d557](https://doi.org/10.5285/a2104adc-e98f-6789-e053-6c86abc0d557)
  (shortDOI `10/ds9g`). *Prochlorococcus*, *Synechococcus*, pico- and
  eukaryotic phytoplankton by FACSort, 200 m → surface.
- Flow-cytometry datasets for the five **validation** cruises used in
  Lange+2020 §2.2 (all BODC, 2020):
  | Cruise | DOI | shortDOI |
  |---|---|---|
  | AMT20 (JC053, Oct–Nov 2010) | [10.5285/a2104adc-e993-6789-e053-6c86abc0d557](https://doi.org/10.5285/a2104adc-e993-6789-e053-6c86abc0d557) | `10/dsrh` |
  | AMT22 (JC079, Oct–Nov 2012) | [10.5285/a2104adc-e991-6789-e053-6c86abc0d557](https://doi.org/10.5285/a2104adc-e991-6789-e053-6c86abc0d557) | `10/dsrm` |
  | AMT23 (JR300, Oct–Nov 2013) | [10.5285/a2104adc-e990-6789-e053-6c86abc0d557](https://doi.org/10.5285/a2104adc-e990-6789-e053-6c86abc0d557) | `10/ds7f` |
  | AMT25 (JR15001, Sep–Nov 2015) | [10.5285/a2104adc-e98e-6789-e053-6c86abc0d557](https://doi.org/10.5285/a2104adc-e98e-6789-e053-6c86abc0d557) | `10/dq4r` |
  | AMT28 (JR18001, Sep–Oct 2018) | [10.5285/a147c314-688b-55e9-e053-6c86abc0dc81](https://doi.org/10.5285/a147c314-688b-55e9-e053-6c86abc0dc81) | `10/dqwc` |
- Parent BODC collection, "Biological and Hydrographic station data collected
  during the AMT programme (1995–)", AMT13→AMT30:
  DOI [10.5285/41479c42-4dfb-4da9-be97-4c532ce13922](https://doi.org/10.5285/41479c42-4dfb-4da9-be97-4c532ce13922)
- Rees, A. P., et al. (2017). "The Atlantic Meridional Transect programme
  (1995–2016)." *Progress in Oceanography* 158, 3–18.
  DOI [10.1016/j.pocean.2017.05.004](https://doi.org/10.1016/j.pocean.2017.05.004)
- Rees, A., Robinson, C., Smyth, T., Aiken, J., et al. (2015). "20 Years of the
  Atlantic Meridional Transect — AMT." *Limnology and Oceanography Bulletin* 24,
  101–107. DOI [10.1002/lob.10069](https://doi.org/10.1002/lob.10069)

### AMT bio-optics — directly useful to IOPtics

- **key** — Jordan, T. M., Dall'Olmo, G., Tilstone, G., Brewin, R. J. W.,
  Nencioli, F., Airs, R., Thomas, C. S., Schlüter, L. (2025). "A compilation of
  surface inherent optical properties and phytoplankton pigment concentrations
  from the Atlantic Meridional Transect." *Earth System Science Data* 17,
  493–516. DOI [10.5194/essd-17-493-2025](https://doi.org/10.5194/essd-17-493-2025).
  Data: Zenodo DOI [10.5281/zenodo.12527954](https://doi.org/10.5281/zenodo.12527954)
  and SeaBASS <https://seabass.gsfc.nasa.gov/archive/PML/AMT>.
  Underway particulate absorption / scattering / attenuation (~400–720 nm,
  ~270k hyperspectral + ~40k multispectral records) plus HPLC pigments for
  **nine cruises including AMT24** — i.e. real in-situ IOPs on the very cruise
  that trained MOANA. This is a candidate new IOPtics dataset in its own right.
- Brewin, R. J. W., Dall'Olmo, G., Pardo, S., van Dongen-Vogels, V., Boss, E. S.
  (2016). "Underway spectrophotometry along the Atlantic Meridional Transect
  reveals high performance in satellite chlorophyll retrievals." *Remote Sensing
  of Environment* 183, 82–97.
  DOI [10.1016/j.rse.2016.05.005](https://doi.org/10.1016/j.rse.2016.05.005)
  — Lange+2020 refs 42/43; the source of the HyperSAS Rrs processing protocol.
- **key** — **Brewin, R. J. W., Pitarch Portero, J. S., Dall'Olmo, G., van der Woerd,
  H. J., Lin, J., Sun, X., Tilstone, G. H. (2023).** "Modern and traditional optical
  measurements, and environmental data, collected on four Atlantic Meridional
  Transect cruises between 2013 and 2018." NERC EDS British Oceanographic Data
  Centre NOC.
  DOI [10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6](https://doi.org/10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6)
  — **in-situ hyperspectral remote-sensing reflectance** plus Secchi depth,
  Forel-Ule colour, Chl-a and diffuse/beam attenuation, at 127 stations on
  **AMT23, AMT25, AMT26 and AMT28**. Found in prompt 7. Does **not** include AMT24,
  but three of these four (23, 25, 28) are Lange+2020's held-out validation cruises
  and already have flow-cytometry DOIs above — so this unblocks validation target
  (ii) with *in-situ* rather than satellite `Rrs`. Formats: Binary + Delimited;
  download is browser-only (BODC PDL).
- **key for noise** — Lin, J., Dall'Olmo, G., Tilstone, G. H., et al. (2022).
  "Derivation of uncertainty budgets for continuous above-water radiometric
  measurements along an Atlantic Meridional Transect." *Optics Express* 30,
  45648. DOI [10.1364/OE.470994](https://doi.org/10.1364/OE.470994)
  — per-band Rrs uncertainty for exactly this instrument/platform, i.e. the
  honest `varRrs` for an AMT fit instead of a flat 5%.
- Alikas, K., Vabson, V., Ansko, I., Tilstone, G. H., et al. (2020). "Comparison
  of Above-Water Seabird and TriOS Radiometers along an Atlantic Meridional
  Transect." *Remote Sensing* 12, 1669.
  DOI [10.3390/rs12101669](https://doi.org/10.3390/rs12101669)
- Pardo, S., Tilstone, G. H., Brewin, R. J. W., et al. (2023). "Radiometric
  assessment of OLCI, VIIRS, and MODIS using fiducial reference measurements
  along the Atlantic Meridional Transect." *Remote Sensing of Environment* 299,
  113844. DOI [10.1016/j.rse.2023.113844](https://doi.org/10.1016/j.rse.2023.113844)
- Brotas, V., Ferreira, A., Veloso, V., et al. (2023). "Assessing phytoplankton
  community composition in the Atlantic Ocean from in situ and satellite
  observations." *Frontiers in Marine Science* 10, 1229692.
  DOI [10.3389/fmars.2023.1229692](https://doi.org/10.3389/fmars.2023.1229692)
  — validates satellite size classes against AMT flow-cytometric cell counts;
  the closest published analogue to what we want to do with MOANA.
- Lange, P. K., Brewin, R. J. W., Dall'Olmo, G., Tarran, G. A., et al. (2018).
  "Scratching Beneath the Surface: A Model to Predict the Vertical Distribution
  of *Prochlorococcus* Using Remote Sensing." *Remote Sensing* 10, 847.
  DOI [10.3390/rs10060847](https://doi.org/10.3390/rs10060847)
  — the same author's precursor; relevant because MOANA is surface-only
  (top 10 m) by construction.

## Coding

Here are guidelines for coding:

- Use Python
- Add inline comments to explain the effort
- Reuse existing code when possible
- Use methods, not classes
- Place MOANA relaated code in the `ioptics/moana/` folder.
- Use matplotlib or seaborn for plotting
- Place import statements at the top of the file.
- Include a description of inputs/outputs in the doc string of all methods

## Prompts

1. **Context**.  Please read all of the files listed in the Context section above.  Then search the NASA PACE websites for any discussion of the MOANA algorithm which was developed for the PACE mission from the Lange+2020 paper. Also search for literature that discusses the AMT24 and its data. If you find any, add them to the Context section above; include DOIs. If you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

2. **Existing code**.  It is possible that the MOANA code is available on the NASA PACE website.  Please search for it and report back in the Reports section below.  You may need to download the (very) large C source code that NASA Ocean Color provides.  If you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

3. **Answers**.  I have answered the questions in the Q&A section below.  Please review them and proceed accordingly.  Use Fable if you can. Log your work.

4. **MOANA Report**.  I have answered the latest round of Q&A; read them.  Now that you have all of these papers and code, generate a report -- `reports/MOANA_Claude_Report.md` with your complete understanding of the MOANA algorithm.  You should include any of the inconsistencies you have identified in your exploration.  Be sure to explain how it works.  And, explain how you'd improve it. Use Fable if you can. Log your work.  If you have any additional questions, please write them in the Q&A section below.

5. **MOANA Report 2**.  I have answered the latest round of Q&A; read them and update the report accordingly.  Use Fable if you can. Log your work.  If you have any additional questions, please write them in the Q&A section below.

6. **Data**.  We are going to reproduce the MOANA algorithm using the AMT24 data.  Please download the data for the AMT24 cruise and put it in the folder `$OS_COLOR/AMT24`.  If you can't find it, provide a list of the places you searched in the Reports section below.  If you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

7. **AMT data**.  You reported back that you could not find the radiometry data for the AMT24 cruise.  How carefully did you search this webpage: `https://amt-uk.org/data/`?  Please try again there and Report back.  Use Fable if you can.  Log your work

8. **AMT24 Radiometry**.  I have finally obtained the radiometry data for the AMT24 cruise.  It is located in `$OS_COLOR/AMT24/Wadiometry`.  Please inspect it and then ask me any questions you have in the Q&A/Radiometry section below.  Log your work.  Use Fable if you can.

9. **More AMT24 Radiometry**.  Please review my answers to the Q&A/Radiometry section below and then ask me any additional questions you have, if any.  Also, revise the prompts that follow (10-13) to reflect the new information as needed. Log your work.  Use Fable if you can.

*(Prompts 10–13 revised 2026-08-15 by prompt 9, to fold in the Q&A #27–32 decisions
and to align file paths with the locked prompt-3 layout — code in `ioptics/moana/`,
tests in `ioptics/tests/` where pytest already discovers them. Original wording:
paths were `code/moana.py`, `tests/moana_tests.py`, `validation/moana_validation.py`.)*

10. **Design**.  Let us generate a design doc for the MOANA algorithm.  We will write it in the file `docs/design/moana_design.md`.  This doc will describe (a) the algorithm itself — retrieval and training, including the operational-vs-ATBD PC-mapping discrepancy and the flag that selects between them; (b) the AMT24 Level-2 → training-matrix pipeline per Q&A #27–34: read the `.sav` files (the ES CSVs are broken), re-derive Rrs from LT/LI/ES with the Lange chain (750–800 nm glint minimisation, tilt/azimuth screening, second-derivative noise filter), resample to 414–660 nm @ 2 nm, median-bin ±15 min onto the flow-cytometry samples; and (c) our approach to coding it in `ioptics/moana/` (methods not classes; raw floats + QC flags rather than NASA's silent clipping).  We will use the files in the Context section above and `reports/MOANA_Claude_Report.md` to help us.  Before proceeding, if you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

11. **More design**.  The design doc looks great.  I have answered your questions 33-35.  Please review them and update the design doc accordingly.  Use Fable if you can. Log your work.  Also create the PML follow-up list as the file `requests/PML_follow_up.md`, as requested in question 35.

12. **Code**.  Let us generate the code for the MOANA algorithm.  We will write it in the `ioptics/moana/` package, split by concern (e.g. `io.py` for the `.sav`/LUT readers, `pipeline.py` for the Level-2 → Rrs chain, `algorithm.py` for standardise → PC scores → regressions, `train.py` for the PCA retraining), reading the LUTs from `ioptics/data/moana/`.  This code will implement the algorithm described in the Design doc.  Before proceeding, if you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

13. **Test**.  Let us generate the tests for the MOANA algorithm.  We will write them in `ioptics/tests/test_moana.py`, following the existing suite's conventions; tests needing the AMT24 tree skip when `$OS_COLOR` is absent, and any Earthdata-dependent test skips when `~/.netrc` is absent (Q&A #16 pattern).  This code will test the algorithm described in the Design doc.  Before proceeding, if you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

14. **Validation**.  Let us generate the validation for the MOANA algorithm.  We will write it in the file `ioptics/moana/validation.py`.  This code will validate the algorithm described in the Design doc, against the three locked targets: (i) reproduce Lange+2020 Tables 1–2 on AMT24 (now data-complete); (ii) apply the model to genuinely held-out cruises — the Brewin et al. 2023 in-situ hyperspectral Rrs for AMT23/25/28 (already downloaded) matched to their flow-cytometry DOIs; (iii) the operational PACE product versus SeaBASS, including the bit-exactness check of one PACE granule (which also settles the Q&A #9 PC-mapping question empirically).  Metrics: log-space bias/MAE/R² with clipped retrievals treated as censored plus a headline unphysical fraction (Q&A #13).  Before proceeding, if you have any questions, please write them in the Q&A section below.  Log your work.  Use Fable if you can.

15. **Validation 2**.  I have answered question 36 and have downloaded the AMT23, AMT25 and AMT28 flow-cytometry data.  Please review the data and then continue with the validation.  Use Fable if you can. Log your work.  Hold off on the Report for now.

16. **Report**.  Let us now update the MOANA report to reflect all of the changes since the last time we touched it.  It is the file `reports/MOANA_Claude_Report.md`.  Have it include a section on "Open items", i.e. the work we still wish to do.  Use Fable if you can. Log your work.

17. **Report additions**.  The report is excellent.  Please also:
    - Our report emphasies the small-scale behaviour (i.e. patchiness) of Synechococcus.  Add to section 9 a discusison of how to reconcile it with the 1km sized pixels of the PACE product.
    - Did we try to validate the product with the other AMT cruises?  If so, what were the results?  If not, can we?
    - The report does not seem to emphasize that the training data are all from 2014 but is being applied now in 2024 and on.  Should we be concerned about this?   For example, how much has the Atlatnic Ocean warmed in these 10 years?  What other changes might we worry about?  Be sure this is discussed in Section 9.
    - Are there other parts of the world where we could train a MOANA like product?  How about the ship-of-opportunity data from Koji Suzuki's program?  Discuss this too

Use Fable if you can. Log your work.

18. **Figures**.  Please generate two new figures for the report:

    - A figure comparing the MOANA values using the wrong coefficients (those in operation today) versus the correct coefficients.   Add this to Section 7.1 and use the one PACE granule
    - A figure showing the analysis in Section 12.3.  The table is nice, but let's include a figure too.

It is possible you will need me to run this on my workstation.  Tell me if so.  Use Fable if you can. Log your work.

19. **PML follow-up**.  I believe that PML has provided the BODC deposit for AMT24.  Please examine those data in the `$OS_COLOR/AMT/AMT24/` folder.  If they are what you need, then proceed with that deferred anlysis and then update the Report and add a comment to the PML follow-up list.  Use Fable if you can. Log your work.

20. **Tom response**.  Tom has written: "the cell counts are the final set of columns in AMT24_JR20140922_AFC_Dataset.xlsx, and these `number codes' map onto the meta data definitions in AMT24_AFC_BODC_sample_metadata_report.html.
e.g. P700A90Z is Abundance of Synechococcus".  Please review this and see what you have for the BODC deposit now.  Use Fable if you can. Log your work.

## Q&A

### Radiometry

*(From prompt 8, 2026-08-15 — Claude. Full inspection findings are in the prompt-8
report below; these questions assume you've skimmed it. Headline: the data are
usable and validation target (i) is unblocked in substance, but what we received is
the upstream calibrated Level-2 stream, not Lange's processed training matrix, and
the CSV export has one real bug.)*

27. **Provenance and documentation.** Who supplied the `Radiometry/level2/` tree —
    PML (Dall'Olmo/Smyth), Lange directly, or BODC — and did it come with any
    README, processing document, or email describing the processing chain? The
    directory holds only the data files. This matters because everything below
    (glint treatment, screening, resampling) is something I had to reverse-engineer
    empirically, and a one-page processing note from the provider would confirm or
    correct it. Also: was `level2` the *only* level they sent, or do higher-level
    (QC-screened / binned / Lange-processed) products exist that we should ask for?
>A. PML (Brewin+).  I will ask them later what the processing chain is.  Do your best for now.

28. **The ES export bug.** `*_ES.dat` is a byte-identical copy of `*_LT.dat` on all
    37 days — the true downwelling irradiance exists only inside the IDL `.sav`
    files (dataset `matrix_es`, which is correct: ~90–175 µW cm⁻² nm⁻¹ with a
    proper solar spectrum shape). `scipy.io.readsav` reads the `.sav` files
    cleanly, so I propose we simply treat the `.sav` as the source of truth and
    ignore the CSVs entirely. OK? And should the bug report — plus the two missing
    days inside the Lange window (see #30) — go on your contact list back to the
    provider?
>A. Yes, I agree.  I will report the bug to PML.

29. **Re-derive Rrs, or use theirs?** The shipped `RRS` is exactly
    `(LT − 0.0280·LI)/ES` — the fixed Mobley ρ, applied at ~0.86 s cadence with
    **no** NIR-residual subtraction (median Rrs(750–796) = 0.0020 sr⁻¹, clearly
    glint-contaminated), **no** azimuth/tilt screening (relative azimuth spans
    −165° to +64°), and **no** binning; ~1.5 % of mid-cruise spectra go negative
    somewhere in 414–660 nm. Lange+2020 instead used their own glint minimisation
    over 750–800 nm, tilt and azimuth screening, and a second-derivative noise
    filter. Since we have LT, LI and ES separately, I propose we **implement the
    Lange chain ourselves from the component streams** and keep the shipped RRS
    only as a cross-check. Confirm? (The alternative — using their fixed-ρ RRS
    as-is — would bake known glint into the PCA training matrix.)
>A. Yes, I agree.  Please implement the Lange chain from the component streams.

30. **Two days are missing inside the training window.** The tree covers DOY
    266–305 but DOY **267, 273 and 298** are absent. Lange's window is
    30 Sep – 1 Nov = DOY 273–305, so **273 (30 Sep) and 298 (25 Oct)** are real
    holes — any flow-cytometry samples on those days lose their Rrs matchup. Do
    you know whether the provider has those days (worth adding to the follow-up),
    or were they lost at sea (instrument down)?
>A. I will ask; good catch.

31. **The spectra are already resampled.** The grid is exactly uniform
    306–796 nm at 3.5 nm (141 bands) — HyperOCR native pixels are ~3.3 nm and
    non-uniform, so the provider has already interpolated once. Getting to MOANA's
    414–660 @ 2 nm therefore means interpolating twice (native → 3.5 → 2), a small
    extra smoothing relative to Lange's native → 2 nm. I'd accept this and move on
    (the effect will be far below the glint/screening choices in #29), but if you
    want fidelity, native-grid files are another item for the provider ask. Which?
>A. Good point.  We will accept this hit for now but will work with PML later to improve.

32. **Binning for the flow-cytometry matchup.** At ~0.86 s cadence there are ~1 M
    spectra; the AMT24 flow-cytometry samples are discrete (CTD stations +
    ~30-minute underway). I propose median-averaging screened spectra within a
    fixed time window around each cell-count sample (±15 min feels right, and the
    median suppresses residual glint spikes), with the spread kept as a per-band
    uncertainty. Happy with that, or do you want a different window/estimator —
    or whatever Lange's text can be pinned down to?
>A. Yes, let's try that.

*(From prompt 9, 2026-08-15 — Claude. All six answers above are reviewed and now
locked; see the prompt-9 log entry for the decision table. Two follow-up questions
arose while turning those decisions into a concrete pipeline — both are choices
that define the training matrix itself, so I'd rather have your call than bury an
assumption in code.)*

33. **Depth cut for the flow-cytometry matchup.** The AMT24 FCM dataset (Tarran &
    Zubkov) spans 200 m → surface, but MOANA is surface-only by construction and
    the ATBD frames the products as "top 10 m". Lange+2020 trained on surface
    samples (the ship's underway supply draws from ~5–7 m, plus surface CTD
    bottles). I propose: keep FCM samples with **depth ≤ 10 m**, and where a
    station has several samples in that range, take the shallowest (not the mean,
    to avoid mixing across the near-surface gradient Lange 2018 documents for
    *Prochlorococcus*). Confirm, or would you rather average 0–10 m?
>A. I confirm.  Take the shallowest sample.

34. **Which SST feeds the Prochlorococcus regression for the AMT24 work?** Three
    candidates on disk: (a) the HSAS ancillary `SST` column (co-sampled with every
    spectrum, but garbage early in the cruise — 0 to 12112 °C on DOY 266); (b)
    `uway_sst` from the Jordan et al. 2025 netCDF (1-minute, QC'd, published,
    8.2–29.3 °C); (c) CTD station temperatures. I propose **(b) as primary** —
    interpolated to each Rrs bin's timestamp — with (a) kept as a cross-check
    (mid-cruise they should agree to ~0.1 °C), because (b) is the QC'd, citable
    record of the same intake stream. Confirm?
>A. Yes, use (b) as the primary SST source with (a) as a cross-check.

*(From prompt 10, 2026-08-16 — Claude. One more, found while writing the design
doc's matchup section.)*

35. **The underway flow-cytometry samples are not in the BODC deposit.** The
    `AMT24_JR20140922_AFC_Dataset.csv` we hold is **CTD-bottle-only**: 814
    samples, 68 stations, all `Gear=CTD`. But Lange+2020 trained on n = 73–78
    including the ~30-minute **underway** FCM sampling across the South Atlantic
    front — and their own sensitivity test shows that dropping the underway
    samples (CTD-only retraining) degrades *Synechococcus* from MAE 1.27 to 1.37.
    So until we obtain the underway FCM, our Tables 1–2 reproduction runs in the
    known-degraded CTD-only configuration (the design doc §5.6 quantifies this
    and pins the fair comparison line). Could you add the underway FCM to the
    PML/Tarran follow-up list? (It may simply be a second BODC deposit we
    haven't found — I will also re-search BODC for an AMT24 underway FCM series
    before prompt 11 if you prefer.)
>A. Yes, I will add the underway FCM to the PML/Tarran follow-up list.  I've somewhat lost track of that PML list.  can you generate it for me as the file `requests/PML_follow_up.md`?

*(From prompt 14, 2026-08-16 — Claude. Validation ran; two data items are all
that separates us from finishing targets (ii) and (iii)b.)*

36. **Two browser downloads to score the remaining targets.** Target (ii)'s
    *retrievals* are done (90 stations on AMT23/25/28 from the Brewin Rrs), but
    *scoring* them needs those cruises' flow-cytometry counts, which are not on
    disk — the same BODC browser-only situation as Q&A #24. Could you grab, into
    `$OS_COLOR/AMT24/` (or a sibling folder — the loader takes a path):
    - AMT23: https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/a2104adc-e990-6789-e053-6c86abc0d557/
    - AMT25: https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/a2104adc-e98e-6789-e053-6c86abc0d557/
    - AMT28: https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/a147c314-688b-55e9-e053-6c86abc0dc81/

    Separately, target (iii)b (PACE vs in-situ counts) has its matchup machinery
    ready but needs the SeaBASS picophytoplankton cell counts, which require a
    SeaBASS account (Q&A #24 left this open). Do you have one, or should that
    stay parked?

>A. I have downloaded those datasets and put each in its own folder under a new folder `$OS_COLOR/AMT/`.  And, yes, I do have a SeaBASS account.  The username is `profx`.  Let me know if you need the password and I'll add it to my .bashrc file for you to read (but not broadcast)

### From prompt 1 (2026-08-01) — Claude

1. **The training PCA loadings are the one hard blocker.** Neither Lange+2020 nor
   the ATBD publishes the loading matrix `V(λ, i)`. Operationally it lives in
   `$OCDATAROOT/common/pca_picophyto.h5` and the regression coefficients in
   `$OCDATAROOT/common/picophyt.json` — two OCSSW *data* files, not source.
   Without `V` we can apply MOANA only in "published-coefficients" mode and only
   if we can get that file. Which route do you want?
   (a) install OCSSW/SeaDAS and pull the two files from `$OCDATAROOT/common`
   (this is prompt 2's territory, so I did not start it);
   (b) email Lange / Cetinić / Minwei Zhang for `pca_picophyto.h5`;
   (c) re-derive the PCA ourselves from AMT24 — scientifically the most
   satisfying, and it makes "reproduce the model" literal, but it needs the
   AMT24 hyperspectral Rrs (see #2) and will not reproduce NASA's coefficients
   exactly (sign/ordering of eigenvectors is arbitrary).
>A. We will install OCSSW/SeaDAS and pull the two files from `$OCDATAROOT/common`.

2. **Where is the AMT24 hyperspectral Rrs?** The *cell counts* have a clean BODC
   DOI. I could **not** confirm a public archive for the AMT24 HyperSAS Rrs
   (414–660 nm at 2 nm) that trained the PCA. Jordan+2025 puts AMT24 *IOPs and
   pigments* on SeaBASS/Zenodo but that is `ap`/`cp`/HPLC, **not** Rrs. Prompt 3
   says "download the data from the AMT24 website" — is it fine if I spend the
   effort digging through SeaBASS `PML/AMT` and BODC for the radiometry, or do you
   already know where the Rrs lives / would you rather just ask PML?
>A. Yes, please spend the effort digging through SeaBASS `PML/AMT` and BODC for the radiometry.

3. **Scope — how does MOANA sit inside IOPtics?** This matters before prompt 4
   (design doc). MOANA is *not* an IOP algorithm: it outputs cell abundances, has
   no `a(λ)`/`bb(λ)`, no BING forward model, no χ² Rrs closure, no free parameters
   to count, so ~none of the existing `ioptics.metrics` battery (§2 closure, §3
   AIC/BIC/ΔBIC, §4 coverage) applies. Its accuracy metrics (log-space bias, MAE,
   R²) *do* already match ours. Do you want it as
   (a) a separate "empirical PCC products" track alongside the IOP track,
   (b) a genuine `AlgorithmSpec` variant that bypasses the BING wrapper, or
   (c) purely exploratory under `papers/moana/` with no `ioptics/` changes yet?
   Prompts 5–7 point at bare `code/`, `tests/`, `validation/` directories that
   don't exist in this repo and sit outside the `ioptics/` package layout — I'd
   like to confirm that's deliberate rather than shorthand.
>A. It is separate from the IOP track, and we will not be evaluating it against the IOP track.  We will be exploring it separatlely and eventually will write a report with our main findings.

4. **What is the primary validation target?** Ranked candidates:
   (i) reproduce Lange+2020 Tables 1–2 on AMT24 (needs #1 and #2);
   (ii) apply the published coefficients to AMT20/22/23/25/28 — genuinely
   held-out cruises, DOIs now in Context, which is exactly what prompt 7 asks;
   (iii) validate the *operational* PACE OCI MOANA L3/L4 product against SeaBASS
   cell counts (note the ATBD admits this has never been done — a real, publishable
   gap);
   (iv) run MOANA on L23 synthetic Rrs. I'd argue (ii) is the best first target
   and (iv) is nearly meaningless — see #5.
>A. We will do (i), (ii) and (iii) as best we can.

5. **SST, and why L23 may be a dead end for MOANA.** The *Prochlorococcus*
   equation needs `log10(SST)`, and operationally that SST is the PACE reference
   field (GHRSST CMC L4). L23 is synthetic and has **no SST at all**, so the Pro
   model cannot be run on L23 without inventing a temperature — and *Synechococcus*
   / picoeukaryotes would be extrapolating a PCA trained on real Atlantic spectra
   onto Hydrolight spectra. Do you agree we skip L23 for MOANA? And for real
   cruises, do we have GHRSST/CMC under `$OS_COLOR`, or should I add an `ocpy`
   loader?
>A. Yes, we will skip L23 for MOANA.

6. **A discrepancy to resolve before coding.** The ATBD writes a *different,
   non-contiguous* PC subset per taxon — Pro: `U1,U2,U6,U17`; Syn:
   `U1,U2,U3,U5,U8,U9,U10,U11,U12,U13`; picoeuk: `U2,U3,U4,U5,U10,U13,U15` — but
   `get_Cpicophyt.c` multiplies coefficient `i` by PC `i` **contiguously** over the
   first `npc` PCs. The only way both are true is if `picophyt.json` zero-pads the
   unused slots. Relatedly, the paper says 14 PCs were retained yet the indices run
   to 17. Do we (a) treat the ATBD table as authoritative and zero-pad, or (b) hold
   until we can actually read `picophyt.json`? I lean (b), since it also
   independently checks the coefficient values.
>A. Let's do (b)

7. **Which wavelength grid?** The ATBD is internally inconsistent: its abstract
   says the algorithm needs valid Rrs over **395–705 nm**, while its Mathematical
   Theory section (and Lange+2020) says **414–660 nm at 2 nm** = 124 bands. I
   propose we implement 414–660 @ 2 nm because that is what the PCA was trained on,
   and treat 395–705 as a mere data-availability screen. Confirm?
>A. We will implement 414–660 @ 2 nm.

8. Minor but worth deciding early: the OCSSW code standardizes each spectrum with
   the **sample** standard deviation (N−1 denominator). Lange+2020 Eq. 3 just says
   "standard deviation". Since the scores scale inversely with this, N vs N−1 shifts
   every retrieval slightly. I'll match OCSSW (N−1) unless you object.
>A. We will match OCSSW (N−1).

### From prompt 2 (2026-08-01) — Claude

*(Questions 1, 6, 7 and 8 above are now settled as fact by the OCSSW data files —
see the prompt-2 report. Q&A #1's blocker is gone entirely: we have the loadings.)*

9. **The one that matters.** The operational `picophyt.json` puts the last
   Prochlorococcus coefficient on **PC7** (ATBD says U17) and the last
   *Synechococcus* coefficient on **PC16** (ATBD says U13); all values are
   identical and picoeukaryotes agree exactly. The ATBD's assignment uses **14
   distinct PCs** across the three taxa — matching Lange+2020's "14 PCs survived
   the significance cut" — while the JSON uses **15**. That is decent evidence the
   ATBD is faithful and the operational LUT has two coefficients in the wrong
   slots, i.e. the shipping PACE products may not be the published algorithm.
   I propose we (a) implement both mappings behind a flag, (b) settle it
   empirically by reproducing a real PACE MOANA granule, and (c) email
   Lange/Cetinić/Zhang once we have that evidence in hand. Do you want to be the
   one to make that contact, and should we hold off publishing the discrepancy
   until they reply?
>A. Yes, let's do (a) and (b) and I will make the contact later.

10. **Bit-exactness check before AMT24?** I'd like to download one PACE OCI L3M Rrs
    granule plus its matching L3M MOANA granule and verify our reimplementation
    reproduces NASA's three fields pixel-for-pixel. It is the cleanest possible
    test of our transcription *and* it settles #9. It needs an Earthdata login for
    `earthaccess` — do you have credentials set up (`~/.netrc` or
    `EARTHDATA_USERNAME`/`_PASSWORD`), or should I plan around that?
>A. I have credentials set up (`~/.netrc`) 

11. **How should clipped zeros be treated?** OCSSW silently sets negative
    Prochlorococcus to 0 rather than flagging it, and casts all three products to
    truncated `int32`. Zeros are fatal to log-space bias/MAE. Do we (a) keep NASA's
    behaviour exactly and mask zeros at the metrics stage, or (b) retain the raw
    float value and record the clip as a QC flag, so we can report how often MOANA
    goes unphysical? I lean (b) for our own runs, with (a) available for
    bit-comparison against NASA.
>A. Let's do (b)

12. **Where do the two LUT files live for us?** They are currently at
    `$OS_COLOR/ocssw/share/common/{pca_picophyto.h5,picophyt.json}` — the canonical
    OCSSW location. Since they are 24 KB and 390 bytes and are the crux of the whole
    exercise, do you want copies committed into the repo (e.g. `ioptics/moana/data/`)
    for reproducibility, or should the code read them from `$OCDATAROOT`/`$OS_COLOR`
    and keep the repo data-free as the design doc prefers?
>A. Yes, copy into the repo (`ioptics/data/moana/`) for reproducibility.

### From prompt 3 (2026-08-01) — Claude

*(All 12 answers above are acted on; see the prompt-3 report. The unanswered half of
#5 — where GHRSST SST comes from — is now resolved: `CMC0.1deg-CMC-L4-GLOB-v3.0` is
reachable via `earthaccess` at 0.1°, matching the MOANA grid. No `ocpy` loader
needed.)*

13. **The 18 % clipping problem may be the most interesting result we have.**
    *Prochlorococcus* is exactly 0 in 18.0 % of real ocean retrievals in the one
    day I examined, because the regression goes negative and OCSSW silently clamps
    it; *Synechococcus* is 0 in another 8.9 % from integer truncation. For
    validation targets (i)–(iii), how should those be treated? Options: (a) exclude
    them (cleanest, but flatters the algorithm by hiding its failures); (b) treat
    them as **censored data / upper limits**, which is defensible and which the BING
    paper already has machinery for; (c) keep them as zeros and report a separate
    "fraction unphysical" statistic alongside bias/MAE. I lean (b) **and** (c)
    together — report the failure rate as a headline number, since a fifth of the
    field going unphysical is arguably more important than the MAE of the
    remainder.
>A. (b) and (c) together.

14. **How hard do we push the product-level defects?** The land-as-254 encoding,
    the `INT32_MIN` pixels, and the unenforced `valid_max` (picoeukaryotes 68× over)
    are defects in a shipping NASA product, not in the science. They belong in the
    final report, but do you want them (a) as a prominent section — they are
    genuinely useful to the community and to OBPG — or (b) as an appendix, keeping
    the report focused on algorithm accuracy? Also: should I report the land-254
    trap to OBPG, or does that go with your Q&A #9 contact?
>A. (b).  Remind me to include this when I reach out to them.

15. **Matchup design for target (iii).** Since there is no L3M MOANA, the
    operational validation has to use L4M composites (daily / 8-day / monthly).
    Bailey & Werdell (2006) — the protocol the ATBD names — is designed around L2
    granules with a 5×5 box and strict time windows, which does not map cleanly onto
    a daily 0.1° composite. Do we (a) use daily L4M and accept the looser matchup,
    (b) run `l2gen`/`l3gen` ourselves from L2 AOP to make true L2 MOANA and follow
    Bailey & Werdell properly, or (c) simply document the mismatch and use daily
    L4M? (b) is the rigorous option and OCSSW is now installed, but it is
    substantially more work.
>A. (c), and we will consider (b) later.

16. Minor: `earthaccess` is now a real dependency of the MOANA track. Add it to
    `requirements.txt`, or keep the MOANA track's dependencies separate since it is
    a side track (answer #3) and CI shouldn't need Earthdata credentials? I'd add it
    to `requirements.txt` but guard the tests with a skip when `~/.netrc` is absent,
    mirroring the existing `$OS_COLOR` skip pattern.
>A. Add it to `requirements.txt`, and guard the tests with a skip when `~/.netrc` is absent, mirroring the existing `$OS_COLOR` skip pattern.

### From prompt 4 (2026-08-01) — Claude

*(#16 is done: `earthaccess` added to `requirements.txt` with a note about the
`~/.netrc` test guard. #13, #14, #15 are reflected in the report — clipped values as
censored data plus a headline unphysical fraction; product defects in Appendix A
rather than the main argument; L4M matchups with the Bailey & Werdell departure
documented.)*

17. **Is §9 a reading list or a work plan?** The report proposes eight improvements.
    Three are cheap and would make a real contribution if we actually implemented and
    benchmarked them against stock MOANA — that would turn the report from "here is
    what MOANA is" into "here is how much better it could be", which is a much
    stronger paper. Specifically: **§9.3's out-of-domain flag** (the truncated-basis
    reconstruction residual) needs **no retraining and no new data** — we could build
    and demonstrate it on PACE granules immediately. **§9.1** (log-link GLM to kill
    the 18 % negatives) and **§9.2** (feed back the two amplitude numbers that
    standardisation discards) both require *retraining*, so they are blocked on the
    same AMT24 training data as validation target (i). Do you want me to scope any of
    these in, or is the report the deliverable and we stay descriptive?
>A. The report is the first deliverable, so we should stay descriptive.  We will then code up MOANA in Python.

18. **Should we ask NASA for the PCA eigenvalues?** They are published nowhere —
    not in the paper, not the ATBD, not even in the LUT, which ships only the
    loadings. Without them we cannot independently confirm the ">96 % on PC1" figure,
    cannot whiten or variance-weight the scores, and cannot reproduce the original
    component-selection decisions. Worth adding to your contact list alongside the
    §7.1 discrepancy?
>A. Yes, add the list of things to ask/tell NASA.

19. **Where should the report live?** It is at `reports/MOANA_Claude_Report.md` as
    instructed, with its figure in `reports/figures/` and its generator in
    `reports/scripts/` — matching the `gloria_fits_report.md` precedent. Do you also
    want it rendered into the Sphinx/readthedocs tree
    (`docs/source/reports/moana/`) like the IOP sweep reports, or does the MOANA
    track stay outside the published site until it is further along?
>A. It stays outside for now.

20. **How many more figures do you want?** The report currently has one (the PC
    loadings basis). Obvious further candidates, all cheap: a map of the PACE MOANA
    product with the land-254 and clipped-zero masks made visible; a histogram of the
    clipped fraction by taxon; and — once we run the bit-exactness test — a
    scatter of our retrieval against NASA's. Say the word and I will add them via
    `reports/scripts/moana_report_figs.py`, which is already structured for it.
>A. I like you suggestions, please add all of them.

### From prompt 5 (2026-08-01) — Claude

*(All four answers are acted on. #17 → §9 stays descriptive and says so explicitly;
#18 → new report §11 consolidates everything to ask/tell NASA, and the §10 reminder
now points at it; #19 → nothing added to the Sphinx tree; #20 → two of the three
figures are in, the third is deferred because it needs the implementation.)*

21. **The clipping turns out to be diagnosable, which changes its status.** Making
    the mask figure exposed something I had not seen from summary statistics: the
    18 % clipped *Prochlorococcus* pixels are strongly spatially organised — median
    |latitude| 47.9° vs 31.7°, 74.6 % poleward of 40°, and 4.5× the picoeukaryote
    abundance where the model fails versus where it succeeds. So the failures are
    **out-of-domain extrapolation into cold productive water**, not noise. That makes
    the §9.3 out-of-domain flag testable *without any retraining*: compute the
    truncated-basis reconstruction residual per pixel and check whether it predicts
    the clipping. Given #17 says stay descriptive, I have only written this up — but
    it is a one-afternoon experiment that would likely be a genuine contribution, and
    it needs no data we do not already have. Worth slotting in after the Python
    implementation?

22. **Should the report get a "how to read the product" section?** Between the
    land-254 trap, the unflagged clipped zeros, the unenforced `valid_max` and the
    `INT32_MIN` pixels, there is now enough for a short, practical "if you are using
    PACE MOANA, mask these four things" recipe. That would be the most immediately
    useful thing in the document for anyone outside this project, and it is
    orthogonal to the accuracy analysis. Add it as a §12, or keep the report
    analytical and save it for a community-facing note later?

23. Minor, and only if you care: the figures are light-mode PNGs, matching the
    existing `gloria_fits_report.md` figures. The house data-viz guidance wants a
    dark-mode variant for anything rendered as a web page. Since #19 keeps this out
    of the Sphinx site, I have not made dark variants. Flag it if the report is ever
    published so I can add them.

### From prompt 6 (2026-08-01) — Claude

24. **Two BODC downloads need your browser.** The Published Data Library is fully
    JavaScript-rendered and exposes no machine-reachable download route (I tried
    `data/`, `download/`, `files/`, `?download=1`, and checked DataCite for a
    `contentUrl` — nothing). Could you grab these two into `$OS_COLOR/AMT24/`?
    - **flow-cytometry cell counts** (the MOANA training truth, and the thing
      validation target (i) cannot proceed without):
      https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/a2104adc-e98f-6789-e053-6c86abc0d557/
    - **CTD profiles** (T, S, fluorescence, PAR), useful for station-matched SST and
      for a Chl sanity check:
      https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/c0152486-cdae-4744-e053-6c86abc050c2/

    The same applies to the five held-out cruises' cell counts (DOIs in Context) when
    we get to validation target (ii). Alternatively, if you have a SeaBASS account,
    tell me and I will structure the code to read `.sb` files as well.

    **Added after prompt 7 — this one is now the highest priority of the three:**
    - **in-situ hyperspectral `Rrs` + Secchi/Forel-Ule/Chl for AMT23, 25, 26, 28**
      (Brewin et al. 2023), which unblocks validation target (ii):
      https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6/

>A. I have manually downloaded and unzipped all of these

25. **Do you want me to draft the PML data request?** The training `Rrs` is not
    published anywhere (see the prompt-6 report for the three independent lines of
    evidence), so it has to be asked for. The cruise report names the team — Giorgio
    Dall'Olmo, Jelizaveta Ross, Rafael Rasse Boada, Tim Smyth at PML — and Priscila
    Lange was aboard, so there are two routes. I would ask specifically for **the
    processed 414–660 nm at 2 nm matrix that was fed to `prcomp`**, not just "the
    AMT24 radiometry": Lange et al. applied their own glint minimisation over
    750–800 nm, tilt/azimuth screening and a second-derivative noise filter, so raw
    or differently-processed HyperSAS data will not reproduce their PCA. I can draft
    the email if useful, or fold it into the §11 list you already have.

>A. I have made the request.

26. **Should the AMT24 IOPs become an IOPtics dataset in their own right?** The
    Zenodo file is a serious dataset — 176-wavelength particulate `ap`/`bp`/`cp`
    *with per-point uncertainties*, at 1-minute resolution, plus HPLC. That is real
    in-situ spectral IOP truth, which is exactly what the IOP track wants and which
    PANGAEA/GLORIA only partly provide. It is unrelated to MOANA (no `Rrs`, so no
    retrieval), but it would be a shame to leave it sitting in `$OS_COLOR/AMT24/`
    unused. Worth a `datasets.py` adapter later, or out of scope?

>A. Yes, let us make it a IOPtics dataset in its own right.

## Reports

### Prompt 8 (2026-08-15) — AMT24 HyperSAS radiometry inspected: the training data exist and are usable, but they are *upstream* Level-2 with one export bug

*(The prompt says `$OS_COLOR/AMT24/Wadiometry`; the folder on disk is
`$OS_COLOR/AMT24/Radiometry` — assumed typo.)*

**Headline: validation target (i) is unblocked in substance.** We now hold the
HyperSAS radiometry from the very cruise that trained MOANA — nearly a million
spectra covering the full transect. But what we received is the *calibrated,
un-screened* Level-2 stream, not the processed 414–660 nm matrix Lange et al. fed
to `prcomp`, so reproducing the training set means re-implementing their QC chain
(details in §4). And one file per day is broken (§3).

#### 1. What is on disk

`$OS_COLOR/AMT24/Radiometry/level2/` — 8.7 GB, 37 day-of-year directories
(`2014266` … `2014305`, i.e. 23 Sep – 1 Nov 2014; DOY 267, 273 and 298 absent).
Each day holds five files, e.g.:

| File | Content |
|---|---|
| `AMT24_HSAS_2014-280_ES.dat` | CSV — **broken: byte-identical copy of the LT file** (all 37 days) |
| `AMT24_HSAS_2014-280_LI.dat` | CSV — sky radiance |
| `AMT24_HSAS_2014-280_LT.dat` | CSV — total upwelling radiance |
| `AMT24_HSAS_2014-280_RRS.dat` | CSV — remote-sensing reflectance |
| `AMT24_HSAS_2014-280.sav` | IDL save — `matrix_es/lt/li/rrs`, the **source of truth** (its `matrix_es` is correct) |

Totals: **993,417 spectra** at ~0.86 s cadence (daylight hours only), tracking the
ship from 52.2°N (Humber approaches) to 50.5°S. Row counts per day range 14k–36k.

#### 2. File layout

Every matrix is 158 columns: **14 ancillary** (ship start time, time [decimal
hour], lat, lon, atmospheric pressure, PAR, wind direction/speed, ship
orientation, salinity, SST, delta-azimuth, azimuth, zenith angle) + **141
wavelengths, 306.0–796.0 nm at exactly 3.5 nm** + 3 trailing ancillary (THS
compass, pitch, roll). The CSV header lists the trailing three *after* the
wavelengths — easy to mis-parse if you assume all ancillary columns lead.
`scipy.io.readsav` reads the `.sav` files cleanly and their RRS/LT/LI match the
CSVs value-for-value; only ES differs (because the CSV ES is wrong, §3).
MOANA's 414–660 nm window is fully covered. Units are consistent with
µW cm⁻² nm⁻¹ (ES; midday max ~175 at 550 nm, correct solar shape) and
µW cm⁻² nm⁻¹ sr⁻¹ (LT, LI), giving Rrs in sr⁻¹.

#### 3. The one real defect: the ES CSVs

`md5sum` proves `*_ES.dat` ≡ `*_LT.dat` on **all 37 days** — the export script
evidently wrote the LT matrix twice. The true ES survives inside every `.sav`
(`matrix_es` ≠ `matrix_lt`, physically sensible). Consequence: **all code must
read the `.sav` files**; the CSVs are at best a convenience view and at worst a
trap. Worth reporting back to the provider (Q&A #28).

#### 4. Processing level, determined empirically

Solving `ρ = (LT − RRS·ES)/LI` per spectrum per band gives **ρ = 0.0280 exactly,
zero spread within and across spectra, zero wind-speed correlation** — i.e. the
shipped RRS is the plain fixed-ρ Mobley correction, nothing more:

- **No NIR-residual/glint-offset subtraction:** median Rrs(750–796) = 0.0020 sr⁻¹.
- **No geometry screening:** relative azimuth spans −165° to +64°; solar zenith
  up to 85° (dawn/dusk included).
- **No binning:** native ~0.86 s samples.
- **Negative Rrs** somewhere in 414–660 nm in ~1.5 % of mid-cruise spectra.

None of Lange+2020's steps (their 750–800 nm glint minimisation, tilt/azimuth
screening, second-derivative noise filter, 2 nm interpolation) has been applied.
So we can — and must — run that chain ourselves from LT/LI/ES (Q&A #29), which is
actually the *better* position: we control the QC and can propagate it.

#### 5. Smaller findings

- **Grid already resampled once:** exactly-uniform 3.5 nm is not HyperOCR-native
  (~3.3 nm, non-uniform), so reaching 414–660 @ 2 nm means a second interpolation
  (Q&A #31).
- **Two holes inside the Lange window** (30 Sep – 1 Nov = DOY 273–305): DOY 273
  and 298 are missing (Q&A #30). DOY 267 is also absent but predates the window.
- **Ancillary quality:** mid-cruise SST/salinity are sensible (e.g. 26.9 °C /
  36.6 at 19°N, tracking the Jordan et al. `uway_sst` climatology) but DOY 266
  (Humber transit) carries garbage (SST 0–12112 °C, salinity 0.01) — that day is
  pre-window and should simply be dropped. Wind speed, pitch/roll and full solar
  geometry are present, i.e. everything Lange-style screening needs.

#### 6. Consequence for the plan

With the flow-cytometry counts (in hand since prompt 6's Q&A #24) and now the
radiometry, **target (i) — reproducing Lange's PCA and Tables 1–2 — is data-complete**,
pending the processing-chain decisions in Q&A #27–32. The natural next step after
those answers: build the screened, binned, 2 nm 414–660 training matrix and
compare our `prcomp` loadings against the operational `pca_picophyto.h5`.

### Prompt 7 (2026-08-01) — Re-searched amt-uk.org properly. Found published AMT hyperspectral Rrs — for the *other* cruises

**First, how carefully I searched it the first time: not carefully at all.** In
prompt 6 I listed `https://amt-uk.org/data/` among the sources "checked and found
nothing". What I actually did was issue an HTTP status probe that returned 200. I
never read the page. Listing it as searched was wrong, and the challenge was fair.

**Second, the good news: searching it properly paid off** — not on the AMT24 `Rrs`,
but on something that materially changes the plan.

#### 1. What amt-uk.org/data/ actually says

The page is a signpost, not a repository. Its substantive content, verbatim:

> "All data sets are avaliable from the **British Oceanographic Data Center (BODC)
> website**. The remaining AMT data sets are **available upon request to BODC**, with
> future web delivery under development."

Plus a link to an "AMT data policy" page
(`bodc.ac.uk/projects/data_management/uk/amt/data_policy/`) which is **dead — HTTP
404**, as is the parent `…/uk/amt/` directory.

To be thorough this time rather than assume: the site is a JavaScript-rendered
WordPress SPA (plain fetches return "Loading..."), so I also queried its REST API
directly — `wp-json/wp/v2/pages?per_page=100` returns **12 pages total**, and the only
ones matching data/optical/radiometry keywords are the Data page above and the
privacy policy. The AMT24 cruise page (`/cruises/amt24/`) carries a personnel table
and two cruise-report links and **no data links at all**. So amt-uk.org genuinely
hosts no data; the whole site defers to BODC.

Two useful details did come out of the cruise page: **Gavin Tilstone (PML)** was
aboard AMT24 — he is a co-author on the AMT radiometry papers — and **Rob Thomas
(BODC)** sailed as the data manager.

#### 2. Following the BODC pointer — the find I had missed

Because amt-uk.org's answer is "ask BODC", I went back to BODC and searched for AMT
*optical* holdings rather than AMT24 holdings. That surfaced a dataset my prompt-6
searches never touched, because every query I ran had been scoped to "AMT24":

> **Brewin, Pitarch Portero, Dall'Olmo, van der Woerd, Lin, Sun & Tilstone (2023).**
> "Modern and traditional optical measurements, and environmental data, collected on
> four Atlantic Meridional Transect cruises between 2013 and 2018."
> NERC EDS BODC. DOI **10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6**

Contents, per its abstract: *"Secchi depth, Forel-Ule colour, Chlorophyll-a
concentration, **hyperspectral remote-sensing reflectance**, diffuse and beam
attenuation, and other auxiliary data … at a total of 127 stations."*
Formats: Binary + Delimited. Collected 2013-10-10 to 2018-10-27.

**The cruises are AMT23, AMT25, AMT26 and AMT28 — AMT24 is not among them.** I want
to be exact about that, because "four AMT cruises between 2013 and 2018" reads at
first glance as though it should include 2014.

#### 3. Why this matters more than it first appears

Lange et al.'s five held-out validation cruises are AMT20, 22, **23**, **25** and
**28**. This dataset covers **three of those five**, with *in-situ hyperspectral*
`Rrs` — and we already hold the matching flow-cytometry cell-count DOIs for all
three (in Context since prompt 1).

That means **validation target (ii) is no longer blocked**, and it can be done
*better than Lange et al. did it themselves*: their held-out test applied the
multispectral model to Aqua-MODIS reflectance (Pro MAE 2.26, bias +75 %), whereas we
could apply the **hyperspectral** coefficients — the ones MOANA actually ships — to
in-situ hyperspectral `Rrs` on genuinely held-out cruises. No published number exists
for that combination.

#### 4. What is still true from prompt 6

The prompt-6 conclusion stands, correctly scoped: **the AMT24 `Rrs` specifically is
not published anywhere.** The DataCite enumeration (7 AMT24 datasets, none
radiometry) and the SeaBASS enumeration (no radiometry on any AMT cruise) both hold.

But one *implication* I left hanging there was wrong and I want to correct it
explicitly: I framed the AMT radiometry as data that "was always intended for cal/val,
just never deposited". That is true of AMT24, but **not** of AMT generally — BODC has
published AMT hyperspectral `Rrs` for four other cruises. So there is a real
precedent, which makes a BODC request for the AMT24 radiometry considerably more
promising than I implied: it is the same data type, the same programme, the same
people (Dall'Olmo, Tilstone), and the same archive that has already done it four
times.

#### 5. Sources searched in this pass

| Source | Result |
|---|---|
| `https://amt-uk.org/data/` (read in full) | signpost only; "remaining data available on request to BODC" |
| `amt-uk.org` WP REST API, all 12 pages | no data holdings |
| `https://amt-uk.org/cruises/amt24/` (full HTML + links) | personnel + cruise reports; no data links |
| `bodc.ac.uk/projects/data_management/uk/amt/data_policy/` | **404** |
| `bodc.ac.uk/projects/data_management/uk/amt/` | **404** |
| DataCite: AMT + reflectance / radiometry / optical underway | **found 10.5285/f3198e10…** (AMT23/25/26/28) |
| DataCite: Brewin / Tilstone / Dall'Olmo + AMT + reflectance | no further datasets |
| BODC PDL download endpoints for 10.5285/f3198e10… | browser-only, as before |

#### 6. Recommended next actions

1. **Download 10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6** (browser needed — added
   to Q&A #24's list). Then confirm from the file itself which wavelengths the `Rrs`
   spans and whether it covers 414–660 nm at usable resolution.
2. **Retarget validation (ii) onto AMT23/25/28** using in-situ hyperspectral `Rrs`.
3. **Keep the AMT24 request open**, now citing the four-cruise dataset as precedent.

### Prompt 6 (2026-08-01) — AMT24 data: what I obtained, and why the training Rrs is not public

> **Partly superseded — read the prompt-7 report above first.** Two corrections:
> (a) my source list below claims `amt-uk.org/data/` was checked, when in fact I only
> ran an HTTP status probe on it and never read it; (b) §4 below implies AMT
> radiometry was never deposited anywhere — untrue in general, since BODC publishes
> AMT hyperspectral `Rrs` for AMT23/25/26/28 (DOI 10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6).
> The central conclusion, that the **AMT24** `Rrs` specifically is unpublished, stands.

**Headline: the AMT24 hyperspectral `Rrs` that trained MOANA is not publicly
archived anywhere.** I can now say that with evidence rather than as a failure to
find it — see §3. Everything else about AMT24 that *is* public is now in
`$OS_COLOR/AMT24/`.

#### 1. What is now on disk

| File | Size | Source |
|---|---|---|
| `amt24_final_with_debiased_chl.nc` | 571 MB | Zenodo [10.5281/zenodo.12527954](https://doi.org/10.5281/zenodo.12527954), CC-BY-4.0 (Jordan et al. 2025) |
| `AMT24_cruise_report.pdf` | 4.7 MB, 170 pp | PML, DOI [10.17031/06j7-s435](https://doi.org/10.17031/06j7-s435) |

Both downloaded openly, no credentials. I opened the netCDF to verify rather than
trusting the description, and it is **richer than the paper title suggests** —
110 variables on a 1-minute underway grid, 32,174 records, 2014-09-26 to 2014-10-31,
latitude −48.3° to +47.8° (the full transect):

- **IOPs:** AC-S and AC-S2 particulate absorption / scattering / beam attenuation
  (`acs_ap`, `acs_bp`, `acs_cp` and `acs2_*`) on **176 wavelengths**, each with an
  uncertainty estimate; AC-9 equivalents at 9 wavelengths; BB3 backscattering at
  3 wavelengths; C-star `cp`.
- **Pigments:** HPLC at 26 stations.
- **Chl:** AC-S-derived, including the debiased variants.
- **Underway ancillaries:** lat, lon, PAR, total irradiance, air temperature,
  humidity, barometric pressure, salinity — and, importantly, **`uway_sst`
  (8.22–29.27 °C)**.

That last one matters more than it looks. The *Prochlorococcus* model needs
`log₁₀(SST)`, and this gives us **in-situ SST at 1-minute underway resolution along
the whole transect** — better matched to Lange et al.'s underway `Rrs` sampling than
the CTD-station temperatures, and it removes any need to fetch a satellite SST field
for the AMT work. Confirmed absent, as expected: nothing resembling `Rrs`,
water-leaving radiance, or sky radiance (`uway_tir*` and `uway_par` are downwelling
meteorological sensors, not water-leaving signals).

#### 2. What exists but needs a manual step or an account

- **AMT24 flow-cytometry cell counts — the MOANA training truth.** BODC,
  DOI [10.5285/a2104adc-e98f-6789-e053-6c86abc0d557](https://doi.org/10.5285/a2104adc-e98f-6789-e053-6c86abc0d557).
  Format "Delimited". The BODC Published Data Library catalogue page is entirely
  JavaScript-rendered and I found **no** machine-reachable download route: `data/`,
  `download/`, `files/` and `?download=1` all return the SPA shell, and the DataCite
  record carries no `contentUrl`. **This needs a browser — a two-click job for you.**
  It is the single most important file for validation target (i).
- **AMT24 CTD profiles** (T, S, fluorescence, transmittance, downwelling PAR), BODC
  DOI [10.5285/c0152486-cdae-4744-e053-6c86abc050c2](https://doi.org/10.5285/c0152486-cdae-4744-e053-6c86abc050c2).
  Worth having: it supplies the **in-situ SST** at the CTD stations, which is what
  Lange et al. used for the *Prochlorococcus* model (SBE 9/11 on the rosette). Same
  browser-only limitation.
- **SeaBASS `PML/AMT/AMT24/archive/`** holds three files —
  `AMT24_HPLC_20140929_20141031_v20240610.sb` and two
  `AMT24_InLine0_ACS*_Particulate_v20240610.sb`. These are the same Jordan et al.
  deposit we already have from Zenodo, in SeaBASS format. Download requires a
  SeaBASS account, which I do not have and did not attempt.

#### 3. Why I am confident the training `Rrs` is not public

Three independent lines of evidence:

**(a) The complete SeaBASS AMT holdings contain no radiometry — on any cruise.** I
enumerated `archive/PML/AMT/<cruise>/archive/` for all nine cruises:

| cruise | files |
|---|---|
| AMT19 | HPLC, AC-9, AC-S particulate |
| AMT22 | HPLC, AC-S particulate, `amt22_seabass.sb` |
| AMT23 | HPLC, AC-S particulate |
| **AMT24** | **HPLC, AC-S, AC-S2 particulate — no radiometry** |
| AMT25 | HPLC, AC-S, AC-S2 particulate |
| AMT26 / AMT27 | HPLC, AC-S particulate |
| AMT28 | HPLC, AC-9, AC-S particulate |
| AMT29 | (empty) |

Every file is dated `v2024xxxx`, i.e. the Jordan et al. 2025 deposit. PML has never
put AMT radiometry in SeaBASS.

**(b) Only seven AMT24 datasets are DOI-published worldwide, and none is
radiometry.** From a DataCite query for "AMT24" (7 hits, all seven inspected): the
cruise report, flow-cytometry counts, micro-molar nutrients, CDOM, CTD profiles,
nitrous oxide, and POC. That is the complete published record for this cruise.

**(c) BODC recorded the measurement but never published it.** The BODC cruise
inventory for JR20140922 explicitly lists *"Optical measurements of surface water
(remote-sensing reflectance). Measurements taken once per minute"* and *"Continuous
underway measurements of total upwelling radiance, downwelling sky radiance, and
total downwelling irradiance."* So the data exist. But of the 167 published BODC
series linked from that report, **all are CTD casts** — none is radiometry (I checked
every series row for radiance/irradiance/reflectance/optics keywords; zero hits).
The inventory documents what was *measured at sea*, not what was *published*.

Also checked and found nothing: `amt-uk.org/data/`, the NERC data catalogue, PANGAEA,
Zenodo beyond the Jordan deposit, and the data-availability statements of Lange et al.
2020 (it cites only the BODC cell counts — there is no Rrs statement) and of the
NTRS preprint.

#### 4. Who actually holds it, and the recommended ask

The cruise report is specific. The instrument was a **Satlantic HyperSAS with
136-channel HyperOCR** radiance and irradiance sensors, and the "Optical properties"
team was **Giorgio Dall'Olmo, Jelizaveta Ross, Rafael Rasse Boada and Tim Smyth
(PML)**. The report states the reflectance *"is to be used for calibration and
validation of satellite ocean colour products"* — i.e. it was always intended for
cal/val, just never deposited. Notably, **Priscila Lange (then Oxford) was aboard
AMT24**, so the MOANA first author has first-hand access.

So the ask is a data request to PML — Dall'Olmo and/or Smyth — and it folds naturally
into the §11 contact list already drafted in the report, where item 11.1.4 asks the
MOANA authors for the training set. Two routes to the same file; Lange is probably
the faster one.

*(One caution worth stating: even if we obtain "the AMT24 Rrs", Lange et al. applied
a specific processing chain — their own glint minimisation over 750–800 nm, tilt and
azimuth screening, the second-derivative noise filter, and interpolation to 2 nm.
Raw or differently-processed HyperSAS data will not reproduce their PCA exactly. We
should ask for the processed 414–660 nm / 2 nm matrix they actually fed to `prcomp`,
not just "the radiometry".)*

#### 5. Consequence for the plan

Validation targets **(i)** reproducing Lange Tables 1–2 and **(ii)** applying the
coefficients to held-out AMT cruises both need `Rrs` we do not have, so both are now
formally blocked on an external request. Target **(iii)** — validating the operational
PACE product — is unaffected and remains fully unblocked, as does the bit-exactness
check. That reinforces the ordering already in the report: implement MOANA, verify it
against a PACE granule, and pursue AMT in parallel with the data request outstanding.

### Prompt 3 (2026-08-01) — Answers reviewed; decisions locked in; access verified

All twelve answers reviewed and acted on. This section records the decisions so
prompts 4–9 don't have to re-litigate them, reports what I verified, and flags what
opening a real NASA granule turned up.

#### Decisions now locked in

| # | Decision | Consequence for later prompts |
|---|---|---|
| 1 | Get the loadings from OCSSW | **Done** (prompt 2). No email needed. |
| 2 | Dig for the AMT24 Rrs in SeaBASS/BODC | **NOT DONE.** A delegated search was started (SeaBASS `PML/AMT`, BODC/NERC, PANGAEA/Zenodo, plus the data-availability statements of Lin+2022 / Pardo+2023 / Brewin+2016) but the agent died before reporting, so **no findings were obtained** — nothing here to carry forward. It got as far as reverse-engineering the SeaBASS `search_results` POST form (`searchType=bio`, `bio_products_custom_text1=<field>`), which is a useful starting point for the retry. To be done under prompt 5 (Data), where it belongs |
| 3 | MOANA is a **separate track** from the IOP work, ending in its own report | No `AlgorithmSpec`, no BING wrapper, no shared metrics table. Code goes in `ioptics/moana/`, data in `ioptics/data/moana/` |
| 4 | Validation targets **(i) AMT24 reproduction, (ii) held-out AMT cruises, (iii) operational PACE product** | (iii) is unblocked *now*; (i) and (ii) wait on #2 |
| 5 | **Skip L23** entirely for MOANA | No Hydrolight path to build |
| 6 | Read `picophyt.json` rather than assume | **Done** (prompt 2) — zero-padding confirmed |
| 7 | Implement **414–660 nm @ 2 nm** | Confirmed by the LUT itself (124 bands) |
| 8 | Standardize with the **sample (N−1)** sd | Matches OCSSW |
| 9 | Implement **both** PC mappings behind a flag; settle empirically; JXP contacts the authors later | Default to the operational mapping so we can reproduce NASA bit-for-bit |
| 10 | Earthdata credentials available | **Verified working** |
| 11 | Keep **raw floats + a QC flag**, not NASA's silent clipping | Our reader must not truncate to int except in bit-comparison mode |
| 12 | Vendor the two LUT files into the repo | **Done** — `ioptics/data/moana/` |

#### Q12 — done

`ioptics/data/moana/` now holds `pca_picophyto.h5` and `picophyt.json`, checksums
matching NASA's manifest, plus a `README.md` documenting provenance, the file
layout, the slot→PC mapping, the known ATBD/LUT discrepancy, and the algorithm
recipe. The repo is now self-sufficient: MOANA can be run without OCSSW installed.

#### Q10 — Earthdata access verified, and the bit-exactness test is de-risked

`earthaccess` 0.18.0 (JXP installed it) authenticates against `~/.netrc`. What I
confirmed about the data we need:

- **NASA's MOANA output:** `PACE_OCI_L4M_MOANA` **v3.2**, daily / 8-day / monthly,
  at 4 km and 0.1°. A daily 0.1° granule is only ~2 MB.
- **There is no L3M MOANA collection in CMR** — only L4M and L4M_NRT — despite the
  `ob-cloud-pace-oci-l3m-moana-3.1` catalog page. Context corrected.
- **The Rrs input:** `PACE_OCI_L3M_AOP` **v3.2**, variable `Rrs`, **172 OCI
  wavelengths from 346.0 to 719.3 nm**. That fully brackets 414–660 nm, so MOANA
  never extrapolates on real OCI data. Spacing is irregular (~1.2–2.5 nm), which is
  exactly why the algorithm interpolates onto its own regular 2 nm grid. A daily
  0.1° granule is ~173 MB, and its header can be read by streaming without
  downloading the whole file.
- **Grids are co-registered.** All 1400 MOANA latitudes and 1100 longitudes appear
  *exactly* in the global AOP grid — MOANA 0.1° is a clean spatial subset
  (lat ±69.95, lon −84.95…24.95, i.e. the Atlantic sector). Alignment is by
  coordinate value; no regridding.
- **SST is available** (this was the unanswered half of Q&A #5):
  `CMC0.1deg-CMC-L4-GLOB-v3.0` — GHRSST CMC L4 at **0.1°, the same resolution as
  the MOANA grid** — is searchable and has granules for our test date. No `ocpy`
  loader needed; `earthaccess` reaches it directly.

**Key de-risking point:** *Synechococcus* and picoeukaryotes need **no SST at all**.
So we can run the full bit-exactness test — interpolation, standardization, scores,
and two of the three regressions — from Rrs alone, and that alone settles the
disputed *Synechococcus* PC index (operational PC16 vs ATBD U13). Only
*Prochlorococcus* needs the SST field, and its disputed index (PC7 vs U17) can be
tested afterwards. We are not blocked on SST for the discrepancy question.

Downloaded one reference granule to `$OS_COLOR/PACE/`:
`PACE_OCI.20250701.L4m.DAY.MOANA.V3_2.0p1deg.nc`.

#### What opening a real granule revealed — four product-level defects

I inspected that granule rather than trusting the catalog. Composition of the
1400×1100 grid:

| category | pixels | share |
|---|---|---|
| land sentinel `254` | 450,592 | 29.3 % |
| `_FillValue` = −32767 | 906,297 | 58.9 % |
| `INT32_MIN` = −2147483648 | 6 | — |
| **real ocean retrievals** | **183,105** | **11.9 %** |

**(a) Land is encoded as `254 cells mL⁻¹`.** Not as `_FillValue`, and worse, `254`
sits *inside* the declared valid range `[0, valid_max]`. Anyone who screens the
honest way — `valid_min ≤ x ≤ valid_max` — silently ingests 450,592 land pixels as
though *Prochlorococcus*, *Synechococcus* and picoeukaryotes were each present at
254 cells mL⁻¹. It nearly caught me: `254` is the modal *and* median value of all
three variables, at 71 % of everything that passes the declared range check. Land
identity confirmed by spot-testing the Congo basin, Sahara, Amazon, Andes and
interior Greenland (all 254) against four open-ocean points (fill or real values).
**Our reader must mask 254 explicitly.**

**(b) Six pixels carry `INT32_MIN`** in `prococcus_moana` — the signature of
`(int32_t)` casting a NaN or out-of-range float, which is undefined behaviour in C.
Not the declared fill value, and ~3500× outside the valid range.

**(c) The declared `valid_max` is not enforced.** On this single day:

| variable | valid_max | actual max | over by | pixels over |
|---|---|---|---|---|
| `prococcus_moana` | 600,000 | 706,915 | 1.2× | 41 (0.02 %) |
| `syncoccus_moana` | 300,000 | 6,594,758 | **22×** | 44 (0.02 %) |
| `picoeuk_moana` | 40,000 | 2,712,984 | **68×** | 930 (0.51 %) |

**(d) Silent clipping is common, not rare.** `prococcus_moana` is *exactly* 0 in
**18.0 %** of real ocean retrievals — that is the negative-clamp in
`get_Cpicophyt.c` firing on nearly a fifth of the product, with no flag to say so.
`syncoccus_moana` is 0 in **8.9 %**, which is the `int32` cast truncating genuine
sub-1 cell mL⁻¹ values to zero. This vindicates answer #11(b): had we copied NASA's
behaviour we would have silently thrown away, or log-transformed, a fifth of the
Prochlorococcus field.

Where the retrievals *are* real they look physically sensible — median 202,098
cells mL⁻¹ for *Prochlorococcus*, 5,145 for *Synechococcus*, 2,963 for
picoeukaryotes — so these are edge-case defects, not a broken product. But (a) and
(d) directly threaten any accuracy assessment, which is the point of this project.

Two smaller documentation mismatches: the in-file variable names are the misspelled
`prococcus_moana` / `syncoccus_moana` / `picoeuk_moana`, **not** the
`prochlorococcus_moana` / `synechococcus_moana` the Earthdata catalog advertises;
and the file's `title` says "OCI Level-3 Standard Mapped Image" while its
`product_name` says `L4m`.

### Prompt 2 (2026-08-01) — Existing MOANA code: located, downloaded, and decoded

**Bottom line: the blocker from Q&A #1 is gone.** The MOANA implementation is
public in full, and the two data files that neither Lange+2020 nor the ATBD
publishes — the PCA loadings and the regression coefficients — are now on disk
with checksums verified against NASA's own manifest. We can reimplement MOANA
bit-for-bit. Decoding those files also turned up a **substantive inconsistency
between NASA's documentation and NASA's operational code** (§4 below).

#### 1. Where the code lives

There is exactly one implementation, in C, inside OCSSW (the Ocean Color Science
Software that also backs SeaDAS). It is not on GitHub — `github.com/nasa/ocssw`
does not exist — and the Ocean Color doxygen browser has been retired (every
`oceancolor.gsfc.nasa.gov/docs/ocssw/…` URL now 301s to the OB.DAAC landing
page). The distribution route is the `install_ocssw` script.

Installed **OCSSW tag V2026.4** (newest of 30 available) to **`$OS_COLOR/ocssw`**:

```
$OS_COLOR/ocssw/
├── install_ocssw, manifest.py       # fetched from oceandata.sci.gsfc.nasa.gov/manifest/
├── ocssw_src/src/l2gen/             # ← the C source, incl. get_Cpicophyt.c
├── share/common/                    # ← the MOANA LUT + coefficients + product.xml
├── opt/                             # 3rd-party sources (pulled in by --src)
└── install.log
```

Command used (`--src` also pulls `opt_src` and `python`; `--common` is the share
bundle holding the MOANA tables):

```
python install_ocssw --tag V2026.4 --install_dir $OS_COLOR/ocssw --src --common
```

~11.5 GB total, of which `share/common` is 10 GB. The first run died partway
through `common` on a read timeout from `oceandata.sci.gsfc.nasa.gov`; re-running
with `--common` resumed and completed. **Note the two files we actually care about
are 24 KB and 390 bytes and are individually addressable** — no 11 GB install is
needed to get them:

```
https://oceandata.sci.gsfc.nasa.gov/manifest/tags/T2023.31/common/pca_picophyto.h5
https://oceandata.sci.gsfc.nasa.gov/manifest/tags/T2023.31/common/picophyt.json
```

(The bundle manifest at `…/tags/V2026.4/common/manifest.json` maps each file to
the tag that last changed it — both MOANA files are pinned at `T2023.31`, i.e.
they have not been touched since 2023 and predate the ATBD.)

MOANA-related source files:

| File | Role |
|---|---|
| `ocssw_src/src/l2gen/get_Cpicophyt.c` | the whole algorithm (279 lines; author Minwei Zhang, created 2023-09-29) |
| `ocssw_src/src/l2gen/l2prod.h` | product catalogue ids: `CAT_prochlorococcus` 353, `CAT_synechococcus` 354, `CAT_autotrophic_picoeukaryotes` 355 |
| `ocssw_src/src/l2gen/prodgen.c:1090` | dispatch — all three ids call `get_Cpicophyt` |
| `ocssw_src/src/l2gen/l12_proto.h:242` | prototype (note: output buffer is `int32_t*`) |
| `share/common/pca_picophyto.h5` | PCA loadings **V** |
| `share/common/picophyt.json` | regression coefficients |
| `share/common/product.xml` | official units / valid ranges / types |

No other implementation exists publicly: no Python or R port, nothing in
`nasa/oceandata-notebooks`, and Lange's original R analysis code was never
released.

#### 2. The two data files, decoded

**`pca_picophyto.h5`** — 24,488 bytes, sha256 `478ace87…deb57`:

- `component` — **(124, 45) float32**: the loading matrix `V[λ, i]`, 45 retained PCs.
  Columns are **orthonormal** (unit norm; `max|offdiag(VᵀV)| = 1.2e-7`), so this is
  a clean truncated eigenvector basis.
- `wavelength` — **(124,) int32 = 414, 416, …, 660 nm**. This settles Q&A #7 as a
  matter of fact, not preference: the operational grid *is* 414–660 nm at 2 nm.
- PC1 is a smooth blue→red ramp (max +0.171 at 414 nm, min −0.097 at 658 nm;
  correlation with wavelength −0.97), matching the paper's description of PC1 as
  backscatter slope plus water absorption.

**`picophyt.json`** — 390 bytes, sha256 `cd6b5c1c…756d4`. Verbatim:

```json
{ "npc": "17",
  "pro_coef":   "-4562144, 770448, 343590, 54975, 0, 0, 0, 290751, -1041404, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0",
  "syn_coef":   "-17.52889, 1.95949, 0.33859, 0.49127, 0, 1.04400, 0, 0, -1.06453, -1.55528, -2.25345, -3.41758, -1.67013, 0, 0, 0,-3.19241, 0",
  "apeuk_coef": "3.31202, 0, 0.21332, 0.09385, 0.72960, 0.97040, 0, 0, 0, 0, 0.67274, 0, 0, -2.25057, 0, -3.45864, 0, 0" }
```

**Q&A #6 is answered: yes, the coefficient vectors are zero-padded.** Their
lengths (19, 18, 18) match the C allocations `malloc(npc+2)`, `malloc(npc+1)`,
`malloc(npc+1)` with `npc=17` exactly, so the padding is by design, and
coefficient slot *k* maps to a fixed PC index.

#### 3. The algorithm as actually implemented

Per pixel, from `get_Cpicophyt.c`:

1. **Drop invalid bands.** Build `(wave_valid, Rrs_valid)` from bands where
   `Rrs != BAD_FLT`. If none are valid → `BAD_INT`.
2. **Linearly interpolate** `Rrs` onto the 124-point 414–660 nm grid. Indices are
   clamped at both ends, so out-of-range points are **extrapolated**, not flagged.
3. **Standardize the spectrum**: `Rrs' = (Rrs − mean)/sd` across those 124 values,
   with `sd` using the **sample (N−1)** denominator (Q&A #8 confirmed — this is
   per-spectrum standardization over wavelength, *not* feature-wise centering on a
   training mean, so **no training mean is needed anywhere**).
4. **Scores** for all 45 PCs: `U_i = Σ_λ V[λ,i]·Rrs'[λ]` — a plain dot product.
5. **Regressions**, using only the first `npc_used = 17` scores:
   ```
   Pro          = pro_coef[0] + pro_coef[1]·log10(SST) + Σ_{i=0..16} pro_coef[i+2]·U_{i+1}
   log10(Syn)   = syn_coef[0]                          + Σ_{i=0..16} syn_coef[i+1]·U_{i+1}
   log10(peuk)  = apeuk_coef[0]                        + Σ_{i=0..16} apeuk_coef[i+1]·U_{i+1}
   ```
   Syn and picoeuk are then exponentiated; **Prochlorococcus is not** (it was fit
   linearly, per the paper).
6. **SST** = `l2rec->sst` if `proc_sst` is on, else the reference field
   `l1rec->sstref` (the ATBD says that reference is GHRSST CMC L4).
7. **Clamp and cast**: a negative Pro becomes 0 when Syn and picoeuk are both
   positive; all three are then **truncated to `int32`**.

A self-consistency check confirms the loading file and the coefficient file belong
together. Because `Rrs'` is standardized, `‖Rrs'‖₂ = √123 = 11.09` exactly, and
because `V`'s columns are orthonormal, **every score is bounded by ±11.09**. Using
that ceiling, Prochlorococcus evaluates to 1.5×10⁵ (SST 15 °C) to 3.6×10⁵ cells/mL
(SST 28 °C) — physically correct for *Prochlorococcus*, and inside the official
`validMax` of 6×10⁵. It only lands in range when `U₁` sits near its ceiling, which
is precisely what >96 % of covariance on PC1 implies. The huge Pro coefficients
(~10⁶) are therefore real, not a transcription error.

#### 4. Inconsistencies found

**(a) The serious one: the ATBD equations and the operational coefficients disagree
about two PC indices.** Decoding `picophyt.json` through the C indexing and
comparing against the ATBD v1.2 equations:

| Taxon | PCs per operational `picophyt.json` | PCs per ATBD v1.2 equations |
|---|---|---|
| Prochlorococcus | 1, 2, 6, **7** | 1, 2, 6, **17** |
| *Synechococcus* | 1, 2, 3, 5, 8, 9, 10, 11, 12, **16** | 1, 2, 3, 5, 8, 9, 10, 11, 12, **13** |
| picoeukaryotes | 2, 3, 4, 5, 10, 13, 15 | 2, 3, 4, 5, 10, 13, 15 ✓ |

Every coefficient *value* matches (to ATBD rounding: `0.33859`→`0.3385`,
`-3.19241`→`-3.1924`). Picoeukaryotes agree completely. But for Prochlorococcus and
*Synechococcus* the **last** coefficient sits on a different principal component in
the file that runs than in the document that describes it.

Evidence on which is right: the union of PCs across all three taxa is **14 distinct
PCs for the ATBD** version versus **15 for the JSON** — and Lange+2020 states that
exactly **14** PCs survived the significance cut for the hyperspectral PCA. The
per-taxon counts (4 / 10 / 7) match the paper's prose either way. So the ATBD
equations reproduce the paper's PC budget exactly and the JSON does not, which
points to `picophyt.json` having two coefficients in the wrong slots — meaning the
**operational PACE MOANA products may not be computing the published algorithm.**
I would not assert that as fact without the authors confirming; the alternative is
that the ATBD subscripts are typos and the paper's "14" counts something slightly
different. Either way it is worth an email, and our implementation should support
both mappings so we can quantify how much the retrievals actually differ.

**(b) The ATBD contradicts itself on spectral range** — abstract says Rrs must be
valid over 395–705 nm, Mathematical Theory says 414–660 nm. The LUT settles it:
414–660 at 2 nm.

**(c) The ATBD says "only 14 PCs were selected"** but its own equations index up to
U17, and the LUT ships 45 components with `npc_used = 17`. These are reconcilable
(14 *used*, 17 *scanned*, 45 *stored*) but nothing in the documentation says so.

**(d) Product short names are misspelled and inconsistent across levels.** L2 uses
`prococcus` / `syncoccus` / `picoeuk`; the L3/L4 files use `prochlorococcus_moana` /
`synechococcus_moana` / `picoeuk_moana`. Any reader we write must handle both.

#### 5. Quirks that will bite a reimplementation

- **Integer output.** `product.xml` declares `<type>int</type>` and the code casts
  with `(int32_t)`, which **truncates rather than rounds**. Sub-1 cell/mL effects
  only, so minor — but it means a bit-exact match to NASA's product requires
  truncating too.
- **Negative Prochlorococcus is silently zeroed, not flagged.** Since Syn and
  picoeuk come from `pow(10, …)` and are always positive when valid, the guard
  reduces to "clip Pro at 0". Zeros cannot be log-transformed, so this materially
  affects any log-space validation statistic — we must decide whether to treat
  clipped zeros as data or as flags.
- **Two dead-code paths**, harmless but confusing when reading: `bindx` is computed
  via `windex` and never used, and the `< 0` tests on `pow(10, …)` results can
  never fire.
- **A latent out-of-bounds read.** With very few valid bands the interpolation
  clamp can select `j = 1` while only index 0 was filled, reading uninitialized
  `wave_valid`/`Rrs_valid`. There is no minimum-valid-band check. Harmless for OCI
  in practice (hundreds of bands) but our version should require a sane minimum.
- **Extrapolation is silent** at the grid edges rather than masked.
- MOANA also runs under `l3gen` (there is an explicit `npix` realloc for it), which
  is consistent with the operational products being L3/L4 mapped rather than L2.

#### 6. Official product metadata (from `share/common/product.xml`)

| L2 name | units | type | validMin | validMax | display |
|---|---|---|---|---|---|
| `prococcus` | cells ml⁻¹ | int | 0 | 6×10⁵ | 1–5×10⁵, linear |
| `syncoccus` | cells ml⁻¹ | int | 0 | 3×10⁵ | 1–1.5×10⁵, log |
| `picoeuk` | cells ml⁻¹ | int | 0 | 4×10⁴ | 1–3×10⁴, log |

All three cite Lange et al. (2020) as the reference.

#### 7. What is still missing

Only one thing, and it is not code: **the AMT24 training data** — the
hyperspectral Rrs in particular (Q&A #2). We can now *run* MOANA exactly as NASA
does; we cannot yet *retrain* it or reproduce Lange's Tables 1–2 without the
training spectra. Also absent: any uncertainty estimate (MOANA emits none), and any
published validation of the satellite products.

#### 8. Recommended next step

Before touching AMT24, do a **bit-exactness check against NASA**: download one PACE
OCI L3M Rrs granule plus the matching L3M MOANA granule and confirm our Python
implementation reproduces NASA's three fields pixel-for-pixel. That validates our
transcription against the operational truth in one shot, and it is the cleanest way
to settle §4(a) empirically — if we reproduce NASA's product with the JSON mapping
but not the ATBD mapping, the operational code is definitively using PC7/PC16.

### Logs

### 2026-09-10 (Prompt 20 — Tom's reply reviewed: mapping confirmed, but it answers a different question)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Tom is right, and it is already what we do.** He wrote that the cell counts are
the final columns of `AMT24_JR20140922_AFC_Dataset.xlsx` and that the codes map
onto `AMT24_AFC_BODC_sample_metadata_report.html`, e.g. `P700A90Z` =
*Synechococcus*. Verified on receipt rather than taken on trust:

- the metadata report defines exactly **six** codes, and the workbook carries
  exactly those six taxon columns — nothing hidden, nothing extra;
- our reader's `FCM_CODES` is identical code for code, including
  `P700A90Z → syn` and `P701A90Z → pro`;
- the values our loader returns match the workbook **cell for cell**
  (`np.allclose` per column; 810–812 finite samples each).

So there is no decoding problem and no column we were failing to read. The request
has been misheard as *"we cannot find the cell counts"* when it is *"we cannot find
the **underway** cell counts"*.

**I turned the gap into numbers, since prose has now failed twice.** Straight from
the deposit: 814 bottles, all `Gear=CTD`, 68 stations, 68 surface (≤ 10 m) samples,
their times bimodal at 03–06 h and 12–15 h UTC — i.e. the two daily rosette casts,
which is exactly what a bottle product should look like. The decisive row is the
front band: **in 25–45° S, where Lange et al. §2.1 and their Fig. 2b describe the
30-minute underway run, this deposit has 16 surface samples with a median spacing
of 14.9 hours (1.11°, ~66 nM).** That is the whole discrepancy in one line, and it
is checkable by anyone with the file.

**Deliverable is a re-ask, worded to stop the loop.** Added to
`requests/PML_follow_up.md` §1: the verification above, a comparison table
(sampling mode / surface count / time of day / front-band coverage / usable n),
and draft text that **opens by confirming we have the columns and read them
correctly** before asking the actual question. The previous phrasing evidently
invited a re-explanation of the codes; leading with "thanks, we do have those and
they match" should prevent a third round. It also repeats the two points most
likely to unstick it: that a bottle deposit *structurally cannot* contain underway
samples (`BODC_bot`/`Rosette_Pos`/`Firing_Seq`), and that **a plain spreadsheet is
sufficient — no archiving needed**.

**One genuine gain.** The code→taxon mapping is now confirmed by the data
originators rather than inferred by us from the metadata document. That is a real
provenance upgrade for every AMT cruise we score, so I recorded it at
`ioptics/moana/io.py:FCM_CODES` with the quote and date. Small, but it is the kind
of thing that is annoying to re-derive later.

**No new analysis.** The deferred underway retraining remains blocked for exactly
the reason logged on 2026-09-07 — nothing in this exchange changed the data on
disk. `pytest -q`: **300 passed** (the io.py change is a comment, but the loader is
on the tested path so it was worth confirming).

**Files changed.** `requests/PML_follow_up.md` (§1 update + re-ask),
`ioptics/moana/io.py` (provenance note on `FCM_CODES`),
`claude_prompts/moana_prompts.md` (this entry).

### 2026-09-07 (Prompt 19 — PML re-delivery examined: same CTD-bottle data, so the deferred analysis stays blocked)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Verdict: the new deposit is not what we need, so I did not run the deferred
analysis.** Prompt 19 was conditional ("if they are what you need, then
proceed"), and the condition is not met. Reporting that rather than forcing the
analysis on data that cannot support it.

**What arrived.** `$OS_COLOR/AMT/AMT24/AMT24_JR20140922_6param/` — `AMT24_AFC.out`,
`AMT24_JR20140922_AFC_Dataset.xlsx`, a BODC sample-metadata report and a flag key.
The "6param" refers to the six taxon columns (Pro, Syn, picoeuk, cryptophytes,
coccolithophores, nanoeukaryotes).

**What I checked, rather than assuming from the filenames.** It is the *same*
data we already had:

- the `.xlsx` holds a single sheet, `botlist`, 814 rows × 29 columns — I opened it
  specifically because a workbook can hide extra sheets, and an underway table
  would plausibly have been one;
- 68 stations, **every row `Gear=CTD` and `ODV_type=b`**;
- joined on the `BODC_bot` bottle id, **all six taxon columns are numerically
  identical** to `AMT24_JR20140922_AFC_Dataset.csv` (`np.allclose` on each);
- only the container and the timestamp format differ (ISO vs `dd/mm/yyyy`).

**The structural point, which is the useful part.** This is not "they sent the
wrong file" — it is that the underway samples *cannot* live in this deposit. Its
own metadata calls it a "Bottle Retrieval Data Report", and every row is keyed by
`BODC_bot` with `Rosette_Pos` and `Firing_Seq`. An underway sample has no bottle,
no rosette position and no firing sequence, so it is unrepresentable in this
series. Re-requesting the DOI will therefore never help; the data need a separate
BODC series or a direct file from Glen Tarran. I re-checked DataCite on the day
(three query phrasings) and there is still exactly **one** AMT24 flow-cytometry
DOI — the bottle one — so there is no second deposit sitting unfound.

**Why it matters, in one number.** Our full-fit *Synechococcus* MAE on the
CTD-only configuration is 1.36, essentially identical to Lange's own CTD-only
figure of 1.37, against 1.27 with the underway samples. The missing data *is* the
gap between our reproduction and theirs — which is a good argument for chasing it
and a bad argument for pretending 30 samples can stand in for 73.

I also made explicit in §12.2 something that was previously only implied: n = 30
is a **ceiling, not a matchup-tuning artefact**. The only AMT24 FCM data in
existence is bottle-only at 68 stations, and fewer than half those casts have
daylight radiometry. No reprocessing on our side can reach Lange's n.

**Updates made.** `requests/PML_follow_up.md` §1 gains a dated update recording
the verification, the structural explanation, and a **sharpened ask**: the
~30-minute underway run across 25–45° S (Lange §2.1, their Fig. 2b), three taxa in
cells mL⁻¹ with UTC timestamp and lat/lon, and an explicit note that **a
spreadsheet is sufficient — it need not be a BODC deposit for us to use it**. That
last point may be what unsticks this, since the previous phrasing could have been
read as asking for a formal archive submission. `reports/MOANA_Claude_Report.md`
§13 item 1 records the re-delivery check so nobody re-treads it, and §12.2 gains
the ceiling note.

**No `pytest` run:** no code touched this pass — the changes are two markdown
documents.

**Files changed.** `requests/PML_follow_up.md`, `reports/MOANA_Claude_Report.md`,
`claude_prompts/moana_prompts.md` (this entry). No new data files; nothing
downloaded.

### 2026-09-07 (Prompt 18 — two new report figures: what the coefficient bug costs, and held-out skill drawn)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**No workstation run needed.** Prompt 18 offered that possibility; it turned out
everything was already local. My first pass looked for the data before JXP had
finished reorganising it and I wrongly concluded two pieces were missing — the
AMT23/25/28 flow-cytometry tables and the L3M AOP Rrs granule. Re-checking after
his prompt found both: the FCM deposits had moved to `$OS_COLOR/AMT/AMT<n>/`
(where `load_fcm` looks for them, so all four cruises load cleanly) and the 181 MB
AOP granule was in the `moana_validation` cache. Everything ran here in ~10 s.

**Figure 1 — `moana_mapping_consequence.png`, added to §7.1.** Deliberately a
*different* question from the existing `moana_mapping_verdict.png`, which asks
which mapping NASA ships. This one asks what the misplacement costs, by running
our own retrieval under both mappings on the same granule so the only difference
is where the two disputed coefficients sit.

The methodological point I was most pleased with: **the *Prochlorococcus*
difference is exactly computable even though we have no GHRSST field.** The
intercept and the `log₁₀(SST)` term are identical under both mappings, so they
cancel in the subtraction and only the relocated PC term survives — meaning a
figure I had assumed would need an ancillary download needed none. I pass a
constant SST purely to exercise the public `run_moana` path, and `nasa_compat=False`
so neither the negative-Pro clamp nor the int32 truncation contaminates the
difference.

Results: *Synechococcus* under the operational mapping is a systematic **0.80×**
the as-published value (median Δlog₁₀ = −0.098, n = 166,733); *Prochlorococcus*
moves by a median **54 × 10³ cells mL⁻¹, which is 21 % of NASA's own retrieved
value**; and **picoeukaryotes come out bit-identical**, which is the useful part —
its assignment is undisputed, so agreement there is a control confirming the
difference is NASA's table and not our code.

The finding worth the section text: the *Prochlorococcus* difference **changes
sign regionally** — negative through the gyres, positive at the high-latitude
margins and the equatorial band. So it is not a calibration offset a user could
absorb into a scale factor; it distorts spatial gradients, across the very fronts
the algorithm exists to detect. I scaled the difference to NASA's own Pro field
rather than leaving it as raw cells mL⁻¹, because "21 % of the value" is
interpretable and "54,000 cells" is not.

**Figure 2 — `moana_heldout_skill.png`, added to §12.3.** Retrieved vs observed,
one panel per taxon, coloured by cruise, with 1:1 and ±3× guides and the
Lange+2020 MODIS held-out numbers quoted in each panel for reference. It draws
the table's story: picoeukaryotes on 1:1, *Synechococcus* fanning, and
*Prochlorococcus* **flat** — retrieved 1–4 × 10⁵ cells mL⁻¹ almost regardless of
what was counted, across nearly four decades of truth. That flatness is the
negative R² made visible, and it is §7.3's clipping seen from the other side: the
model cannot express low Pro, so when the truth is low it neither tracks it nor
flags it.

**One small library change.** `validate_heldout_cruises` computed the matchup
pairs and then discarded them, keeping only the metrics. Rather than duplicate the
±3 h matching logic in a report script — which would have been a second place for
it to drift — I added a `pairs` DataFrame to its return. Backwards-compatible
(new key only), documented, and the metrics it reproduces are unchanged: the run
matched §12.3's published table exactly, cruise by cruise, before I touched
anything.

**Mistakes I made and fixed.** Two worth recording. First, my margin patches for
the held-out figure kept landing on `plot_mapping_consequence` instead, because I
anchored the replacements on a `tight_layout(rect=...)` string that appears in
several functions and only asserted that *some* replacement had happened. Three
rounds of "fixed it" produced byte-identical output before I checked the function
body and found the real cause; I then re-anchored by line number and verified
every layout call against its owning function. The lesson is that
`assert s != orig` is too weak a guard when a patch makes several edits — it
passes if any one lands. Second, the equal-aspect log panels genuinely do defeat
`tight_layout`'s `rect` (and `bbox_inches="tight"` did not rescue it), so that
figure now uses explicit `subplots_adjust`.

I also fixed a **pre-existing bug** while regenerating: running
`python reports/scripts/moana_report_figs.py` directly failed with
`ModuleNotFoundError: ioptics`, because `sys.path[0]` is the script's own
directory, not the repo root. It had only ever worked when driven from an
interpreter already rooted in the repo. The script now inserts `_REPO` on
`sys.path`, so the whole set — including the older verdict figure, which had the
same latent problem — regenerates from one command.

**Verification.** Loaded the house data-viz guidance before writing chart code.
Magnitude stays on the single-hue blue sequential ramp; the two signed-difference
maps use the palette's diverging blue↔red with a neutral grey midpoint and
symmetric limits at the 95th percentile of |Δ| (the red arm is stepped to mirror
the blue arm's lightness, since `palette.md` tabulates only the blue); cruise
identity uses the first three categorical slots, which are documented as
validated all-pairs — the right gate for a scatter, where every pair is
co-visible. Still no JS runtime here, so the validator could not be run; I relied
on the documented validation of that exact subset rather than eyeballing, as
before. Rendered and inspected every figure, which is what caught the axis and
legend faults. `pytest -q`: **300 passed**.

**Files changed.** `reports/scripts/moana_report_figs.py` (two figure functions,
diverging ramp, `sys.path` fix, wired into `main()`),
`reports/MOANA_Claude_Report.md` (§7.1 and §12.3 additions),
`ioptics/moana/validation.py` (`pairs` in the return),
`claude_prompts/moana_prompts.md` (this entry). New:
`reports/figures/moana_mapping_consequence.png`,
`reports/figures/moana_heldout_skill.png`.

### 2026-08-17 (Prompt 17 — report rev. 4: scale, decadal drift, new-basin training, and the AMT accounting)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

All four requested additions are in `reports/MOANA_Claude_Report.md` (rev. 4):

1. **Patchiness vs pixels → new §9.9.** The scale contradiction stated
   plainly: the Syn evidence rests on ~300 m observation footprints (1.6 mL
   bottles matched to 1-minute Rrs), the application is ~1.2 km OCI pixels
   composited to 4–11 km — and because two models are nonlinear in Rrs,
   retrieval-of-mean ≠ mean-of-retrievals (Jensen), with our §7.1 compositing
   result as one-level-up proof that the ordering matters. Three
   reconciliations proposed; the first (measure subpixel retrieval variance
   directly from the AMT24 1-minute stream, data already in hand) became Open
   item 7.

2. **Other AMT cruises → §12.3 gains the complete accounting table.** Did we
   try? Yes — AMT23/25/28 are done (that *is* §12.3). AMT26 is possible
   pending its FCM (Open item 6); AMT20/22 are impossible without a PML data
   rescue (their radiometry was never archived — verified in the prompt-6/7
   sweeps); and **no AMT cruise can validate the PACE-era product** because
   the archived FCM series ends in 2018 while PACE launched 2024-02 — product
   validation is the SeaBASS route.

3. **2014 training vs 2024+ application → new §9.10.** The concern is real
   and now has magnitudes: the Pro model's SST term gives
   dPro/dSST = 770448/(T·ln10) ≈ 16,700 cells mL⁻¹/°C at 20 °C (~7–8 % of a
   typical retrieval per degree). Mean Atlantic warming 2014→PACE-era
   (~0.2–0.3 °C) is a tolerable ~2 % inflation, but the 2023–24 North
   Atlantic heatwave (~+1 °C basin-wide) manufactures ~8 % (locally tens of
   %) of Pro signal from temperature alone — in exactly the years PACE
   observes. Ecological drift (poleward range shifts, ecotype changes, CDOM
   covariate drift) is discussed alongside; fixes: periodic retraining, a
   training-epoch attribute, and the §9.3 out-of-domain flag. A within-decade
   hint from our own §12.3 (AMT25 beats AMT28) is noted with the appropriate
   small-n caution. Quantifying the drift per-pixel became Open item 8.

4. **Other basins / Suzuki's ship-of-opportunity program → new §9.11.** The
   ingredient list (co-located hyperspectral Rrs + Pro/Syn/peuk FCM + SST at
   underway resolution), why the North Pacific SOOP route fits (repeat
   transects across the subarctic/Kuroshio/gyre boundaries; Suzuki's group
   already runs the HPLC+FCM toolkit — grounded in Suzuki et al. 2005), and
   the three questions that decide feasibility (radiometer aboard? taxa
   resolved? which lines repeat?). Complementary resources named: SeaFlow
   (27 NE-Pacific underway-FCM cruises, Ribalet et al. 2019), HOT/BATS, and
   simply *newer AMT cruises with radiometry archived*. Became Open item 9.
   Both new citations verified via Crossref before inclusion.

Header bumped to rev. 4; §13 renumbered (old item 7 → 10). I could not verify
the specifics of the Suzuki program's instrumentation from the open web (two
searches; his 2005 subarctic-Pacific FCM/HPLC work is solid ground, the SOOP
details are not), so §9.11 frames those as the questions to ask rather than
asserting them — flagging that here so it isn't mistaken for reticence.

**No `pytest` run:** report-only changes, no code touched. Nothing committed.
**Files changed:** `reports/MOANA_Claude_Report.md` (rev. 4),
`claude_prompts/moana_prompts.md` (this entry).

### 2026-08-16 (Prompt 16 — report rev. 3: results absorbed, §7.1 marked SETTLED, the promised verdict figure added, new §12 results + §13 Open items)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**`reports/MOANA_Claude_Report.md` updated to rev. 3.** The changes, in the
order a reader meets them:

- **Header/status** — the algorithm is now independently reimplemented and
  validated; pointers to §12 (results) and §13 (open items).
- **§7.1 gains a SETTLED box** with the granule-experiment numbers, the new
  figure, and the two structural corollaries (bit-exactness from L3M is
  impossible because NASA retrieves at L2 then composites; that ordering
  likely explains the missing L3M MOANA collection).
- **New §7.6** — the paper-vs-ATBD/LUT disagreement on the picoeukaryote SST
  term (found in prompt 10, previously only in the design doc).
- **§10 rewritten as a status board** — items 1–3 done (with their caveats),
  item 4 pending data; pointer added to `requests/PML_follow_up.md`.
- **§11 updated** — item 1 reframed (we now *know* what ships; the question
  is what was intended), item 4 updated for the PML delivery, new items 8
  (SST-for-peuk) and 9 (confirm the L2-composite ordering).
- **New §12 "Our reproduction and validation"** — the training-data
  characterisation (§12.1), the Table-1 reproduction with the basis-recovery
  result and the published-model transfer experiment (§12.2), the first-ever
  held-out hyperspectral skill table (§12.3), and the granule experiment
  summary (§12.4).
- **New §13 "Open items"** (the prompt's explicit ask): underway FCM,
  target iii-b (SeaBASS credentials pending), the mapping-intent question,
  the §9.3 out-of-domain experiment, Lange-strict + Eq.-7 sensitivity passes,
  AMT26 as a bonus cruise, and publication hygiene (dark-mode figures, the
  Q&A #22 product-reading recipe — both still open).
- **References** — the §12 data sources (PML delivery, four BODC FCM DOIs,
  Brewin 2023, Jordan 2025, the PACE granule pair).

**The promised third figure now exists** — `moana_mapping_verdict.png`, added
to `reports/scripts/moana_report_figs.py` (house palette, single-hue
sequential density, 1:1 baseline): two log-log density panels of our
Synechococcus against NASA's under each mapping. The operational panel sits on
the 1:1 line (median Δlog₁₀ −0.005); the ATBD panel is visibly displaced
(+0.082). The script's granule paths now also resolve the prompt-14
`earthaccess` cache, and all four figures regenerated cleanly.

**Also updated:** `ioptics/data/moana/README.md` — its "unconfirmed as of
2026-08-01" mapping caveat now records the empirical settlement.

**`pytest -q`: 300 passed.** Nothing committed. **Files changed:**
`reports/MOANA_Claude_Report.md` (rev. 3),
`reports/scripts/moana_report_figs.py` (+`plot_mapping_verdict`, path fixes),
`reports/figures/moana_mapping_verdict.png` (new; others regenerated),
`ioptics/data/moana/README.md`, `claude_prompts/moana_prompts.md` (this entry).

### 2026-08-16 (Prompt 15 — target (ii) scored: first-ever held-out numbers for hyperspectral MOANA; peuk transfers, pro/syn do not)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Data reviewed (Q&A #36 answer).** The three deposits are on disk under the
new `$OS_COLOR/AMT/AMT<c>/` layout — AMT23: 730 bottles / 54 stations, AMT25:
915 / 73, AMT28: 755 / 70 — all CTD-bottle-only, all sharing the AMT24 file
layout. The taxon↔code mapping was re-verified against **each cruise's own
metadata document** (not assumed from AMT24): P700A90Z = *Synechococcus*,
P701A90Z = *Prochlorococcus* everywhere. AMT28 carries two extra taxa codes
(ignored). One surprise handled: the reorganisation **moved** the old
`$OS_COLOR/AMT24/` tree into `$OS_COLOR/AMT/AMT24/`, which broke every loader
path — all AMT24 paths now resolve through a tolerant `_amt24_dir()` helper
(both layouts accepted), the conftest probe likewise, and `load_fcm` gained a
`cruise=` parameter.

**Target (ii) scored** — the combination with no published number:
in-situ hyperspectral Rrs (Brewin 2023) through the published NASA
coefficients, against flow-cytometry truth (±3 h, ≤10 m), 66–71 matchups
pooled over AMT23/25/28 (bias/MAE/R²):

| taxon | AMT23 | AMT25 | AMT28 | pooled | Lange held-out (MODIS) |
|---|---|---|---|---|---|
| pro | 1.20/1.59/−1.8 | 2.47/2.90/0.15 | 3.17/4.07/−0.4 | **2.19/2.75/−0.14** | 1.75/2.26/0.54 |
| syn | 1.26/1.77/0.71 | 1.08/2.05/0.54 | 0.48/3.90/−10.6 | **0.85/2.46/−3.0** | 0.93/2.20/0.40 |
| peuk | 0.79/1.40/0.76 | 1.16/1.58/0.77 | 1.24/1.43/0.80 | **1.07/1.47/0.78** | 1.05/1.53/0.60 |

The reading: **picoeukaryotes transfer** — pooled MAE 1.47 / R² 0.78, *better*
than Lange's own MODIS-based held-out line. **Prochlorococcus and
Synechococcus do not**, even with best-case hyperspectral in-situ radiometry —
pro is +119 % biased with negative R², and syn collapses specifically on AMT28
(bias 0.48, R² −10.6). Since the radiometry here is in-situ and hyperspectral,
the transferability failure is in the *model* (single-cruise training, SST
climatology, PC fragility), not in satellite reflectance error — the same
per-taxon pattern as target (i)'s published-model column. Unphysical fractions
are small on these oligotrophic-dominated transects (0–4 %).

**Mechanics:** pooled metrics folded into `validate_heldout_cruises` itself
(per-cruise + `'pooled'`); one new tier-2 test runs the scored path end-to-end.
**SeaBASS (target (iii)b):** noted your username `profx` — yes, I will need
the password for the archive's authenticated download; please add it to
`~/.bashrc` as offered and (iii)b proceeds next pass. Nothing else blocks.

**`pytest -q`: 300 passed.** Nothing committed. **Files changed:**
`ioptics/moana/io.py` (`_amt24_dir`, `load_fcm(cruise=)`),
`ioptics/moana/validation.py` (pooled metrics + verbose table),
`ioptics/tests/{conftest,test_moana}.py`, `claude_prompts/moana_prompts.md`
(this entry).

### 2026-08-16 (Prompt 14 — validation implemented AND run: Q&A #9 settled — NASA ships the operational mapping; Lange Table 1 reproduced in the CTD-only configuration)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Deliverable:** `ioptics/moana/validation.py` (the stub replaced) — Seegers
log-space metrics with censored handling and the headline unphysical fraction
(Q&A #13); `validate_amt24` (target (i): published model + retrained full-fit +
80/20 bootstrap CV, with Lange's Table 1 reference values shipped in the module
for diffing); `validate_heldout_cruises` (target (ii), on a new
`load_brewin2023` reader in `io.py`); `bitexact_pace` and `match_pace_to_insitu`
(target (iii), via `earthaccess` — installed into ocean14, already in
`requirements.txt`). Plus 4 new tests (metrics + two `needs_amt24` regressions);
suite now **299 passed**.

**Headline 1 — Q&A #9 is settled empirically: the shipping PACE product uses
the operational (`picophyt.json`) mapping.** On 100k ocean pixels of the
2025-07-01 daily 0.1° granule pair, our Synechococcus (no SST needed) tracks
NASA's at **median Δlog₁₀ = −0.005 (MAD 0.007)** under the operational mapping,
versus a **+0.084 median offset (MAD 0.090)** under the ATBD mapping — and the
two mappings genuinely differ on 98 % of pixels. Combined with the prompt-4
inference that the ATBD's 14-distinct-PC assignment is what the paper trained,
this now *demonstrates* the operational product does not implement the
published algorithm. One structural discovery en route: **true bit-exactness
from L3M inputs is impossible** — NASA retrieves at L2 and then composites
(retrieve-then-average ≠ average-then-retrieve; also the likely reason no L3M
MOANA exists), so the test's verdict runs in log space; only ~12 % of pixels
(presumably single-overpass) match exactly. Both findings belong in the report
(prompt 15) and the §11 NASA list.

**Headline 2 — target (i): Lange Table 1 reproduced, within the CTD-only
handicap.** n = 30 matchups (Lange: 73–78 incl. underway FCM). Our retrained
full-fit vs Lange Table 1 (bias/MAE/R²): pro 1.01/1.22/0.77 vs 1.08/1.31/0.82;
**syn 1.00/1.36/0.91 vs 1.00/1.27/0.92 — our MAE 1.36 lands on Lange's own
CTD-only sensitivity figure of 1.37**, exactly the predicted degradation;
peuk 1.00/1.19/0.97 vs 1.00/1.21/0.95. The bootstrap CV degrades sharply
(pro MAE 1.78, syn unusable) — honest small-n stepwise behaviour at n = 30.
The *published* NASA model applied to our matchups: peuk transfers beautifully
(1.03/1.23/0.94), pro (+96 %) and syn (−39 %) do not — consistent with those
two taxa riding the processing-chain details (glint treatment, SST source)
that differ between our chain and Lange's.

**Target (ii): retrievals done, scoring blocked on data.** 90 stations on
AMT23/25/28 retrieved from the Brewin 2023 in-situ Rrs (plain, non-BRDF family)
with physically sensible values (Pro medians 1.9–2.6×10⁵ cells mL⁻¹, 3–7 %
negative-Pro). The flow-cytometry counts for those cruises are **not on disk**
— posted as Q&A #36 with the three BODC links, alongside the SeaBASS-account
question for target (iii)b.

**One bug found and fixed during the runs:** `load_brewin2023` initially swept
all 453 `Rrs(`-prefixed columns into one "spectrum" — the file interleaves
*three* 151-column families (plain Rrs, BRDF-corrected Rrs, uncertainty-%),
which drove Pro to −5×10⁶ and Syn to 0. The loader now selects one family
explicitly (default: plain, matching the training geometry; `brdf=True`
available) and returns the uncertainties; a regression test pins it. Also two
of my new metric tests initially used degenerate inputs (constant obs → R²
undefined); fixed the tests, not the code. PACE granules cache under
`$OS_COLOR/PACE/moana_validation/` (173 MB AOP + 2 MB MOANA per day).

**`pytest -q`: 299 passed.** Nothing committed. **Files changed:**
`ioptics/moana/validation.py` (implemented), `ioptics/moana/io.py`
(`load_brewin2023`), `ioptics/tests/test_moana.py` (4 tests),
`claude_prompts/moana_prompts.md` (Q&A #36, this entry).

### 2026-08-16 (Prompt 13 — `ioptics/tests/test_moana.py`: 31 tests, full suite 295 green)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Deliverable:** `ioptics/tests/test_moana.py` — 31 tests implementing the
design doc's §8 plan — plus two new conftest markers following the house
two-tier pattern: `needs_amt24` (the AMT24 tree under `$OS_COLOR`) and
`needs_netrc` (Earthdata credentials; defined now per Q&A #16, first used by
validation target (iii) in prompt 14).

**Tier 1 (27 tests, run everywhere):** the vendored LUTs against their pinned
sha256 (a silent NASA LUT change now fails loudly) and the 414:2:660 grid;
the §3.1 invariants (mean 0 and norm √123 exactly, |Uᵢ| bounded, orthonormal
basis); standardisation checked against an independent hand computation of
Lange Eq. 3; interpolation flags; SST flagging; **the PC-mapping flag** (moves
exactly the two disputed slots, and the two mappings agree identically on
spectra built orthogonal to the disputed PCs); **nasa_compat** (clamp,
truncate-not-round); the reconstruction residual separating in-basis from
out-of-basis spectra; every pipeline stage on fabricated streams with planted
answers — geometry screens cutting exactly at Lange's 5°/10°/80°/50°/170°,
1-min selection taking the NIR-darkest, **the glint fit recovering planted
(ρ_sky, L_NIR) to 2×10⁻³**, each QC screen firing on its own planted defect,
the ±15 min median rejecting a planted glint spike where Lange-strict mode
(correctly) takes it; the training recipe recovering a planted linear model,
the saturation guard, the starved-taxon error, the per-taxon SST switch, and
basis self-comparison; and a **regression test for the pandas ≥2 time-
resolution trap** found in prompt 12.

**Tier 2 (4 tests, `needs_amt24`):** the real `.sav` layout + the ES≢LT guard
+ the fixed-ρ=0.0280 provider-Rrs identity; the BODC FCM mapping — including a
science check that would catch a P700/P701 swap instantly (tropical surface
Pro ≫ Syn) — and the shallowest-bottle rule; one real day through the full
chain (attrition monotone, no surviving negatives, ρ in bounds); the Jordan
SST series physical and monotone.

**One fix to my own tests:** the first version of the per-taxon-SST assertion
was a tautology (`'logSST' in {…, 'logSST'}`); rewritten to assert SST enters
the candidate pool when switched on and never appears under the operational
default.

**`pytest -q`: 295 passed** (264 existing + 31 new, ~2 min). With `$OS_COLOR`
unset: 27 passed, 4 skipped in 1.4 s — the tier-2 guards work. No new Q&A
questions. Nothing committed. **Files changed:**
`ioptics/tests/test_moana.py` (new), `ioptics/tests/conftest.py` (two
markers), `claude_prompts/moana_prompts.md` (this entry).

### 2026-08-16 (Prompt 12 — `ioptics/moana/` implemented; NASA's PC1/PC2 recovered from our own retraining at |cos| ≥ 0.998)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Deliverable:** the `ioptics/moana/` package, per the design doc — five modules,
methods not classes, every Lange threshold carried in one serialisable config
dict (`DEFAULT_PIPELINE`) with citations:

- `io.py` — LUT loader; `.sav`-only day reader (raises if the ES≡LT bug ever
  reaches the `.sav`s); BODC FCM reader; Jordan `uway_sst` reader. The BODC
  parameter-code mapping was transcribed from the deposit's own metadata
  document and it vindicates the design rule against guessing: **P700A90Z is
  *Synechococcus* and P701A90Z is *Prochlorococcus*** — the opposite of what
  the code strings suggest.
- `algorithm.py` — the retrieval: interpolate → standardise (N−1) → project →
  three regressions; raw floats + QC bitmask (Q&A #11), `pc_mapping=
  'operational'|'atbd'` (Q&A #9), `nasa_compat` clamp/truncate mode for the
  bit-exactness test, and the reconstruction-residual out-of-domain score
  (free, report §9.3).
- `pipeline.py` — stages 1–5: geometry screens, 1-min minimum-NIR selection,
  per-spectrum L1 glint fit (bounded 1-D search: for fixed ρ the optimal L1
  offset is a median), Lange QC, resampling, and the FCM matchup with the
  confirmed ≤10 m/shallowest + ±15 min median rules (plus Lange-strict
  `window=0`).
- `train.py` — `prcomp`-equivalent PCA (centring configurable, design §4.4),
  Lange's sd cut and backward-stepwise-AIC selection, per-taxon `use_sst`
  (design §4.3), and the sign/order-invariant basis comparison.
- `validation.py` — documented stub for prompt 14; `__init__.py` re-exports.

**Smoke-tested end-to-end on the real data**, not just synthetics:

- Invariants exact: standardised norm √123 to 1e-6, LUT orthonormal to 2e-8.
- Day 2014-280 (19°N): 24,364 raw → 20,495 geometry → 546 one-per-minute →
  312 QC'd spectra; retrieved Pro ≈ 3.1×10⁵ cells mL⁻¹ (right magnitude for
  the tropical Atlantic), ρ_sky median 0.047.
- Full cruise in 16 s: 993,417 → 513,519 → 15,363 → 6,423 QC'd 1-min spectra;
  **30 matched stations** of 68 (CTD casts are often pre-dawn — no daylight
  radiometry — on top of the known CTD-only handicap, Q&A #35).
- SST cross-check: Jordan `uway_sst` vs HSAS ancillary, median |Δ| = 0.01 °C.
- **The headline:** retraining the PCA on those 30 matchups (no centring)
  recovers NASA's basis in order — **PC1 |cos| = 0.999, PC2 0.998**, PC3–5
  0.96–0.98 — from a different processing chain and a third of the sample
  count. The methodology reproduces; and early evidence already points to
  `center=False` as Lange's `prcomp` convention (design §4.4).

**Two bugs found by the smoke and fixed:** (1) a saturated regression when the
candidate PC pool approached the sample count (all-NaN p-values) — now a
documented dof cap plus a stepwise guard; (2) a pandas ≥2 resolution trap where
FCM timestamps parse as *microsecond* datetimes while the netCDF's are
nanoseconds, so a bare `astype('int64')` mixed units and silently voided every
SST matchup — now a `_time_ns()` helper normalises all time axes.

**`pytest -q`: 264 passed** (full suite, ~2 min; no MOANA-specific tests yet —
that is prompt 13). No new Q&A questions: nothing arose that the design doc's
§9 does not already track. Nothing committed. **Files changed:**
`ioptics/moana/{__init__,io,algorithm,pipeline,train,validation}.py` (new),
`claude_prompts/moana_prompts.md` (this entry).

### 2026-08-16 (Prompt 11 — design doc rev. 2; PML follow-up list created)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Answers #33–35 reviewed.** All three confirm the proposals: FCM matchup uses
depth ≤ 10 m with the **shallowest** bottle per station (#33); the Jordan et al.
`uway_sst` is the primary SST with the HSAS ancillary column as cross-check
(#34); the underway FCM goes on the PML follow-up list (#35).

**Design doc updated** (`docs/design/moana_design.md`, rev. 2): §5.6's two
"pending" markers replaced with the confirmed rules (the values stay
parameterised, but these are now the defaults used for every reported result,
and the SST cross-check difference is stored as a QC diagnostic); the §9 open-
items table closes item 1 and points items 2 and 5 at the new follow-up file.
With #33–35 settled, **no open decisions remain on the matchup definition** —
prompt 12 (Code) is unblocked.

**Created `requests/PML_follow_up.md`**, consolidating every PML-bound item
scattered through Q&A #25, #27–31 and #35 into one email-ready list, priority-
ordered: (1) the underway FCM samples (→ Tarran; blocks full-fidelity
retraining — the CTD-only handicap is quantified from Lange's own Table 2);
(2) the ES≡LT CSV export bug; (3) the missing DOY 273/298 inside the training
window; (4) a processing-chain note for the Level-2 delivery, with our
empirical findings listed for confirmation; (5) native-grid spectra;
(6) Lange's processed `prcomp` input matrix if PML holds it; (7) the suggestion
to deposit the AMT24 radiometry with BODC, citing the AMT23/25/26/28 precedent.
The file explicitly separates PML items from the NASA-bound questions (those
stay in the report's §11).

**No `pytest` run:** documentation only, no package code. Nothing committed.
**Files changed:** `docs/design/moana_design.md` (rev. 2),
`requests/PML_follow_up.md` (new), `claude_prompts/moana_prompts.md` (this
entry).

### 2026-08-16 (Prompt 10 — design doc written; two new discrepancies found; papers/ restored)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Deliverable:** `docs/design/moana_design.md` — the retrieval spec (§3, with the
three free invariants any implementation must satisfy), the discrepancy register
(§4), the AMT24 Level-2 → training-matrix pipeline with every Lange threshold
pinned to the paper text (§5), the retraining and basis-comparison plan (§6), the
`ioptics/moana/` layout (§7), the test strategy (§8), and open items (§9).
Q&A #33–34 remain unanswered; their proposals are carried as configurable
defaults, not blockers.

**Side-task first: the Context PDFs were gone.** `papers/lange2020.pdf` and
`papers/moana_atbd.pdf` — both cited by the Context section and needed for the
pipeline thresholds — did not exist on disk (`papers/*.pdf` is gitignored, so the
machine migration silently dropped them). Restored: the ATBD from
oceancolor.gsfc.nasa.gov, and Lange+2020 as the **NTRS accepted manuscript**
(citation 20205006419; the Optica publisher PDF is bot-gated). Same content and
equation numbers, different pagination — noted in the design doc.

**With the paper back, the Lange chain is now pinned verbatim** (design §5.2–5.4):
tilt < 5°; 10° < θ₀ < 80°; 50° < Δϕ < 170°; 1-minute minimum-`Lt(750–800)`
selection; per-spectrum L1 fit of (ρ_sky, L_NIR) over 750–800 nm; local time
09:00–17:00; no negatives in 400–700 nm; |d²Rrs/dλ²| ≤ 2×10⁻⁴ in 610–660 nm.
Also confirmed: matchups were by exact date-hour-minute (so the design carries a
`window=0` Lange-strict mode alongside our ±15 min median), the PCA was `prcomp`
with cuts at 0.1 % of PC1's sd (→20 PCs) then p > 0.05 (→14), and selection was
backward stepwise by AIC.

**Two new findings, both now in the design doc:**

1. **The paper gives picoeukaryotes an SST term; the ATBD and the operational
   LUT do not** (Lange Eq. 7 *and* its prose: "SST was used as an additional
   predictor for Prochlorococcus and picoeukaryotes"). Unresolvable from
   documents since the paper tabulates no coefficients. `train.py` gets a
   per-taxon `use_sst` switch; the item belongs on the report's §11 NASA list at
   its next revision (not edited into the report now — not asked).
2. **The BODC flow-cytometry deposit is CTD-bottle-only** (814 rows, 68 stations,
   130 samples ≤ 10 m) — Lange's underway FCM samples are absent, which is
   exactly the configuration Lange showed degrades *Synechococcus*. Posted as
   Q&A #35 with the PML/Tarran follow-up ask.

Also recorded (§4.4): the `prcomp` **centring ambiguity** — R's default column-
centring is inconsistent with how the operational code applies the loadings, so
`train.py` supports both conventions and the basis comparison decides empirically.

**No `pytest` run:** design only, no package code. Nothing committed. **Files
changed:** `docs/design/moana_design.md` (new), `papers/lange2020.pdf` +
`papers/moana_atbd.pdf` (restored, gitignored), `claude_prompts/moana_prompts.md`
(Q&A #35, this entry).

### 2026-08-15 (Prompt 9 — radiometry answers locked; prompts 10–13 revised; two new questions)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Answers #27–32 reviewed; decisions now locked:**

| # | Decision | Consequence |
|---|---|---|
| 27 | Provider is PML (Brewin+); no processing note yet — "do your best for now" | The empirically-determined processing level (fixed ρ = 0.0280, nothing else) stands as our working truth until PML says otherwise |
| 28 | `.sav` files are the source of truth; CSVs ignored | JXP reports the ES≡LT export bug to PML |
| 29 | **We implement the Lange chain ourselves** from LT/LI/ES | The shipped fixed-ρ RRS is kept only as a cross-check; pipeline scope grows accordingly |
| 30 | Missing DOY 273 & 298 go on the PML follow-up | Matchups on those days are lost for now |
| 31 | Accept the double interpolation (native → 3.5 nm → 2 nm) | Revisit with PML later; no native-grid ask yet |
| 32 | ±15 min median binning onto FCM samples, spread kept as per-band uncertainty | Defines the training-matrix construction |

**Two new questions posted** (Q&A #33–34): the FCM depth cut for the matchup
(propose ≤10 m, shallowest sample) and which SST feeds the *Prochlorococcus*
regression (propose Jordan et al. `uway_sst` as primary, HSAS ancillary as
cross-check). Both shape the training matrix itself, so they warrant your call
rather than a buried assumption.

**Prompts 10–13 revised**, as instructed. Two kinds of change: (a) folded the
#27–32 decisions into the design/code/test/validation scopes — the pipeline
(.sav reader, Lange chain, 2 nm resampling, ±15 min binning) is now explicit in
prompt 10, and prompt 13 names the three locked validation targets including the
Brewin 2023 AMT23/25/28 data; (b) fixed the file paths, which still pointed at
bare `code/`, `tests/`, `validation/` directories that don't exist — they now
match the locked prompt-3 layout (`ioptics/moana/`, `ioptics/tests/test_moana.py`,
LUTs from `ioptics/data/moana/`). The original paths are preserved in a note
above prompt 10.

**One self-correction during this pass:** I initially wrote placeholder ">A."
answer lines under my own new questions #33–34 — those lines are yours to write,
not mine — and removed them in the same pass. Flagging it here since the diff
touched the Q&A answer format.

**No `pytest` run:** no package code touched. Nothing committed. **Files
changed:** `claude_prompts/moana_prompts.md` only — Q&A #33–34 added, prompts
10–13 revised, this entry.

### 2026-08-15 (Prompt 8 — AMT24 radiometry inspected; upstream L2 with a broken ES export)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**What I did.** The prompt's `Wadiometry` folder is `Radiometry` on disk (typo
assumed). Enumerated the tree (37 DOY dirs, 8.7 GB), parsed the CSVs with pandas
and the `.sav` files with `scipy.io.readsav` (ocean14 env), checksummed every
day's files, ran a cruise-wide sweep of row counts / positions / ancillaries
(993,417 spectra, 52.2°N → 50.5°S), and determined the processing level
empirically by inverting the glint equation: the shipped RRS is exactly
`(LT − 0.0280·LI)/ES` — fixed Mobley ρ, no NIR-residual subtraction, no
geometry screening, no binning. Full findings in the prompt-8 report; six
questions posted as Q&A #27–32 (provenance, the ES bug, re-deriving Rrs, the
two missing days in the Lange window, double-interpolation, matchup binning).

**Two of my own wrong turns, both caught in-pass.** (1) My first "glint" estimate
(ρ ≈ 0.06, spectrally variable) was an artifact of the ES bug itself — I was
unknowingly computing with LT in place of ES; only after `md5sum` exposed
ES ≡ LT and I switched to the `.sav`'s true `matrix_es` did the clean fixed
ρ = 0.0280 appear. (2) I briefly "found" garbage in the last three spectral
columns of the `.sav` matrices — that was my own column bookkeeping (I assumed
all 17 ancillary columns lead; actually 3 of them trail the wavelengths). The
CSV header order is faithful; no such defect exists. Recorded here because the
same trap awaits the eventual reader code.

**The one genuine defect** — `*_ES.dat` byte-identical to `*_LT.dat`, verified by
checksum on all 37 days — means the reader must use the `.sav` files.

**No `pytest` run:** no package code touched; this pass was data inspection only.
Nothing committed. **Files changed:** `claude_prompts/moana_prompts.md` — Q&A
Radiometry section filled (#27–32), new prompt-8 report, this entry.

### 2026-08-01 (Prompt 7 — re-searched amt-uk.org; found published AMT hyperspectral Rrs for AMT23/25/26/28)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**The honest answer to "how carefully did you search that page": I didn't.** In
prompt 6 I put `https://amt-uk.org/data/` in a list of sources "checked and found
nothing". What I actually ran against it was a `curl -o /dev/null -w "%{http_code}"`
probe in a batch of five URLs, got 200, and never looked at the content. Listing it as
searched was a misrepresentation of my own work, and the challenge was warranted. I
have annotated the prompt-6 report to say so rather than quietly fixing it.

**Reading it properly.** The page is a signpost: "All data sets are avaliable from the
British Oceanographic Data Center (BODC) website. The remaining AMT data sets are
available upon request to BODC, with future web delivery under development." Its
"AMT data policy" link is dead (404), as is the parent BODC AMT project directory.
Because the site is a JavaScript-rendered WordPress SPA — plain fetches return
"Loading..." — I also queried `wp-json/wp/v2/pages?per_page=100` to be sure nothing
was hidden behind the renderer: 12 pages total, nothing else data-bearing. The AMT24
cruise page has personnel and cruise reports and no data links. So amt-uk.org hosts
no data at all; it defers entirely to BODC.

**Which is where the actual find came from.** Taking the page's instruction seriously
meant searching BODC for AMT *optical* holdings instead of AMT24 holdings — and every
query I had run in prompt 6 was scoped to the string "AMT24". Rescoping surfaced
**10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6** (BODC, 2023; Brewin, Pitarch
Portero, Dall'Olmo, van der Woerd, Lin, Sun, Tilstone): Secchi depth, Forel-Ule
colour, Chl-a, **hyperspectral remote-sensing reflectance**, and diffuse/beam
attenuation at 127 stations on **AMT23, AMT25, AMT26 and AMT28**.

That is the search-design lesson worth recording: my prompt-6 searches were thorough
*within* a scope that was too narrow, and the negative result they produced was
therefore stronger-sounding than it deserved. "Not published for AMT24" is what I had
evidence for; I let it read as "not published for AMT".

**Why it matters.** AMT24 is not in the dataset — the "2013 to 2018" range invites the
opposite assumption, so I checked the abstract's explicit cruise list rather than
inferring. But Lange et al.'s five held-out cruises are AMT20, 22, **23**, **25**,
**28**, so this covers **three of the five**, and we already hold the matching
flow-cytometry DOIs for all three. Validation target (ii) is therefore unblocked — and
can be done *better* than the original paper managed it, since Lange et al. tested
held-out skill by pushing the multispectral model through Aqua-MODIS reflectance
(Pro MAE 2.26), whereas we can put in-situ **hyperspectral** `Rrs` through the
hyperspectral coefficients MOANA actually ships. No published number exists for that
combination.

**Correction issued.** Prompt 6's §4 implied AMT radiometry was intended for cal/val
but "never deposited". True of AMT24, false of AMT generally. I added a superseded
banner to the prompt-6 report rather than editing the claim away, since the reasoning
trail is the useful part. The narrow conclusion — AMT24's `Rrs` is unpublished —
survives both the DataCite enumeration and the SeaBASS enumeration unchanged. And the
new dataset strengthens the case for asking: same data type, same programme, same
people, same archive that has already published it four times.

**Still browser-only.** The new dataset has the same BODC PDL limitation (base,
`data/` and `download/` all return the SPA shell), so it goes on the Q&A #24 list —
now as the highest-priority item of the three, ahead of the AMT24 cell counts, because
it unblocks a whole validation target rather than one cruise.

**No `pytest` run:** no code touched. Nothing downloaded this pass (everything found
is browser-gated). Nothing committed.

**Files changed.** `claude_prompts/moana_prompts.md` — new prompt-7 report, a
superseded banner on the prompt-6 report, the new DOI added to Context, Q&A #24
extended, and this entry.

### 2026-08-01 (Prompt 6 — AMT24 data: 576 MB obtained; training Rrs shown to be unpublished)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

The searched-sources list and evidence are in the prompt-6 report above; this is the
work record.

**Obtained (`$OS_COLOR/AMT24/`, 576 MB total).** The Jordan et al. 2025 underway
dataset from Zenodo (571 MB, CC-BY-4.0) and the 170-page AMT24 cruise report from PML
(4.7 MB). Both open, no credentials needed.

I opened the netCDF rather than trusting its title, which was worth doing: it is
richer than "surface IOPs and pigments" suggests — 110 variables on a 1-minute
underway grid, 32,174 records covering the full transect (−48.3° to +47.8°,
26 Sep – 31 Oct 2014), with 176-wavelength AC-S/AC-S2 particulate `ap`/`bp`/`cp`
*each carrying an uncertainty estimate*, AC-9 and BB3 channels, HPLC pigments, and a
full set of underway ancillaries. The find that matters for MOANA: **`uway_sst`
(8.22–29.27 °C)** — in-situ SST at underway resolution, which supplies the
*Prochlorococcus* model's `log₁₀(SST)` term directly and removes any need for a
satellite SST field on the AMT side. Verified absent, as expected: no `Rrs`, no
water-leaving or sky radiance (`uway_tir*`/`uway_par` are downwelling meteorological
sensors).

**The negative result, established three independent ways.** I did not want to report
"could not find it", so I went after proof:

1. **Enumerated the entire SeaBASS `PML/AMT` tree** — all nine cruises,
   AMT19 through AMT29. Every published file is HPLC pigments or inline AC-S/AC-9
   particulate IOPs, all dated `v2024xxxx` (the Jordan deposit). PML has never put AMT
   radiometry in SeaBASS, on any cruise.
2. **Queried DataCite for "AMT24"** — exactly 7 DOI-published datasets exist
   worldwide: cruise report, flow-cytometry counts, nutrients, CDOM, CTD profiles,
   N₂O, POC. None is radiometry. That is the complete published record.
3. **Read the BODC cruise inventory.** It explicitly records *"Optical measurements
   of surface water (remote-sensing reflectance), once per minute"* and the
   underway `Lt`/`Lsky`/`Ed` measurements — so the data exist — but of the 167
   published BODC series linked from that page, **all are CTD casts**; I keyword-checked
   every series row for radiance/irradiance/reflectance/optics and got zero hits. The
   inventory documents what was measured at sea, not what was published.

So the conclusion is "measured, documented, never deposited" rather than "I couldn't
find it", which is a materially different thing to tell JXP.

**Who holds it.** The cruise report was the payoff here: a Satlantic HyperSAS with
136-channel HyperOCR sensors, run by Dall'Olmo, Ross, Rasse Boada and Smyth at PML,
with the report stating the reflectance was intended for satellite cal/val. And
Priscila Lange — MOANA's first author — was aboard AMT24, so there are two routes to
the same file.

**A caveat I flagged rather than glossed.** Even with "the AMT24 Rrs" in hand we would
not reproduce Lange's PCA unless we get *their processed product*: they applied their
own glint minimisation over 750–800 nm, tilt and relative-azimuth screening, a
second-derivative noise filter, and interpolation to 2 nm. The right ask is the
414–660 nm / 2 nm matrix that actually went into `prcomp`.

**What I could not do.** BODC's Published Data Library is entirely
JavaScript-rendered and exposes no machine-reachable download route — I tried `data/`,
`download/`, `files/`, `?download=1`, and checked the DataCite record for a
`contentUrl`. So the flow-cytometry cell counts and CTD profiles need a browser;
raised as Q&A #24. SeaBASS file downloads need an account, which I have and did not
attempt to obtain, and I did not try to authenticate anywhere.

**Consequence for sequencing.** Validation targets (i) and (ii) are now formally
blocked on an external data request. Target (iii) and the bit-exactness check are
unaffected, which reinforces the existing order: implement MOANA, verify against a
PACE granule, chase AMT in parallel.

**No `pytest` run:** no package code touched. Nothing committed.

**Files changed.** `claude_prompts/moana_prompts.md` (prompt-6 report, Q&A #24–26,
this entry). New outside the repo: `$OS_COLOR/AMT24/amt24_final_with_debiased_chl.nc`
and `$OS_COLOR/AMT24/AMT24_cruise_report.pdf`.

### 2026-08-01 (Prompt 5 — report revised: NASA contact list, two new figures, and a new result about *where* MOANA fails)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

Read the four new answers and revised `reports/MOANA_Claude_Report.md` accordingly.

**#17 — stay descriptive.** §9 now says so explicitly: the report is the first
deliverable, none of the eight proposed improvements has been implemented or
benchmarked, and the next step is a Python implementation of MOANA *as it stands*.

**#18 — new §11, "Consolidated list for the NASA contact."** Split into questions
only the authors can answer (which PC assignment is correct; the unpublished
eigenvalues; `Cov(β̂)` so the product could carry an error bar; whether the AMT24
training `Rrs` can be shared; the 395–705 vs 414–660 contradiction; whether
HyperSAS→OCI instrument transfer was considered; whether an OCI validation is
planned) and defects worth reporting (the seven Appendix A items, ordered by how
likely they are to corrupt someone else's analysis — land-as-254 first). Written so
it can be lifted straight into an email. The §10 reminder JXP asked for now points
at it.

**#20 — two of three figures added.** Loaded the house data-viz guidance before
writing any chart code this time, having skipped that step in prompt 4. Followed the
documented palette: single-hue blue sequential ramp for abundance magnitude, the
first three categorical slots for taxon identity (that subset is documented as
validated all-pairs for colour-vision deficiency, which matters because these are map
and histogram forms where every pair is co-visible), the reserved `critical` red used
*only* for the clipped-zero defect and never as a series, and recessive greys for
land and chrome. Could not run the palette validator — there is no JS runtime on this
machine — so I relied on the documented validation of that exact subset rather than
eyeballing, and said so. Every non-data colour is named in a legend, so nothing is
carried by colour alone; the aqua slot carries direct value labels as its contrast
relief.

Rendered and inspected both figures, which caught two layout faults the code review
would not have: a horizontal colourbar overlapping the "longitude" axis label (fixed
by dropping the axis labels — degree ticks are self-explanatory on a map), and a
left-panel title so long it ran off the canvas (moved into the panel's empty
lower-right as an annotation). The third figure JXP approved — our retrieval
scattered against NASA's — genuinely cannot be made yet, since it *is* the
bit-exactness test and needs the implementation from prompt 8. I noted that in §10
rather than quietly dropping it.

**A new result, which is the most valuable thing in this revision.** Drawing the mask
map exposed something the summary statistics had hidden: the clipped *Prochlorococcus*
pixels are strongly spatially organised. Quantified on the granule — median
|latitude| 47.9° versus 31.7° for successful retrievals, 74.6 % poleward of 40°
versus 27.3 %, and 4.5× the picoeukaryote abundance (9,202 vs 2,035 cells mL⁻¹) where
the model fails. So the 18 % is **out-of-domain extrapolation into cold, productive,
picoeukaryote-dominated water** — exactly the regime absent from a training set with
Chl never above 1 mg m⁻³ — and not random failure. That is a far more useful
diagnosis than a bare failure rate, and it makes the §9.3 out-of-domain detector
empirically testable with no retraining and no new data.

I also wrote the fair counter-argument into §7.3 rather than only the damning
version: *Prochlorococcus* genuinely is rare poleward of 40°, so a negative
prediction there is the model straining to say "very little" with an estimator that
cannot express near-zero. The abundances are not wildly wrong. What the clamp
destroys is the *information* — absent, 10 cells mL⁻¹, and "the model broke" all
collapse to the same 0.

**Provenance of the figure script.** `reports/scripts/moana_report_figs.py` now has
the basis figure (data-free, runs anywhere) plus two product figures that skip with a
message if the granule is absent, mirroring the repo's Tier-1/Tier-2 convention. It
prints the ocean-only statistics as it goes, and those reproduced the report's
numbers exactly (18.0 / 8.9 / 0.0 % zeros; 0.02 / 0.02 / 0.51 % over `valid_max`) —
an independent regeneration rather than a restatement, since the script recomputes
from the granule.

**No `pytest` run:** no package code touched. Nothing committed.

**Files changed.** Modified: `reports/MOANA_Claude_Report.md` (new §11, §7.3 spatial
result, §9.3 empirical support, two figures embedded, scope and status notes),
`reports/scripts/moana_report_figs.py` (two new figures + palette constants),
`claude_prompts/moana_prompts.md` (Q&A #21–23 and this entry). New:
`reports/figures/moana_product_masks.png`, `reports/figures/moana_clipping.png`.

### 2026-08-01 (Prompt 4 — MOANA report written; earthaccess added to requirements; loadings figure generated)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**Read the new answers first, as instructed.** #13 → clipped values treated as
censored data *and* the unphysical fraction reported as a headline number; #14 →
product defects go in an appendix, and JXP wants reminding to raise them with NASA
(that reminder is now written into the report body so it cannot get lost); #15 →
document the Bailey & Werdell departure and use daily L4M, revisit true L2 later;
#16 → `earthaccess` added to `requirements.txt` with a comment recording the
`~/.netrc` test-skip convention. That last one was the only code-adjacent action and
is done.

**Wrote `reports/MOANA_Claude_Report.md`.** Structured so it is self-contained — the
equations, all 25 coefficients with their PC assignments, the file formats, and the
data-source collections are all in it, so an implementation can be built from the
report alone without re-reading the paper, the ATBD or the C. Ten sections plus two
appendices: what MOANA is; why the approach can work; the five algorithm steps; the
coefficient table; a provenance matrix of which constant comes from which source;
published performance; inconsistencies; an assessment; eight proposed improvements;
and the plan. Product-level defects went to Appendix A per #14, file formats to
Appendix B.

Two things I was careful about. First, **separating verified statements from
inferences** — behaviour read out of the C and the data files is verified, whereas
"the ATBD is probably right and the operational LUT probably has two coefficients
misplaced" is an inference from a single coincidence (the 14-PC union count) and is
labelled as such. Second, **not over-claiming on mechanism**: the PC5 loading has
narrow structure near 450/510/605 nm which would be a tidy phycobiliprotein story
for why hyperspectral helps *Synechococcus* most, and I wrote it down explicitly as
untested speculation rather than as a finding, because post-hoc pigment attribution
of PCA components is easy to get wrong.

**The improvements section is where I put the most thought.** Ranked by value per
unit effort, the substantive ones are: replace the identity-link linear
*Prochlorococcus* model with a **log-link GLM** so negatives become structurally
impossible (this is a better answer than either option Lange et al. considered,
since fitting `log₁₀(y)` by least squares carries a retransformation bias that is
plausibly what "significantly reduced the performance" actually measured); **feed
back the two amplitude numbers that standardisation throws away**
(`log₁₀(mean(Rrs))` and `log₁₀(sd(Rrs))`) as extra predictors, since they are already
computed and cost nothing; add **per-pixel uncertainty**, which is nearly free for a
linear model on fixed components; and add an **out-of-domain flag** from the
truncated-basis reconstruction residual `‖Rrs' − VVᵀRrs'‖`, which needs no
retraining and replaces silent extrapolation with a principled novelty score. Also
argued for retraining the basis at OCI's real spectral resolution (the current basis
inherits the 2014 HyperSAS's ~10 nm resolution, so it cannot contain the fine
structure PACE exists to measure), for regularisation instead of stepwise AIC at
n = 73, for **spatially blocked** cross-validation instead of random splits of
autocorrelated underway data, and for shipping an Rrs-only *Prochlorococcus* variant
so that studies of temperature–*Prochlorococcus* coupling are not circular.

**Made one figure.** `reports/scripts/moana_report_figs.py` draws the first six PCA
loadings from the vendored LUT, annotated with which taxa use each component. It
follows the repo's existing `reports/{figures,scripts}` convention and needs no
satellite or in-situ data, so it runs anywhere the repo is checked out. Useful side
effect: the script re-derives the PC assignments independently of my earlier
analysis and got the same answer (Pro 1/2/6/7, Syn 1/2/3/5/8/9/10/11/12/16, picoeuk
2/3/4/5/10/13/15), which is a genuine check on the decode rather than a restatement
of it. It also prints the invariants the report leans on — orthonormality to
1.2×10⁻⁷ and `‖Rrs'‖₂ = 11.091`.

**Corrections I caught while writing.** I had mixed up the multispectral and
hyperspectral numbers in Lange et al.'s CTD-versus-underway retraining experiment
(the −29 % bias is the multispectral case; hyperspectral is −16 % with MAE 1.37) and
fixed it against the paper. I also had to fix a matplotlib title that used `$…$` and
was silently mathtext-mangling `pca_picophyto.h5` into a subscript.

**Deliberately not done:** I did not run the bit-exactness test, even though it
would have let the report settle §7.1 definitively. It needs a 173 MB granule
download and interpolation code, which are prompts 5 and 7; the report specifies the
experiment precisely instead and flags it as the highest-value next step.

**No `pytest` run:** no package code was touched (the figure script is a report
artefact under `reports/scripts/`, matching the existing precedent). Nothing
committed.

**Files changed.** New: `reports/MOANA_Claude_Report.md`,
`reports/scripts/moana_report_figs.py`, `reports/figures/moana_pc_loadings.png`.
Modified: `requirements.txt` (added `earthaccess`),
`claude_prompts/moana_prompts.md` (Q&A #17–20 and this entry).

### 2026-08-01 (Prompt 3 — answers reviewed and acted on; LUTs vendored; PACE access verified; product defects found)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

Reviewed all twelve answers and acted on everything actionable without stepping into
prompts 4–9 (report, data download, design, code). The decision table and the full
findings are in the prompt-3 report above; this is the work record.

**Q12 — vendored the reference tables.** Created `ioptics/data/moana/` with
`pca_picophyto.h5` and `picophyt.json` copied from
`$OS_COLOR/ocssw/share/common/`, checksums re-verified after the copy. Wrote a
`README.md` there covering provenance (upstream tag, direct URLs, SHA-256s), the
dataset layout, the coefficient-slot→PC mapping, the known ATBD/LUT discrepancy, and
the algorithm recipe. The repo can now run MOANA with no OCSSW install.

**Q10 — verified access and de-risked the bit-exactness test.** `earthaccess` 0.18.0
authenticates from `~/.netrc`. Established that NASA's MOANA output is
`PACE_OCI_L4M_MOANA` v3.2 (there is **no** L3M MOANA collection in CMR, contrary to
the catalog page I cited in prompt 1 — Context corrected), that the Rrs input is
`PACE_OCI_L3M_AOP` v3.2 with 172 OCI bands spanning 346–719 nm (so 414–660 is fully
covered and MOANA never extrapolates on real data), and that the two grids are
exactly co-registered — all 1400 MOANA latitudes and 1100 longitudes appear verbatim
in the global AOP grid, so alignment is by coordinate value with no regridding.

Also closed the half of Q&A #5 that went unanswered: GHRSST CMC L4 SST is
`CMC0.1deg-CMC-L4-GLOB-v3.0`, reachable through `earthaccess` at 0.1° — the same
resolution as the MOANA grid — so no `ocpy` loader is needed.

The useful realisation here: *Synechococcus* and picoeukaryotes need no SST, so the
bit-exactness test can run from Rrs alone and by itself settles the disputed
*Synechococcus* PC index (operational PC16 vs ATBD U13). We are not blocked on SST
for the question that matters.

**Four defects in the shipping product.** I downloaded one daily 0.1° MOANA granule
(2 MB) rather than trusting the catalog, and it was worth doing. The most important
find is that **land is encoded as `254 cells mL⁻¹`** — not as `_FillValue`, and
inside the declared valid range, so the honest screen `valid_min ≤ x ≤ valid_max`
swallows 450,592 land pixels as real abundance. This nearly fooled me: my first pass
reported a median of 254 for all three variables and I only caught it because an
identical median across three independent regressions is impossible. Spot-testing
the Congo basin, Sahara, Amazon, Andes and interior Greenland (all 254) against four
open-ocean points settled it. Also found six `INT32_MIN` pixels (undefined behaviour
from the `(int32_t)` cast), `valid_max` unenforced (picoeukaryotes reach 68× their
declared ceiling), and — most consequential — **`prococcus_moana` is exactly zero in
18.0 % of real ocean retrievals**, which is the silent negative-clamp firing on
nearly a fifth of the field, plus 8.9 % zeros in *Synechococcus* from integer
truncation. Correct ocean-only statistics are tabulated in the report.

That last one vindicates answer #11(b). Had we mirrored NASA's behaviour we would
have log-transformed or silently discarded a fifth of the Prochlorococcus field
without noticing.

**Method note for my own future reference:** I twice caught myself about to record a
number that was an artefact — the "median 254" above, and a spurious 0.1° lat offset
that came from my own `searchsorted` misuse on descending coordinate arrays rather
than from the data. Both were resolved by re-checking with a different method
(raw unscaled reads; exact set membership) before writing anything down.

**No `pytest` run:** this step added a data directory and a README, no package code.
Nothing was committed.

**Files changed.** `claude_prompts/moana_prompts.md` — corrected the product/variable
names in Context, added the prompt-3 report, Q&A #13–16, and this log entry. New:
`ioptics/data/moana/{pca_picophyto.h5,picophyt.json,README.md}`. New outside the
repo: `$OS_COLOR/PACE/PACE_OCI.20250701.L4m.DAY.MOANA.V3_2.0p1deg.nc`.

### 2026-08-01 (Prompt 2 — OCSSW installed; MOANA loadings and coefficients recovered and decoded)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

The full findings are in the prompt-2 report above; this is the work record.

**Searching.** Confirmed there is exactly one public MOANA implementation, the C
routine inside OCSSW. `github.com/nasa/ocssw` does not exist, and the Ocean Color
doxygen browser that hosted the readable source in prompt 1 has been retired
(everything under `/docs/ocssw/` now 301s to the OB.DAAC landing page). No Python
or R port exists — not in `nasa/oceandata-notebooks`, not in the PACE Hackweek
material — and Lange's original R analysis code was never released. So OCSSW is the
only authority, which is why installing it mattered.

**Installing.** `install_ocssw` and `manifest.py` come from
`oceandata.sci.gsfc.nasa.gov/manifest/` (the `/ocssw/install_ocssw` path 404s).
Listed 30 tags and took the newest, **V2026.4**, installing the `--src` and
`--common` bundles into `$OS_COLOR/ocssw` (~11.5 GB; `share/common` alone is 10 GB).
The first run died partway through `common` on a read timeout from NASA; re-running
with `--common` resumed cleanly and finished. Disk left: ~473 GB.

Worth recording for anyone repeating this: **the 11.5 GB install was not actually
necessary.** The bundle manifest at `…/manifest/tags/V2026.4/common/manifest.json`
lists every file with its size, sha256 and originating tag, and individual files are
directly fetchable. Both MOANA tables — 24,488 bytes and 390 bytes — came down in
under a second from `…/manifest/tags/T2023.31/common/`, and I verified both against
the manifest checksums before trusting them. (The `T2023.31` tag also tells us these
tables have not changed since 2023 and therefore predate the Dec-2024 ATBD.)

**What the files contain.** `pca_picophyto.h5` holds `component`, a (124, 45)
float32 loading matrix, and `wavelength`, int32 414–660 nm in 2 nm steps. I checked
that the columns are genuinely orthonormal (`max|offdiag(VᵀV)| = 1.2e-7`), so it is
a proper truncated eigenvector basis. `picophyt.json` holds `npc = 17` and three
comma-separated coefficient strings of length 19/18/18 — exactly matching the C
allocations `npc+2`, `npc+1`, `npc+1`, which confirms the zero-padding hypothesis I
raised in Q&A #6.

The thing I most wanted to be sure of: **there is no training mean to recover.** The
standardization is per-spectrum across wavelength (Lange Eq. 3), not feature-wise
centering, and the scores are a plain dot product with the loadings. So these two
small files plus the code are genuinely everything needed to run MOANA exactly as
NASA runs it. Q&A #1 is closed.

**Reading the code.** `get_Cpicophyt.c` in the installed tree is 279 lines against
the 256 of the 2022 snapshot I read in prompt 1. The differences matter: the output
buffer changed from `float*` to `int32_t*` (so abundances are now **truncated** to
whole cells/mL), a valid-band filter was added ahead of the interpolation, and an
`npix` realloc was added so the algorithm can run under `l3gen` — consistent with
the operational products being L3/L4 mapped. Also catalogued the surrounding
plumbing (`l2prod.h` ids 353–355, the `prodgen.c` dispatch) and pulled the official
units and valid ranges from `share/common/product.xml`.

**Verification I did.** Rather than hand-count coefficient positions, I decoded
`picophyt.json` with a script replicating the C indexing exactly
(`pro_coef[i+2]·pc[i]`, `syn_coef[i+1]·pc[i]`, `apeuk_coef[i+1]·pc[i]`) and
diffed the result against the ATBD equations. Picoeukaryotes match perfectly;
Prochlorococcus and *Synechococcus* each have their final coefficient on a
different PC than the ATBD states, with identical values. Separately I
re-extracted the ATBD coefficient table in `pypdf` layout mode to rule out a
column-mixing artefact in my prompt-1 transcription — the numbers were correct.

I also ran a self-consistency check on the two files together. Because a
standardized spectrum has `‖Rrs'‖₂ = √123 = 11.09` and `V`'s columns are
orthonormal, every score is bounded by ±11.09; feeding that ceiling through the
Prochlorococcus regression gives 1.5–3.6×10⁵ cells/mL across SST 15–28 °C, which is
physically right and sits inside the official `validMax` of 6×10⁵. That the huge
~10⁶ coefficients only produce sane abundances when `U₁` is near its ceiling — which
is exactly what >96 % of covariance on PC1 implies — is good evidence the loading
file and coefficient file belong to each other and that my index mapping is right.

**The substantive finding.** The union of PCs used across the three taxa is 14
distinct PCs under the ATBD's assignment and 15 under the operational JSON's, and
Lange+2020 says exactly 14 survived the significance cut. That is a real argument
that `picophyt.json` has two coefficients in the wrong slots, i.e. the shipping PACE
MOANA products may not implement the published algorithm. I have deliberately *not*
asserted this as established — the alternative is that the ATBD subscripts are
typos — and have proposed settling it empirically by reproducing a real PACE granule
before contacting the authors (Q&A #9, #10).

**No `pytest` run:** this step wrote no package code. Nothing was committed.

**Files changed.** `claude_prompts/moana_prompts.md` — added the prompt-2 report,
four new Q&A items (#9–#12), and this log entry. New on disk outside the repo:
`$OS_COLOR/ocssw/` (OCSSW V2026.4).

### 2026-08-01 (Prompt 1 — context assembled; MOANA ATBD and OCSSW implementation located)

**Model:** run as Claude Fable 5, per "Use Fable if you can".

**What I read.** `docs/design/IOPtics_design.md` (v0.15) and
`docs/design/IOPtics_implementation.md` (v0.22) in full, and all 24 pages of
`papers/lange2020.pdf` (text extracted with `pypdf` in `ocean14`).

**Lange et al. 2020 in brief.** *Optics Express* 28(18), 25682–25705,
DOI `10.1364/OE.398127`. Empirical principal-component regression predicting
surface cell abundances of *Prochlorococcus*, *Synechococcus* and autotrophic
picoeukaryotes. Trained on AMT24 (Sep–Nov 2014): underway HyperSAS Rrs
interpolated to 2 nm over 414–660 nm, each spectrum standardized (Eq. 3) to
suppress amplitude and expose shape, PCA via R `prcomp`, then backward stepwise
multilinear regression on the PC scores with AIC-based selection, optionally
adding SST. 20 PCs survived the variance cut, 14 the significance cut. PC1 alone
carries >96 % of the covariance (backscatter slope + water absorption), PC2 ~3.5 %
(Chl absorption), PC3 ~0.16 % (accessory pigments). *Prochlorococcus* is fit
**linearly** (it is normally distributed); the other two in log10. SST helps Pro
markedly (hyperspectral MAE 1.49→1.31, R² 0.42→0.82) and does nothing for Syn.
Hyperspectral beats 7-band MODIS throughout — the paper's headline argument for
PACE. Skill degrades on satellite Rrs (Syn MAE 2.04) and degrades further on the
five held-out cruises (Pro MAE 2.26, bias 1.75). One methodological lesson worth
carrying: the fine 30-minute underway sampling across the South Atlantic front is
what gave *Synechococcus* its dynamic range — retraining on CTD casts alone drove
its bias from ~0 to −29 %.

**What I found on the NASA side.** Three things, in increasing order of value:

1. MOANA is a *current operational PACE product*, not just a paper. Level-3 and
   Level-4 Regional Mapped MOANA (v3.1, and an unregistered v3.2 + NRT), variables
   `prochlorococcus_moana` / `synechococcus_moana` / `picoeuk_moana` in cells mL⁻¹,
   2024-03-05 onward. It is part of the first PACE OCI PCC suite and was used for
   the PACE "first light" image. The ATBD concedes the *satellite* products **have
   not been validated**, which is a real opening for us.
2. **The ATBD (v1.2, 16 Dec 2024)** — downloaded to `papers/moana_atbd.pdf`,
   DOI `10.5067/0AV5267G5C77`. It confirms MOANA = "Multiple Ordination ANAlysis"
   and, crucially, **tabulates all 25 regression coefficients** that Lange+2020
   omits, along with which PC index each multiplies:
   - `Pro       = p1 + p2·log10(SST) + p3·U1 + p4·U2 + p5·U6 + p6·U17`  (linear)
   - `log10(Syn) = s1 + s2·U1 + s3·U2 + s4·U3 + s5·U5 + s6·U8 + s7·U9 + s8·U10 + s9·U11 + s10·U12 + s11·U13`
   - `log10(picoeuk) = a1 + a2·U2 + a3·U3 + a4·U4 + a5·U5 + a6·U10 + a7·U13 + a8·U15`
   with `U_i = Σ_λ Rrs'(λ)·V(λ,i)` over 414–660 nm. SST is the PACE reference
   field (GHRSST CMC L4). Note the hyperspectral picoeukaryote model needs **no**
   SST, consistent with Table 1 of the paper.
3. **The OCSSW implementation**, `src/l2gen/get_Cpicophyt.c` (Minwei Zhang,
   2023-09-29). The Ocean Color doxygen tree has been retired — every
   `oceancolor.gsfc.nasa.gov/docs/ocssw/…` URL now 301s to the OB.DAAC landing
   page — so I recovered it from a 2025-07-12 Wayback snapshot. Reading it settles
   several things the prose leaves vague: Rrs is **linearly** interpolated onto the
   PCA grid; standardization uses the **sample** sd (N−1); scores are a plain dot
   product of the standardized spectrum with a loading column (no re-centering on a
   training mean, no singular-value scaling); Syn and picoeuk are exponentiated,
   Pro is not; and a negative Pro is clamped to 0 when the other two are positive,
   else flagged `BAD_FLT`. It also reveals where the two missing data files live:
   loadings in `$OCDATAROOT/common/pca_picophyto.h5` (`component` [nwave × npc],
   `wavelength`), coefficients in `$OCDATAROOT/common/picophyt.json`. I did not
   pursue the OCSSW download itself — that is prompt 2 — but flagged it there.

**AMT24 and its data.** The training cell counts have a citable BODC DOI, as do
all five Lange+2020 validation cruises (AMT20/22/23/25/28) — the full table is now
in Context, every DOI resolved against DataCite and `doi.org`. The most valuable
incidental find is **Jordan et al. 2025** (*ESSD* 17, 493–516): underway
particulate IOPs (~400–720 nm) and HPLC pigments for nine AMT cruises **including
AMT24**, on SeaBASS and Zenodo. That is in-situ IOP truth on the exact cruise that
trained MOANA, and it is a plausible new IOPtics dataset independent of this MOANA
work. Second most valuable: **Lin et al. 2022** (*Opt. Express* 30, 45648), an
uncertainty budget for continuous above-water AMT radiometry — the principled
`varRrs` for an AMT fit, in place of the flat-5 % fallback `ioptics.prep` currently
uses for in-situ data without per-band errors.

**Gap I could not close.** I found no public archive for the AMT24 hyperspectral
Rrs itself, and no published PCA loading matrix. Both are prerequisites for a true
reproduction; see Q&A #1 and #2.

**Verification.** Every DOI written into Context was checked: resolved through
`doi.org` where the publisher permits it, and via the Crossref or DataCite REST
APIs for the six that return 403 to automated clients (Wiley, MDPI, AGU) or that
are DataCite-only. Two DOIs advertised on Earthdata — `…/L4M/MOANA/3.2` and
`…/L3M/MOANA/3.2` — are **not yet registered** and do not resolve; Context cites
v3.1 and says so. No `pytest` run: this step touched no package code.

**Files changed.** `claude_prompts/moana_prompts.md` (Context, Q&A, this log);
added `papers/moana_atbd.pdf` (6.9 MB, downloaded from OB.DAAC).
