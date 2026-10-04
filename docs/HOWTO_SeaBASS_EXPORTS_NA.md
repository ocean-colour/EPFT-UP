# HOWTO: download the EXPORTS North Atlantic (2021) files from SeaBASS

**Purpose.** Execution #9 of `claude_prompts/2nd_derivative_prompts.md` (the
independent hold-out). We need the HPLC pigments for the 17 EXPORTS-NA
spectra that Kramer et al. (2024) added, and radiometry metadata that lets
each spectrum be tied to a time. SeaBASS refuses scripted access
(`cgi-bin/file_search.cgi` returns 403/444, and archive pages list files only
via JavaScript), so this is a manual, logged-in browser download.

> **Revision 2 (2026-10-03).** The first version told you to search for the
> file-name pattern `exports_na_*C-OPS_CAST*`. That returns nothing, because
> the File Search's **"Keyword Search Filters" box matches only affiliation,
> investigator, experiment or cruise names, not file names**. This version
> searches on the cruise name `EXPORTSNA` and narrows the results with the
> product, wavelength and "Experimental" settings. It also changes the
> radiometry to fetch: Kramer's 17 test spectra are 1 nm hyperspectral, and
> the hyperspectral EXPORTS-NA Rrs on SeaBASS is NASA GSFC's **HyperSAS** on
> RRS Discovery (DY131), stored per station. The UCSB C-OPS casts are
> multispectral profiler data and are now optional.

## What to download

Everything belongs to SeaBASS **cruise `EXPORTSNA`** (experiment EXPORTS,
May 2021, Porcupine Abyssal Plain; RRS James Cook = JC214, RRS Discovery =
DY131).

| Set | Files | Where (archive directory) | Needed? |
|---|---|---|---|
| **A. HPLC** | `EXPORTS_EXPORTSNA_JC214_process_rosette_HPLC_20230202_R1.sb`, `EXPORTS_EXPORTSNA_DY131_survey_rosette_HPLC_20230202_R1.sb`, `EXPORTS_EXPORTSNA_DY130_rosette_HPLC_20230202_R1.sb` | `UCSB/CRSEO/EXPORTS/EXPORTSNA/archive/` (Siegel group, i.e. Kramer's lab) | **required** (3 files) |
| **B. HyperSAS station Rrs** | the 34 files `EXPORTS_EXPORTSNA_DY131_HyperSAS_<yyyymmdd>_<hhmmss>_L2_Rrs_STATION_<n>_0_R1.sb` (13 stations, 2–29 May 2021) | `NASA_GSFC/EXPORTS/EXPORTSNA/archive/` | **required** (34 files; listed in the Appendix) |
| C. extra HPLC | `EXPORTS_EXPORTSNA_DY131_HPLC-inline_survey_R1.sb`, `EXPORTS_EXPORTSNA_DY131_HPLC-pump_survey_R1.sb` | `BOWDOIN/ROESLER/EXPORTS/EXPORTSNA/archive/` | optional |
| D. C-OPS casts | 99 `exports_na_{dy131,jc214}_<date>_<time>_C-OPS_CAST_R0.sb` | `UCSB/CRSEO/EXPORTS/EXPORTSNA/archive/` | optional (fallback for timing) |

Set B is small text, a few MB in total. If it's easier, you may take *all*
the DY131 HyperSAS **Rrs** files (about 180, non-station ones included)
instead of picking out the 34; I will ignore the rest. Do **not** fetch the
`…_L2_Es_…` files.

## Route 1: the archive browser (simplest; recommended)

You know the directories and file names, so the archive browser is the most
direct route. You need **three directories**.

1. Log in at <https://seabass.gsfc.nasa.gov/login/> with your NASA
   Earthdata account (the same one in `~/.netrc`).
2. **HPLC.** Open
   <https://seabass.gsfc.nasa.gov/archive/UCSB/CRSEO/EXPORTS/EXPORTSNA/archive/>.
   Wait for the table to fill (it is built by JavaScript), type `rosette_HPLC`
   in the table's **Search/filter box** (top right of the table), and
   download the 3 files by clicking each name. If a click opens the file as
   text, use the browser's *Save Page As…* or right-click → *Download Linked
   File*.
3. **HyperSAS station Rrs.** Open
   <https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/>,
   type `L2_Rrs_STATION` in the filter box (34 rows should remain), and
   download each one.
4. Optional: in
   <https://seabass.gsfc.nasa.gov/archive/BOWDOIN/ROESLER/EXPORTS/EXPORTSNA/archive/>,
   filter `HPLC`; for C-OPS, back in the UCSB/CRSEO directory, filter
   `C-OPS_CAST`.

The **"Download All"** button on a directory page makes one `.tgz` of the
whole directory. That's fine for `NASA_GSFC/…` if you'd rather not click
34 files (it also contains Es files, which I'll ignore), but avoid it for
`UCSB/CRSEO/…`, which also holds CTD, ADCP, nutrients, POC and more.

## Route 2: the File Search form (what you tried, with the right settings)

At <https://seabass.gsfc.nasa.gov/search/> → **File Search**, run two
searches.

**Search 1: HPLC.**
- *Measured between*: `2021-04-25` and `2021-06-05`.
- *Within the coordinates*: N `50`, S `48`, W `-16`, E `-13` (all
  EXPORTS-NA stations are near 49°N, 15°W).
- *Keyword Search Filters → Search String*: `EXPORTSNA` (a cruise name),
  with *Any*.
- *Data Use Warnings?*: **Yes** (include files carrying warnings).
- *Experimental?*: tick it, to be safe.
- *Products?*: choose "Find files containing any of the selected products"
  and tick **HPLC** under Grouped Products.
- Click **Search**. In the results, select the three
  `…_rosette_HPLC_20230202_R1.sb` files (and the two Bowdoin `HPLC-…`
  files if you like) and use the results page's download/bundle option.

**Search 2: hyperspectral Rrs.** The same dates, coordinates, keyword and
warnings settings, then:
- *Wavelength Options?*: **Hyperspectral**.
- *Products?*: "Find files containing any of the selected products", tick
  **AOP**.
- In the results, select the files whose names contain `L2_Rrs_STATION`
  (34), or all DY131 HyperSAS `…_L2_Rrs_…` files, and download them.

If a search still returns nothing, clear the coordinates (back to ±90/±180)
and the dates, keep only the keyword `EXPORTSNA`, and add the filters back
one at a time.

## Where to put the files

Unpack and place every file **unrenamed** under

    $OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA/

(`$OS_COLOR` = `/Users/xavier/Projects/Oceanography/data/Color/`).
Sub-folders are fine, e.g. `EXPORTS_NA/HPLC/` and `EXPORTS_NA/HyperSAS/`.
Note the download date if you can; I will record it with every file's
SHA-256.

Quick check:

    find $OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA -name '*rosette_HPLC*.sb' | wc -l        # expect 3
    find $OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA -name '*L2_Rrs_STATION*.sb' | wc -l      # expect 34
    head -40 "$(find $OS_COLOR/PANGAEA/Kramer2022/EXPORTS_NA -name '*JC214_process_rosette_HPLC*.sb')"

Each `.sb` file is plain text: a `/begin_header … /end_header` block, then
data columns named in `/fields=`. If a "file" turns out to be an HTML page
(it starts with `<!DOCTYPE html>`), the click saved the web page instead of
the data. Re-download it while logged in.

## What happens next (Claude, Execution #9)

1. Pin the files: a provenance record with SHA-256s and a small SeaBASS
   reader.
2. **Identify Kramer's 17 spectra:** compare each test spectrum (and its
   position) with the HyperSAS station Rrs. Each match gives the spectrum a
   date and time. Spectra that match more than one file, or none, mean I
   stop and ask. If HyperSAS turns out not to be the source, the C-OPS
   casts (set D) are the fallback for timing.
3. Match each timed spectrum to surface (≤ 7 m) rosette HPLC within ±2 h
   and nearby (Kramer 2022's rule), averaging replicates. Cross-check
   against the HPLC chlorophyll already in Kramer's test file, which must
   agree. Ambiguities again mean stop and ask.
4. Run the full hold-out with no retraining: the reproduced PCR, the Tchla
   nulls and the shared W with intervals, on 13 pigments and 12 ratios,
   reporting the processing caveat.

## Appendix: exact archive URLs (sets A and B)

Use these in a logged-in browser (a scripted `curl` is refused).

    https://seabass.gsfc.nasa.gov/archive/UCSB/CRSEO/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY130_rosette_HPLC_20230202_R1.sb
    https://seabass.gsfc.nasa.gov/archive/UCSB/CRSEO/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_survey_rosette_HPLC_20230202_R1.sb
    https://seabass.gsfc.nasa.gov/archive/UCSB/CRSEO/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_JC214_process_rosette_HPLC_20230202_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210502_142906_L2_Rrs_STATION_0_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210502_150235_L2_Rrs_STATION_0_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210502_160202_L2_Rrs_STATION_0_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210505_125609_L2_Rrs_STATION_1_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210505_133834_L2_Rrs_STATION_1_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210505_140133_L2_Rrs_STATION_1_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210513_123102_L2_Rrs_STATION_6_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210513_130233_L2_Rrs_STATION_6_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210513_140231_L2_Rrs_STATION_6_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210516_113528_L2_Rrs_STATION_8_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210516_120312_L2_Rrs_STATION_8_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210516_130238_L2_Rrs_STATION_8_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210517_115600_L2_Rrs_STATION_9_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210519_114959_L2_Rrs_STATION_10_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210521_143257_L2_Rrs_STATION_11_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210521_150229_L2_Rrs_STATION_11_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210521_160237_L2_Rrs_STATION_11_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210522_101943_L2_Rrs_STATION_12_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210522_110233_L2_Rrs_STATION_12_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210522_120240_L2_Rrs_STATION_12_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210524_125514_L2_Rrs_STATION_13_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210524_130227_L2_Rrs_STATION_13_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210524_140256_L2_Rrs_STATION_13_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210525_133319_L2_Rrs_STATION_14_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210525_141233_L2_Rrs_STATION_14_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210525_150734_L2_Rrs_STATION_14_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210527_104011_L2_Rrs_STATION_16_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210527_110245_L2_Rrs_STATION_16_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210527_120237_L2_Rrs_STATION_16_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210528_112323_L2_Rrs_STATION_17_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210528_120232_L2_Rrs_STATION_17_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210528_130235_L2_Rrs_STATION_17_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210529_113714_L2_Rrs_STATION_18_0_R1.sb
    https://seabass.gsfc.nasa.gov/archive/NASA_GSFC/EXPORTS/EXPORTSNA/archive/EXPORTS_EXPORTSNA_DY131_HyperSAS_20210529_120234_L2_Rrs_STATION_18_0_R1.sb
