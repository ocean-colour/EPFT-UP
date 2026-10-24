"""New AMT data #2: stage Tom Jordan's ``AMT_Rrs_HPLC.csv`` on the laptop and the shared drive.

Placement decided in the prompt doc (Q&A → New AMT data, Q1–Q2):

* laptop: ``$OS_COLOR/AMT/Rrs_HPLC/`` (copied, not moved, from ~/Downloads);
* Google shared drive *AI_for_Ocean*: ``data/Color/AMT/Rrs_HPLC/`` (mirrors
  ``$OS_COLOR``), with ``--drive``.

Writes next to the CSV:

* ``README.md``: source, checksum, columns, cruises, known issues, the
  processing defaults EPFT-UP adopted, and sharing terms;
* ``AMT_Rrs_HPLC_flags.csv``: one row per input row with issue flags
  (shared spectrum, non-positive Rrs, implausible salinity, missing MVchlb).

The CSV is never modified. An existing copy with a different checksum stops
the script (new versions from Tom get a dated suffix instead).

Run with::

    conda run -n ocean14 python scripts/sdp/stage_amt.py [--drive]
"""
import argparse
import datetime
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))
from epft_up.sdp import data as sdpdata  # noqa: E402

SRC = Path.home() / 'Downloads' / 'AMT_Rrs_HPLC.csv'
NAME = 'AMT_Rrs_HPLC.csv'
LOCAL = Path(os.environ['OS_COLOR']) / 'AMT' / 'Rrs_HPLC'
DRIVE = (Path.home() / 'Library' / 'CloudStorage' / 'GoogleDrive-jxp@ucsc.edu' /
         'Shared drives' / 'AI_for_Ocean' / 'data' / 'Color' / 'AMT' / 'Rrs_HPLC')
SAL_MIN_OPEN_OCEAN = 30.0
#: Files written by scripts/sdp/amt_harmonise.py, mirrored with --drive when present
DERIVED = ('AMT_Rrs_HPLC_harmonised.csv', 'AMT_Rrs_HPLC_harmonised.json',
           'AMT_Kramer_folds.csv', 'woa23_salinity_fix.csv', 'woa23_salinity_fix.json')


def flags(df):
    rcols = [c for c in df.columns if c.startswith('Rrs_')]
    R = df[rcols].to_numpy(float)
    key = pd.Series([hash(r.tobytes()) for r in np.round(R, 10)], index=df.index)
    gid = key.map({k: i for i, k in enumerate(key.unique())})
    n_in = gid.map(gid.value_counts())
    F = pd.DataFrame({
        'Sample_ID': df['Sample_ID'], 'Original_ID': df['Original_ID'],
        'Cruise': df['Cruise'],
        'spectrum_group': gid,
        'shared_spectrum': n_in > 1,
        'n_nonpositive_Rrs': (R <= 0).sum(axis=1),
        'nonpositive_Rrs': (R <= 0).any(axis=1),
        'bad_salinity': df['Sal'] < SAL_MIN_OPEN_OCEAN,
        'MVchlb_missing': df['MVchlb'].isna(),
    })
    F['any_flag'] = F[['shared_spectrum', 'nonpositive_Rrs', 'bad_salinity']].any(axis=1)
    return F


def readme(df, F, sha):
    t = pd.to_datetime(df['DateTime'], format='%m/%d/%Y %H:%M')
    g = df.assign(t=t).groupby('Cruise')
    rows = '\n'.join(
        f"| {c} | {len(d)} | {d['t'].min():%Y-%m-%d} – {d['t'].max():%Y-%m-%d} | "
        f"{'/'.join(sorted(set(d['HPLC_lab'])))} | {d['Lat'].min():.1f} to {d['Lat'].max():.1f} |"
        for c, d in g)
    pairs = F[F['shared_spectrum']].groupby('spectrum_group')['Sample_ID'].apply(
        lambda s: '/'.join(map(str, s)))
    bad_s = ', '.join(map(str, F.loc[F['bad_salinity'], 'Sample_ID']))
    nonpos = ', '.join(map(str, F.loc[F['nonpositive_Rrs'], 'Sample_ID']))
    return f"""# AMT_Rrs_HPLC — AMT above-water Rrs + HPLC matchups

**Source:** compiled and shared by Tom Jordan; received by J. X. Prochaska on
2026-10-24 (`~/Downloads/AMT_Rrs_HPLC.csv`).
**File:** `{NAME}`,
SHA-256 `{sha}`.
**Staged by:** EPFT-UP `scripts/sdp/stage_amt.py` on
{datetime.date.today().isoformat()}. The CSV is unmodified.
**Sharing:** for the AI_for_Ocean group. Do not redistribute beyond it, and
check with J. X. Prochaska / Tom Jordan before publishing results. Citation
and credit terms are pending (asked in EPFT-UP
`requests/AMT_Rrs_HPLC_questions.md`).

## Contents

- **{len(df)} rows**, one per HPLC sample matched to an above-water
  hyperspectral Rrs spectrum.
- **Meta columns:** `Sample_ID`, `Original_ID`, `Cruise`, `Data_link`,
  `HPLC_lab`, `Rad_method`, `DateTime` (MM/DD/YYYY HH:MM, UTC assumed),
  `Lat`, `Lon`, `Temp` [°C], `Sal`.
- **24 HPLC pigments** [mg m⁻³]: `Tchla Tchlb Tchlc ABcaro But Hex Allo
  Diadino Diato Fuco Perid Zea MVchla DVchla MVchlb DVchlb Chlc12 Chlc3 Lut
  Neo Viola Pras Psc Ppc`. `But` = 19′-butanoyloxyfucoxanthin, `Hex` =
  19′-hexanoyloxyfucoxanthin; Psc/Ppc are photosynthetic/photoprotective
  carotenoids. Missing values are the string `nan`.
- **Rrs:** `Rrs_400` … `Rrs_700` [sr⁻¹], 1 nm.

| Cruise | rows | dates | HPLC lab | latitude |
|---|---|---|---|---|
{rows}

## Known issues (see `AMT_Rrs_HPLC_flags.csv`, one row per input row)

1. **Shared spectra:** {int(F['shared_spectrum'].sum())} rows ({len(pairs)} pairs)
   reuse one Rrs spectrum with two different HPLC samples, giving
   {F['spectrum_group'].nunique()} unique spectra. Pairs (Sample_ID):
   {', '.join(pairs)}. Treat each pair as one observation, or keep both rows
   in the same CV fold; otherwise train/test leakage follows.
2. **Non-positive Rrs** (red, 600–700 nm) in {int(F['nonpositive_Rrs'].sum())}
   rows: {nonpos}.
3. **Implausible salinity** ({df['Sal'].min():.2f}, open South Atlantic):
   Sample_ID {bad_s}.
4. **MVchlb and DVchlb missing** for every cruise except AMT29/AMT30 (NASA
   GSFC HPLC).
5. **Below-detection handling differs from Kramer et al. (2022):** almost
   no zeros, values down to ~1e-4 mg m⁻³. Sub-LOD values appear to be
   reported as measured; per-lab LODs are not included.
6. **Radiometry is unsmoothed and unrounded.** Several cruises show a 4 nm
   periodicity in the 2nd derivative, which suggests resampling from a
   coarser native grid. Instruments and processing are not documented in
   the file.
7. **Four HPLC labs** (PML, DHI, NASA GSFC, FURG): lab offsets are
   confounded with cruise.

## EPFT-UP processing defaults (until Tom's answers arrive)

- 5 nm moving mean before the GSM residual, as for Kramer's deposit.
- One row per spectrum, HPLC averaged over shared pairs.
- Below-LOD values zeroed using Kramer's deposit LOD proxy.
- MVchlb: Tchlb as a flagged proxy; MVchlb scored on AMT29/30 only.
- Pair Sample_ID 107/108 dropped (one spectrum, Tchla 0.073 vs 0.165).
- Bad salinity replaced by WOA23 climatology.
- Rows with non-positive Rrs excluded.
- CV groups: each cruise is a group; also a leave-all-AMT-out test.

The rules are implemented in EPFT-UP `epft_up/sdp/amt.py::load_amt` and
specified in `reports/provenance/AMT_benchmark_protocol.md` §3–§4.

## Derived files (written by EPFT-UP `scripts/sdp/amt_harmonise.py`)

- `AMT_Rrs_HPLC_harmonised.csv`: one row per usable spectrum after the rules
  above. It holds meta, pigments under EPFT-UP names (HexFuco, ButFuco),
  flags, the as-delivered `Rrs_400…700`, and Tchla from a GSM-like inversion
  and from OC4v6, for both 5 nm-smoothed (`_smooth5`) and as-delivered
  (`_asis`) spectra. `sample_id` = `AMT:<smallest Sample_ID>`, and
  `sample_ids` lists every collapsed row.
- `AMT_Kramer_folds.csv`: `sample_id, set, group, fold` for the Kramer et
  al. (2022) deposit (`K22:<row>`, 8 campaigns, folds 0–7) and AMT (8
  cruises, folds 8–15).
- `AMT_Rrs_HPLC_harmonised.json`: provenance, rows per cruise after each
  rule, diagnostics, checksums.
- `woa23_salinity_fix.csv` / `.json`: the WOA23 salinity used for the
  implausible values.

The inspection is documented in the EPFT-UP report
(`reports/SDP_Claude_Report.md`, Appendix A) and in
`reports/figures/sdp/amt_inspect_summary.json`.
"""


def stage(dest, src, sha, F, text):
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / NAME
    if target.exists():
        if sdpdata.sha256_of(target) != sha:
            raise SystemExit(f'{target} exists with a different checksum; '
                             'stage the new version under a dated name')
    else:
        shutil.copy2(src, target)
    if sdpdata.sha256_of(target) != sha:
        raise SystemExit(f'checksum mismatch after copy: {target}')
    F.to_csv(dest / 'AMT_Rrs_HPLC_flags.csv', index=False)
    (dest / 'README.md').write_text(text, encoding='utf-8')
    if dest != LOCAL:
        for name in DERIVED:
            if (LOCAL / name).exists():
                shutil.copy2(LOCAL / name, dest / name)
                if sdpdata.sha256_of(dest / name) != sdpdata.sha256_of(LOCAL / name):
                    raise SystemExit(f'checksum mismatch after copy: {dest / name}')
                print(f'  + {name}')
    print(f'staged {target}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--drive', action='store_true', help='also copy to the shared drive')
    a = ap.parse_args()
    src = SRC if SRC.exists() else LOCAL / NAME
    sha = sdpdata.sha256_of(src)
    df = pd.read_csv(src, na_values=['nan'])
    F = flags(df)
    text = readme(df, F, sha)
    stage(LOCAL, src, sha, F, text)
    if a.drive:
        if not DRIVE.parent.exists():
            raise SystemExit(f'shared drive folder not found: {DRIVE.parent}')
        stage(DRIVE, LOCAL / NAME, sha, F, text)
    print(f"sha256 {sha}; {int(F['any_flag'].sum())} flagged rows; "
          f"{F['spectrum_group'].nunique()} unique spectra")


if __name__ == '__main__':
    main()
