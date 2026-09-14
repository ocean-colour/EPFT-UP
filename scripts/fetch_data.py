"""Fetch and verify the public inputs of the MOANA baseline under ``$OS_COLOR``.

Rather than re-depositing datasets that already have DOIs, EPFT-UP pulls them
from their archives and pins their checksums here (Q&A/RTD S2, 2026-09-14):

* **Jordan et al. (2025)** AMT24 netCDF (underway SST) — Zenodo, CC-BY-4.0;
  downloaded directly and verified.
* **PACE OCI** daily 0.1° L3M AOP + L4M MOANA granules — NASA OB.DAAC via
  ``earthaccess`` (needs Earthdata credentials in ``~/.netrc``).
* **BODC** deposits (AMT23/24/25/28 flow cytometry; Brewin et al. 2023 Rrs) —
  NERC Open Data Licence; BODC downloads are browser-only, so this script
  prints the DOI and the expected filename for any that are missing, and
  verifies the ones that are present.

The one input this script cannot fetch is the PML AMT24 HyperSAS Level-2
radiometry (private communication, unpublished); see the Data page.

Run with::

    export OS_COLOR=/path/to/data/root
    conda run -n ocean14 python scripts/fetch_data.py            # fetch + verify
    conda run -n ocean14 python scripts/fetch_data.py --verify   # checksums only
    conda run -n ocean14 python scripts/fetch_data.py --no-pace  # skip Earthdata
"""
import argparse
import hashlib
import os
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

#: Public inputs: path under $OS_COLOR, SHA-256 as measured on 2026-09-14 from
#: the copies used for every number in reports/MOANA_Claude_Report.md, and
#: where they come from.
PUBLIC_FILES = {
    'AMT/AMT23/AMT23_JR20131005_AFC_Dataset.csv': {
        'sha256': '06928683259b74903ba2a34405de6c76788a52a0aa643c9509aa9f2a69db3aab',
        'source': 'BODC doi:10.5285/a2104adc-e990-6789-e053-6c86abc0d557 (Tarran & Zubkov 2020, AMT23)',
        'download': None},
    'AMT/AMT24/AMT24_JR20140922_AFC_Dataset.csv': {
        'sha256': 'f955813d2cd5c58c9583634033da226bd53bdad04da4ec2da8eb3bd6590521cd',
        'source': 'BODC doi:10.5285/a2104adc-e98f-6789-e053-6c86abc0d557 (Tarran & Zubkov 2020, AMT24)',
        'download': None},
    'AMT/AMT25/AMT25_JR15001_AFC_Dataset.csv': {
        'sha256': '62ef207f7274b93b41b6de8aa3414ca0300a89d05dc8d14f2dee1c29aa9f5c39',
        'source': 'BODC doi:10.5285/a2104adc-e98e-6789-e053-6c86abc0d557 (Tarran, Lange & Zubkov 2020, AMT25)',
        'download': None},
    'AMT/AMT28/AMT28_JR18001_AFC_Dataset.csv': {
        'sha256': '53df2349f26cb735f708556900f21322b3e16e242d741d9412ff074aba422f85',
        'source': 'BODC doi:10.5285/a147c314-688b-55e9-e053-6c86abc0dc81 (Tarran & May 2020, AMT28)',
        'download': None},
    'AMT/AMT24/BODC_data_AMT_modern_and_historical_optical_observations.csv': {
        'sha256': '448139d806d396aab016e2911b2b2786e2a14d4e0d10c3defef47d942a11db2a',
        'source': 'BODC doi:10.5285/f3198e10-faf3-1525-e053-6c86abc0d2f6 (Brewin et al. 2023)',
        'download': None},
    'AMT/AMT24/amt24_final_with_debiased_chl.nc': {
        'sha256': '1538c327eefdd1ca1bf3487830043c91badfa30a6f9c9479d71405408666d758',
        'source': 'Zenodo doi:10.5281/zenodo.12527954 (Jordan et al. 2025), CC-BY-4.0, 571 MB',
        'download': 'https://zenodo.org/api/records/12527954/files/amt24_final_with_debiased_chl.nc/content'},
}

#: The granule day every report figure and the target (iii-a) experiment use.
PACE_DATE = '2025-07-01'


def sha256_of(path, chunk=1 << 22):
    """SHA-256 hex digest of a file, streamed."""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    """Stream ``url`` to ``dest`` with a coarse progress line."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + '.part')
    with urllib.request.urlopen(url) as resp, open(tmp, 'wb') as out:
        total = int(resp.headers.get('Content-Length') or 0)
        done = 0
        for block in iter(lambda: resp.read(1 << 22), b''):
            out.write(block)
            done += len(block)
            if total:
                print(f'\r  {dest.name}: {100 * done / total:5.1f} %', end='', flush=True)
    print()
    tmp.replace(dest)


def fetch_public_files(root, verify_only=False):
    """Download (where possible) and verify every entry of :data:`PUBLIC_FILES`.

    Returns
    -------
    (n_ok, n_bad, missing) — verified count, checksum failures, and the
    relative paths that are absent and cannot be downloaded here.
    """
    n_ok, n_bad, missing = 0, 0, []
    for rel, meta in PUBLIC_FILES.items():
        path = root / rel
        if not path.is_file():
            if meta['download'] and not verify_only:
                print(f'fetching {rel}\n  from {meta["download"]}')
                download(meta['download'], path)
            else:
                missing.append(rel)
                print(f'MISSING  {rel}\n  -> {meta["source"]}'
                      + ('' if meta['download'] else '  [browser download from BODC]'))
                continue
        got = sha256_of(path)
        if got == meta['sha256']:
            n_ok += 1
            print(f'ok       {rel}')
        else:
            n_bad += 1
            print(f'CHECKSUM MISMATCH {rel}\n  expected {meta["sha256"]}\n  got      {got}')
    return n_ok, n_bad, missing


def fetch_pace(date=PACE_DATE):
    """Fetch (or confirm cached) the PACE granule pair via the package's own helper."""
    from epft_up.moana.validation import fetch_pace_pair
    aop, moana = fetch_pace_pair(date)
    for p in (aop, moana):
        print(f'ok       PACE/moana_validation/{p.name}  ({p.stat().st_size / 1e6:.0f} MB)')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--root', default=os.environ.get('OS_COLOR'),
                        help='data root (default: $OS_COLOR)')
    parser.add_argument('--verify', action='store_true',
                        help='verify checksums only; download nothing')
    parser.add_argument('--no-pace', action='store_true',
                        help='skip the PACE granules (no Earthdata access)')
    args = parser.parse_args()
    if not args.root:
        sys.exit('set $OS_COLOR or pass --root')
    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    print(f'data root: {root}')

    n_ok, n_bad, missing = fetch_public_files(root, verify_only=args.verify)
    if not args.no_pace and not args.verify:
        try:
            fetch_pace()
        except Exception as exc:                      # no ~/.netrc, offline, ...
            print(f'PACE granules not fetched: {exc}')

    print(f'\n{n_ok} verified, {n_bad} checksum failures, {len(missing)} missing')
    if missing:
        print('BODC files must be downloaded in a browser from the DOIs above and '
              'placed at the paths shown; then re-run with --verify.')
    print('Not fetchable here: the PML AMT24 HyperSAS Level-2 radiometry '
          '(AMT/AMT24/Radiometry/level2/) — available on request from PML.')
    return 1 if n_bad else 0


if __name__ == '__main__':
    sys.exit(main())
