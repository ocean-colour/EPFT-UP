"""Fetch and verify the Kramer et al. (2022) inputs under ``$OS_COLOR``.

* PANGAEA doi:10.1594/PANGAEA.937536 (CC-BY-4.0) -> downloaded as the PANGAEA
  tab "textfile" export and checked against the pinned SHA-256.
* The two reference repositories (``sashajane19/Rrs_pigments``,
  ``max-danenhower/rrs-SDP-pigments``) -> cloned into ``ref/`` and checked out
  at the pinned commits. Nothing from them is copied into EPFT-UP (Q&A #23).

Pins live in :mod:`epft_up.sdp.data` (``PANGAEA_SHA256``, ``REF_REPOS``).

Run with::

    conda run -n ocean14 python scripts/sdp/fetch_kramer2022.py            # fetch + verify
    conda run -n ocean14 python scripts/sdp/fetch_kramer2022.py --verify   # verify only
"""
import argparse
import os
import subprocess
import sys
import urllib.request

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                os.pardir, os.pardir)))

from epft_up.sdp import data as sdpdata  # noqa: E402


def fetch_pangaea(dest, verify_only=False):
    """Download the deposit if missing; return True if it matches the pin."""
    if not dest.is_file():
        if verify_only:
            print(f'MISSING  {dest}')
            return False
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f'GET      {sdpdata.PANGAEA_URL}')
        tmp = dest.with_suffix('.part')
        urllib.request.urlretrieve(sdpdata.PANGAEA_URL, tmp)
        tmp.rename(dest)
    sha = sdpdata.sha256_of(dest)
    ok = sha == sdpdata.PANGAEA_SHA256
    print(f'{"OK" if ok else "MISMATCH":8s} {dest.name}  sha256={sha}')
    return ok


def fetch_repo(name, info, ref_dir, verify_only=False):
    """Clone ``name`` if missing and check out its pinned commit."""
    path = ref_dir / name
    if not (path / '.git').exists():
        if verify_only:
            print(f'MISSING  {path}')
            return False
        ref_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(['git', 'clone', '-q', info['url'], str(path)], check=True)
    if not verify_only:
        subprocess.run(['git', '-C', str(path), 'fetch', '-q', 'origin'], check=False)
        subprocess.run(['git', '-C', str(path), 'checkout', '-q', info['sha']],
                       check=True)
    status = sdpdata.ref_repo_status(ref_dir)[name]
    print(f'{"OK" if status["matches_pin"] else "MISMATCH":8s} {name}  '
          f'HEAD={status["sha"]}  pin={info["sha"]}')
    return status['matches_pin']


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--verify', action='store_true',
                   help='only verify what is present; download nothing')
    args = p.parse_args(argv)

    root = sdpdata.kramer2022_dir()
    ok = fetch_pangaea(root / sdpdata.PANGAEA_FILE, verify_only=args.verify)
    for name, info in sdpdata.REF_REPOS.items():
        ok &= fetch_repo(name, info, root / 'ref', verify_only=args.verify)
    print('all inputs verified' if ok else 'SOME INPUTS MISSING OR MISMATCHED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
