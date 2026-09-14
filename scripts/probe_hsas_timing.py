"""Time one AMT24 day through the MOANA Level-2 pipeline (``process_day``).

Written for prompt Move/3 (2026-09-13) to put a measured number on the cost
of the full re-derivation that Q14 asks for as the Move's acceptance bar.
Imports ``epft_up.moana``; before the Move (2026-09-13) it fell back to the
IOPtics copy, which is where the 0.13 s/day figure in Q18 was measured.

Run with::

    conda run -n ocean14 python scripts/probe_hsas_timing.py

Needs the AMT24 HyperSAS Level-2 tree under ``$OS_COLOR``.
"""
import time

import os
import sys

# Run-as-a-script support: make the repo root importable without an install.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from epft_up.moana.io import hsas_day_paths          # noqa: E402
from epft_up.moana.pipeline import process_day       # noqa: E402


def main():
    """Process one mid-cruise day and print its attrition and wall-clock."""
    paths = hsas_day_paths()
    print(f"{len(paths)} delivery days")
    path = paths[len(paths) // 2]        # a mid-cruise day, not the odd DOY 266
    t0 = time.perf_counter()
    out = process_day(path)
    dt = time.perf_counter() - t0
    print(f"{path.name}: {out['attrition']}  ->  {dt:.2f} s")
    print(f"naive whole-cruise estimate: {dt * len(paths):.1f} s")


if __name__ == "__main__":
    main()
