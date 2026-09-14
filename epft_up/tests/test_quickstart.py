"""The README quickstart (``scripts/quickstart_moana.py``) keeps working.

Tier 2: needs the 2025-07-01 PACE granule pair cached under
``$OS_COLOR/PACE/moana_validation/``. The network step (``fetch_pace_pair``)
is replaced by the cached paths so the test never touches Earthdata.
"""
import importlib.util
import os
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT = _ROOT / 'scripts' / 'quickstart_moana.py'


def _cached_pair():
    root = os.getenv('OS_COLOR')
    if root is None:
        return None
    d = Path(root) / 'PACE' / 'moana_validation'
    aop = d / 'PACE_OCI.20250701.L3m.DAY.AOP.V3_2.0p1deg.nc'
    moana = d / 'PACE_OCI.20250701.L4m.DAY.MOANA.V3_2.0p1deg.nc'
    return (aop, moana) if aop.is_file() and moana.is_file() else None


needs_pace_cache = pytest.mark.skipif(
    _cached_pair() is None,
    reason='requires the cached 2025-07-01 PACE granule pair under $OS_COLOR/PACE')


@needs_pace_cache
def test_quickstart_retrieves_a_clean_pixel(capsys):
    spec = importlib.util.spec_from_file_location('quickstart_moana', _SCRIPT)
    qs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(qs)
    qs.fetch_pace_pair = lambda date='2025-07-01': _cached_pair()   # no network

    qs.main()
    out = capsys.readouterr().out
    lines = {l.split()[0]: l for l in out.strip().splitlines()}
    assert lines['flags'].split()[1] == '0', out
    for taxon in ('pro', 'syn', 'apeuk'):
        value = float(lines[taxon].split()[1])
        assert value > 0, out
    # The README quotes this pixel; keep the two in step.
    assert lines['pixel'].split()[1:3] == ['31.55°N', '60.05°W'], out
