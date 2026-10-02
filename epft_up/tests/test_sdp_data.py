"""Tests for :mod:`epft_up.sdp.data` (Kramer 2022 PANGAEA ingest).

Tier 1 runs on a tiny synthetic PANGAEA-format file; Tier 2 (``needs_kramer2022``)
checks the real deposit against Kramer et al. (2022) Table 1.
"""

import numpy as np
import pytest

from epft_up.sdp import data as sdpdata
from epft_up.tests.conftest import needs_kramer2022


def _write_synthetic_pangaea(path, n=3):
    """Write a minimal PANGAEA-style tab export with ``n`` rows."""
    pig_cols = [f'{p} [µg/l] (High Performance Liquid Chrom...)'
                for p in {**sdpdata.PIGMENTS_13, **sdpdata.PIGMENTS_EXTRA}.values()]
    rrs_cols = [f'Rrs_{int(w)} [1/sr] (Hyperspectral radiometer)' for w in sdpdata.WAVE]
    cols = (['Campaign', 'URL ref', 'PI', 'Date/Time', 'Latitude', 'Longitude',
             'Depth water [m] (min)', 'Depth water [m] (max)'] + pig_cols
            + ['Temp [°C] (CTD)', 'Sal (PSU, CTD)'] + rrs_cols)
    camps = ['ANT-XXIV/4 (PS71)', 'Tara Mediterranean', 'EXPORTS']
    lines = ['/* DATA DESCRIPTION:', 'Citation:\tsynthetic', '*/', '\t'.join(cols)]
    for i in range(n):
        row = [camps[i % 3], 'https://example.org', 'PI', '2010-01-01T00:00:00',
               str(10.0 + i), str(-20.0 - i), '0', '7']
        row += [f'{0.1 * (i + 1):.3f}'] * len(pig_cols)
        row += ['15.0', '35.0']
        row += [f'{0.005 * np.exp(-(w - 400) / 150):.6f}' for w in sdpdata.WAVE]
        lines.append('\t'.join(row))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return path


def test_parse_synthetic(tmp_path):
    f = _write_synthetic_pangaea(tmp_path / 'synth.tab')
    d = sdpdata.parse_pangaea(f)
    assert d.n == 3
    assert d.Rrs.shape == (3, sdpdata.WAVE.size)
    np.testing.assert_array_equal(d.wave, sdpdata.WAVE)
    assert list(d.pigments.columns[:13]) == list(sdpdata.PIGMENTS_13)
    assert d.pigments13.shape == (3, 13)
    assert list(d.meta['campaign']) == ['ANT', 'Tara Med', 'EXPORTS']
    assert d.meta['time'].dt.tz is not None
    np.testing.assert_allclose(d.meta['temp'], 15.0)


def test_checksum_mismatch_raises(tmp_path):
    f = _write_synthetic_pangaea(tmp_path / 'synth.tab')
    with pytest.raises(ValueError, match='SHA-256'):
        sdpdata.load_kramer2022(f, verify=True)
    d = sdpdata.load_kramer2022(f, verify=False)
    assert d.provenance['sha256'] == sdpdata.sha256_of(f)
    assert d.provenance['sha256_verified'] is False


def test_unmapped_campaign_raises(tmp_path):
    f = _write_synthetic_pangaea(tmp_path / 'synth.tab')
    f.write_text(f.read_text(encoding='utf-8').replace('EXPORTS', 'NEWCRUISE'),
                 encoding='utf-8')
    with pytest.raises(ValueError, match='unmapped campaigns'):
        sdpdata.parse_pangaea(f)


def test_campaign_map_covers_table1():
    assert set(sdpdata.CAMPAIGN_MAP.values()) == set(sdpdata.TABLE1.index)
    assert sdpdata.TABLE1['n'].sum() == 145


def test_second_derivative_qc_synthetic():
    wave = sdpdata.WAVE
    smooth = 0.002 * np.exp(-(wave - 400) / 100)
    spiky = smooth.copy()
    spiky[np.searchsorted(wave, 640)] += 5e-4   # d2 ~ 1e-3 >> 2e-4
    rej, mx = sdpdata.second_derivative_qc(wave, np.vstack([smooth, spiky]))
    assert rej.tolist() == [False, True]
    assert mx[1] > 2e-4 > mx[0]
    # resampling path (Lange's 2 nm grid) still flags the spike
    rej2, _ = sdpdata.second_derivative_qc(wave, np.vstack([smooth, spiky]), step=2.0)
    assert rej2.tolist() == [False, True]


def test_second_derivative_qc_quadratic_exact():
    wave = sdpdata.WAVE
    a = 3e-6
    _, mx = sdpdata.second_derivative_qc(wave, a * (wave - 500.0)**2)
    np.testing.assert_allclose(mx, 2 * a, rtol=1e-6)


@needs_kramer2022
def test_real_deposit_matches_table1():
    d = sdpdata.load_kramer2022()
    assert d.n == 145
    assert d.provenance['sha256'] == sdpdata.PANGAEA_SHA256
    assert np.isfinite(d.Rrs).all() and (d.Rrs > 0).all()
    assert not d.pigments.isna().any().any()
    t1 = sdpdata.table1_comparison(d)
    assert t1['n_ok'].all()
    # Table 1 Tchla statistics, to the paper's quoted precision
    for col in ('chl_min', 'chl_max', 'chl_median', 'chl_mean'):
        np.testing.assert_allclose(t1[col], t1[col + '_paper'], atol=0.0051)
    assert d.pigments['Tchla'].min() == pytest.approx(0.019)
    assert d.pigments['Tchla'].median() == pytest.approx(0.110)
