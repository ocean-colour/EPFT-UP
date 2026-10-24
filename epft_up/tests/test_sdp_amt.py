"""Tests for the AMT loader and harmonisation (Test on new AMT data #1)."""
import numpy as np
import pandas as pd
import pytest

from epft_up.sdp import amt, spectral
from epft_up.sdp import validate as V

WAVE = np.arange(400, 701)


def _row(sid, oid, cruise, rrs, sal=35.5, **pig):
    base = {'Sample_ID': sid, 'Original_ID': oid, 'Cruise': cruise,
            'Data_link': 'x', 'HPLC_lab': 'PML', 'Rad_method': 'Above-water',
            'DateTime': '10/09/2013 12:00', 'Lat': 10.0, 'Lon': -20.0, 'Temp': 20.0,
            'Sal': sal}
    p = {c: 0.01 for c in amt.PIG_MAP}
    p.update({'Tchla': 0.1, 'MVchlb': np.nan, 'DVchlb': np.nan, 'Tchlb': 0.02})
    p.update(pig)
    base.update(p)
    base.update({f'Rrs_{w}': v for w, v in zip(WAVE, rrs)})
    return base


def _spec(scale):
    return scale * (0.004 * np.exp(-(WAVE - 400) / 120.0) + 1e-4)


@pytest.fixture
def raw(tmp_path):
    rows = [
        _row(1, 'a', 'AMT23', _spec(1.0)),
        _row(2, 'b1', 'AMT23', _spec(2.0), Tchla=0.2, Fuco=0.004),
        _row(3, 'b2', 'AMT23', _spec(2.0), Tchla=0.4, Fuco=0.008),   # shares spectrum with 2
        _row(107, 'm1', 'AMT29', _spec(3.0), Tchla=0.073),            # mismatched pair
        _row(108, 'm2', 'AMT29', _spec(3.0), Tchla=0.165),
        _row(9, 'neg', 'AMT29', np.r_[_spec(1.5)[:250], -1e-5 * np.ones(51)]),
        _row(132, '13', 'AMT30', _spec(4.0), sal=22.23, MVchlb=0.005, DVchlb=0.001),
        _row(5, 'lod', 'AMT30', _spec(5.0), Perid=0.003, Neo=0.0, MVchlb=0.004,
             DVchlb=0.0),
        _row(6, 'dhi', 'AMT27', _spec(6.0), DVchla=0.0),
    ]
    p = tmp_path / 'AMT_Rrs_HPLC.csv'
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_checksum_refusal(raw):
    with pytest.raises(ValueError, match='SHA-256'):
        amt.load_amt(path=raw)
    with pytest.raises(ValueError, match='SHA-256'):
        amt.load_amt(path=raw, raw=True)


def test_raw_keeps_everything(raw):
    d = amt.load_amt(path=raw, raw=True, verify=False)
    assert d.n == 9
    assert 'HexFuco' in d.pigments and 'ButFuco' in d.pigments


def test_rules(raw):
    d = amt.load_amt(path=raw, verify=False, sal_fix={132: 36.1})
    ids = set(d.meta['sample_id'])
    # 107/108 dropped, non-positive excluded, 2/3 collapsed to AMT:2
    assert ids == {'AMT:1', 'AMT:2', 'AMT:132', 'AMT:5', 'AMT:6'}
    m = d.meta.set_index('sample_id')
    assert m.loc['AMT:2', 'n_hplc'] == 2 and m.loc['AMT:2', 'sample_ids'] == '2;3'
    pig = d.pigments.set_index(d.meta['sample_id'])
    assert pig.loc['AMT:2', 'Tchla'] == pytest.approx(0.3)
    assert pig.loc['AMT:2', 'Fuco'] == pytest.approx(0.006)
    # salinity replaced
    assert m.loc['AMT:132', 'sal'] == pytest.approx(36.1)
    assert m.loc['AMT:132', 'sal_source'] == 'WOA23'
    # steps record counts per rule
    assert d.steps['total'].tolist() == [9, 7, 6, 5, 5, 5, 5, 5]


def test_missing_salinity_fix_raises(raw):
    with pytest.raises((ValueError, FileNotFoundError)):
        amt.load_amt(path=raw, verify=False)


def test_lod_zeroing_and_mvchlb_proxy(raw):
    d = amt.load_amt(path=raw, verify=False, sal_fix={132: 36.1})
    pig = d.pigments.set_index(d.meta['sample_id'])
    fl = d.flags.set_index(d.meta['sample_id'])
    # Perid 0.003 < LOD 0.004 -> 0; Neo already 0 counts as zeroed
    assert pig.loc['AMT:5', 'Perid'] == 0.0 and fl.loc['AMT:5', 'lod_zeroed:Perid']
    assert fl.loc['AMT:5', 'lod_zeroed:Neo']
    # values at/above LOD survive
    assert pig.loc['AMT:1', 'Fuco'] == pytest.approx(0.01)
    # MVchlb proxy from Tchlb where missing, flagged, not scorable
    assert fl.loc['AMT:1', 'MVchlb_proxy'] and not fl.loc['AMT:1', 'MVchlb_scorable']
    assert pig.loc['AMT:1', 'MVchlb'] == pytest.approx(0.02)
    # measured MVchlb on AMT30 is scorable
    assert fl.loc['AMT:132', 'MVchlb_scorable'] and not fl.loc['AMT:132', 'MVchlb_proxy']


def test_ratio_half_lod_rule(raw):
    d = amt.load_amt(path=raw, verify=False, sal_fix={132: 36.1})
    pig = d.pigments.set_index(d.meta['sample_id'])
    lod = amt.LOD_PROXY_DEPOSIT['Perid']
    r = np.log10(V.replace_zeros(np.array([pig.loc['AMT:5', 'Perid']]), lod)
                 / pig.loc['AMT:5', 'Tchla'])
    assert r[0] == pytest.approx(np.log10(0.5 * lod / 0.1))


def test_lod_proxy_matches_deposit():
    from epft_up.sdp import data as sdpdata
    try:
        k = sdpdata.load_kramer2022()
    except (FileNotFoundError, KeyError, OSError):
        pytest.skip('Kramer deposit not available')
    assert V.lod_proxy(k.pigments13) == pytest.approx(amt.LOD_PROXY_DEPOSIT)


def test_moving_mean_shrink_edges():
    y = np.arange(10.0)
    m = spectral.moving_mean(y, 5, edges='shrink')
    assert m[0] == pytest.approx(1.0) and m[1] == pytest.approx(1.5)
    assert m[-1] == pytest.approx(8.0)
    np.testing.assert_allclose(m[2:-2], spectral.moving_mean(y, 5)[2:-2])
    with pytest.raises(ValueError):
        spectral.moving_mean(y, 5, edges='pad')


def test_dvchla_unseparated_cruise_is_missing(raw):
    d = amt.load_amt(path=raw, verify=False, sal_fix={132: 36.1})
    pig = d.pigments.set_index(d.meta['sample_id'])
    fl = d.flags.set_index(d.meta['sample_id'])
    assert np.isnan(pig.loc['AMT:6', 'DVchla'])
    assert not fl.loc['AMT:6', 'DVchla_scorable']
    assert not fl.loc['AMT:6', 'lod_zeroed:DVchla']
    assert fl.loc['AMT:1', 'DVchla_scorable']
