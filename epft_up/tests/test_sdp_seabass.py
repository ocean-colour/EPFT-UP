"""Tests for the SeaBASS reader (Execution #9)."""
import numpy as np

from epft_up.sdp import seabass

SB = """/begin_header
/fields=date,time,lat,lon,depth,Tot_Chl_a,Zea
/units=yyyymmdd,hh:mm:ss,degrees,degrees,m,mg/m^3,mg/m^3
/missing=-9999
/below_detection_limit=-8888
/delimiter=comma
/end_header
20210505,15:24:40,49.1,-14.5,8.375,0.960,-8888
20210505,15:24:40,49.1,-14.5,8.375,1.036,-9999
"""


def test_read_sb(tmp_path):
    f = tmp_path / 'a.sb'
    f.write_text(SB)
    h, df, bdl = seabass.read_sb(f)
    assert h['missing'] == '-9999'
    assert np.isclose(df['Tot_Chl_a'].mean(), 0.998)
    assert bdl['Zea'].tolist() == [True, False]
    assert df['Zea'].isna().all()
    assert str(df['datetime'].iloc[0]) == '2021-05-05 15:24:40+00:00'


def test_index_dir_dedup(tmp_path):
    (tmp_path / 'a.sb').write_text(SB)
    (tmp_path / '0123456789_a.sb').write_text(SB)
    idx = seabass.index_dir(tmp_path)
    assert len(idx) == 1 and idx['n_copies'].iloc[0] == 2
