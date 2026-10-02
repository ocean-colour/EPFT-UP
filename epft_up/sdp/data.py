"""Ingest of the Kramer et al. (2022) HPLC + hyperspectral Rrs matchups.

The only data source is the PANGAEA deposit that accompanies the paper:

    Kramer, S.J., Siegel, D.A., Maritorena, S., Catlett, D. (2021). Global
    surface ocean HPLC phytoplankton pigments and hyperspectral remote sensing
    reflectance. PANGAEA, doi:10.1594/PANGAEA.937536 (CC-BY-4.0).

It holds the **145 quality-controlled** samples of Kramer 2022 (the 33
spectra removed by visual inspection are not deposited), with Rrs on a 1 nm
400-700 nm grid, 27 HPLC pigment columns, and in-situ temperature/salinity.

Two public code repositories are used as references only (never vendored;
fetched by ``scripts/sdp/fetch_kramer2022.py`` and pinned in
:data:`REF_REPOS`). Neither has a LICENSE file, but the README of
``sashajane19/Rrs_pigments`` states that "all code and data in this repository
are freely available for use by anyone for any and all applications".

Everything lives under ``$OS_COLOR/PANGAEA/Kramer2022/`` (Q&A #17, #23)::

    PANGAEA_937536.tab          the deposit (PANGAEA "textfile" export)
    ref/Rrs_pigments/           sashajane19/Rrs_pigments (MATLAB original)
    ref/rrs-SDP-pigments/       max-danenhower/rrs-SDP-pigments (Python port)

No processing happens here beyond renaming and typing: the Rrs are returned
exactly as deposited, and pigment zeros (below detection, set to 0 by the
authors) are left as zeros. Zero handling for log-ratios lives elsewhere
(Q&A #20).
"""

from __future__ import annotations

import datetime
import hashlib
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

#: DOI of the deposit.
PANGAEA_DOI = '10.1594/PANGAEA.937536'
#: Export URL used by the fetch script (PANGAEA tab-delimited "textfile").
PANGAEA_URL = f'https://doi.pangaea.de/{PANGAEA_DOI}?format=textfile'
#: File name under the Kramer2022 directory.
PANGAEA_FILE = 'PANGAEA_937536.tab'
#: SHA-256 of the export as downloaded on 2026-10-02 (every number in
#: reports/SDP_Claude_Report.md derives from this copy).
PANGAEA_SHA256 = '63da2981596ef8c3d856e0ae12eab4cb715a9ec535e92d83458a7a7f5f910f3b'

#: Reference repositories, pinned to the commits checked out on 2026-10-02.
REF_REPOS = {
    'Rrs_pigments': {
        'url': 'https://github.com/sashajane19/Rrs_pigments.git',
        'sha': 'b3e366225d16c4fe85c7406cd3d57e554b8cb57f',
        'note': 'MATLAB original (Kramer); A,B table, aw, GSM, rrsModelTrain.m'},
    'rrs-SDP-pigments': {
        'url': 'https://github.com/max-danenhower/rrs-SDP-pigments.git',
        'sha': 'fb17c2f580dea9488b6ba11c690fe8c959e0a0fe',
        'note': 'Python port (Danenhower); trained A/C coefficient sheets'},
}

#: Native wavelength grid of the deposit [nm].
WAVE = np.arange(400, 701, dtype=float)

#: The 13 pigments modelled in Kramer 2022 (paper abbreviations, paper order)
#: -> PANGAEA column-name prefix.
PIGMENTS_13 = {
    'Tchla': 'TChl a',
    'HexFuco': 'Hex-fuco',
    'ButFuco': 'But-fuco',
    'Allo': 'Allo',
    'Fuco': 'Fuco',
    'Perid': 'Perid',
    'Zea': 'Zea',
    'DVchla': 'DV chl a',
    'MVchlb': 'MV chl b',
    'Chlc12': 'Chl c1+c2',
    'Chlc3': 'Chl c3',
    'Neo': 'Neo',
    'Viola': 'Viola',
}

#: The other deposited pigments (not modelled in Kramer 2022) -> PANGAEA prefix.
#: Abbreviations follow Kramer & Siegel (2019).
PIGMENTS_EXTRA = {
    'Tchlb': 'Chl b + DV Chl b',
    'Tchlc': 'Chl c1+c2+c3',
    'ABcaro': 'a-Car + b-Car',
    'Diadino': 'Diadino',
    'Diato': 'Diato',
    'MVchla': 'MV chl a',
    'Chllide': 'Chlide a',
    'DVchlb': 'DV chl b',
    'Lut': 'Lut',
    'Phytin': 'Phaeophytin',
    'Phide': 'Phaeopho a',
    'Pras': 'Pras',
}

#: PANGAEA "Campaign" value -> the campaign label of Kramer 2022 Table 1.
#: The three Polarstern legs are pooled as "ANT" in the paper.
CAMPAIGN_MAP = {
    'ANT-XXIV/4 (PS71)': 'ANT',
    'ANT-XXV/1 (PS73)': 'ANT',
    'ANT-XXVI/4 (PS75)': 'ANT',
    'NAAMES': 'NAAMES',
    'RemSensPOC': 'RemSensPOC',
    'SABOR': 'SABOR',
    'Tara Oceans': 'Tara Oceans',
    'Tara Mediterranean': 'Tara Med',
    'BIOSOPE': 'BIOSOPE',
    'EXPORTS': 'EXPORTS',
}

#: Kramer 2022 Table 1, transcribed: retained samples and Tchla statistics
#: [mg m^-3] per campaign.
TABLE1 = pd.DataFrame(
    [('ANT', 26, 0.033, 4.15, 0.232, 0.648),
     ('NAAMES', 11, 0.094, 0.987, 0.496, 0.540),
     ('RemSensPOC', 27, 0.049, 1.09, 0.090, 0.173),
     ('SABOR', 9, 0.070, 1.31, 0.252, 0.471),
     ('Tara Oceans', 16, 0.021, 0.950, 0.168, 0.194),
     ('Tara Med', 29, 0.026, 0.170, 0.055, 0.064),
     ('BIOSOPE', 23, 0.019, 1.47, 0.069, 0.326),
     ('EXPORTS', 4, 0.172, 0.292, 0.224, 0.228)],
    columns=['campaign', 'n', 'chl_min', 'chl_max', 'chl_median', 'chl_mean'],
).set_index('campaign')


def kramer2022_dir():
    """Directory holding the Kramer 2022 inputs: ``$OS_COLOR/PANGAEA/Kramer2022``.

    Returns
    -------
    Path
    Raises
    ------
    RuntimeError : if ``$OS_COLOR`` is not set.
    """
    root = os.environ.get('OS_COLOR')
    if root is None:
        raise RuntimeError(
            "The SDP loaders need the $OS_COLOR data tree; the variable is not set.")
    return Path(root) / 'PANGAEA' / 'Kramer2022'


def sha256_of(path, chunk=1 << 22):
    """SHA-256 hex digest of a file, streamed.

    Parameters
    ----------
    path : str or Path
    chunk : int, optional
        Read size in bytes.

    Returns
    -------
    str
    """
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def ref_repo_status(ref_dir=None):
    """Commit SHA of each checked-out reference repo, and whether it is the pin.

    Parameters
    ----------
    ref_dir : str or Path, optional
        Defaults to ``kramer2022_dir() / 'ref'``.

    Returns
    -------
    dict
        ``{name: {'path', 'sha', 'pinned_sha', 'matches_pin'}}``; ``sha`` is
        ``None`` when the clone is missing.
    """
    ref_dir = Path(ref_dir) if ref_dir is not None else kramer2022_dir() / 'ref'
    out = {}
    for name, info in REF_REPOS.items():
        path = ref_dir / name
        sha = None
        if (path / '.git').exists():
            try:
                sha = subprocess.run(
                    ['git', '-C', str(path), 'rev-parse', 'HEAD'],
                    capture_output=True, text=True, check=True).stdout.strip()
            except (OSError, subprocess.CalledProcessError):
                sha = None
        out[name] = {'path': str(path), 'sha': sha, 'pinned_sha': info['sha'],
                     'matches_pin': sha == info['sha']}
    return out


def _header_length(path):
    """Number of lines in the PANGAEA ``/* ... */`` metadata block (inclusive)."""
    with open(path, encoding='utf-8') as fh:
        for i, line in enumerate(fh, start=1):
            if line.startswith('*/'):
                return i
    raise ValueError(f'{path}: no PANGAEA "*/" header terminator found')


def _column(df, prefix):
    """The single column of ``df`` whose name starts with ``prefix + ' ['``."""
    hits = [c for c in df.columns if c.startswith(prefix + ' [')]
    if len(hits) != 1:
        raise KeyError(f'expected one PANGAEA column for {prefix!r}, found {hits}')
    return hits[0]


@dataclass
class KramerData:
    """The Kramer 2022 matchups as deposited.

    Attributes
    ----------
    wave : ndarray, shape (n_wave,)
        Wavelength [nm] (1 nm, 400-700).
    Rrs : ndarray, shape (n_samples, n_wave)
        Above-water remote-sensing reflectance [sr^-1], as deposited.
    meta : DataFrame
        One row per sample: ``campaign`` (Table 1 label), ``sub_campaign``
        (PANGAEA label), ``pi``, ``url_ref``, ``time`` (UTC), ``lat``, ``lon``,
        ``depth_min``, ``depth_max`` [m], ``temp`` [degC], ``sal`` [PSU].
    pigments : DataFrame
        HPLC pigments [mg m^-3 = ug L^-1]; the 13 of :data:`PIGMENTS_13`
        first, then :data:`PIGMENTS_EXTRA`. Zeros mean below detection.
    provenance : dict
        Source, checksum, reference-repo commits, loader version, load time.
    """
    wave: np.ndarray
    Rrs: np.ndarray
    meta: pd.DataFrame
    pigments: pd.DataFrame
    provenance: dict = field(default_factory=dict)

    @property
    def n(self):
        """Number of samples."""
        return self.Rrs.shape[0]

    @property
    def pigments13(self):
        """The 13 pigments modelled in Kramer 2022 (DataFrame view)."""
        return self.pigments[list(PIGMENTS_13)]


def parse_pangaea(path):
    """Parse a PANGAEA 937536 tab export into a :class:`KramerData`.

    No checksum or provenance; see :func:`load_kramer2022`.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    KramerData
    """
    path = Path(path)
    df = pd.read_csv(path, sep='\t', skiprows=_header_length(path))

    rrs_cols = [c for c in df.columns if c.startswith('Rrs_')]
    wave = np.array([float(c.split('_')[1].split()[0]) for c in rrs_cols])
    if not np.array_equal(wave, WAVE):
        raise ValueError(f'{path}: unexpected Rrs wavelength grid '
                         f'{wave[0]}..{wave[-1]} (n={wave.size})')
    Rrs = df[rrs_cols].to_numpy(dtype=float)

    unknown = sorted(set(df['Campaign']) - set(CAMPAIGN_MAP))
    if unknown:
        raise ValueError(f'{path}: unmapped campaigns {unknown}')
    meta = pd.DataFrame({
        'campaign': df['Campaign'].map(CAMPAIGN_MAP),
        'sub_campaign': df['Campaign'],
        'pi': df['PI'],
        'url_ref': df['URL ref'],
        'time': pd.to_datetime(df['Date/Time'], utc=True),
        'lat': df['Latitude'].astype(float),
        'lon': df['Longitude'].astype(float),
        'depth_min': df['Depth water [m] (min)'].astype(float),
        'depth_max': df['Depth water [m] (max)'].astype(float),
        'temp': df[_column(df, 'Temp')].astype(float),
        'sal': df[[c for c in df.columns if c.startswith('Sal (')][0]].astype(float),
    })

    pig_map = {**PIGMENTS_13, **PIGMENTS_EXTRA}
    pigments = pd.DataFrame({short: df[_column(df, prefix)].astype(float)
                             for short, prefix in pig_map.items()})

    return KramerData(wave=wave, Rrs=Rrs, meta=meta, pigments=pigments)


def load_kramer2022(path=None, verify=True):
    """Load the Kramer 2022 PANGAEA deposit with a provenance record.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to ``kramer2022_dir() / PANGAEA_FILE``.
    verify : bool, optional
        If True (default), require the file's SHA-256 to equal
        :data:`PANGAEA_SHA256`.

    Returns
    -------
    KramerData

    Raises
    ------
    FileNotFoundError : the file is missing (run scripts/sdp/fetch_kramer2022.py).
    ValueError : ``verify`` and the checksum differs from the pin.
    """
    path = Path(path) if path is not None else kramer2022_dir() / PANGAEA_FILE
    if not path.is_file():
        raise FileNotFoundError(
            f'{path} not found; run scripts/sdp/fetch_kramer2022.py')
    sha = sha256_of(path)
    if verify and sha != PANGAEA_SHA256:
        raise ValueError(
            f'{path}: SHA-256 {sha} != pinned {PANGAEA_SHA256}. PANGAEA may '
            f'have re-rendered the export; inspect it and re-pin deliberately.')

    data = parse_pangaea(path)

    from epft_up import __version__
    data.provenance = {
        'source': f'PANGAEA doi:{PANGAEA_DOI} (Kramer et al. 2021), CC-BY-4.0',
        'file': str(path),
        'sha256': sha,
        'sha256_pinned': PANGAEA_SHA256,
        'sha256_verified': bool(verify),
        'n_samples': int(data.n),
        'wave_range_nm': [float(data.wave[0]), float(data.wave[-1])],
        'ref_repos': ref_repo_status(path.parent / 'ref'),
        'loader': 'epft_up.sdp.data.load_kramer2022',
        'epft_up_version': __version__,
        'loaded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec='seconds'),
    }
    return data


def table1_comparison(data):
    """Per-campaign counts and Tchla statistics vs Kramer 2022 Table 1.

    Parameters
    ----------
    data : KramerData

    Returns
    -------
    DataFrame
        Indexed by campaign; columns ``n``, ``chl_min``, ``chl_max``,
        ``chl_median``, ``chl_mean`` from the data, the same with suffix
        ``_paper`` from :data:`TABLE1`, and ``n_ok`` (exact count match).
    """
    chl = data.pigments['Tchla']
    g = chl.groupby(data.meta['campaign'])
    ours = pd.DataFrame({'n': g.size(), 'chl_min': g.min(), 'chl_max': g.max(),
                         'chl_median': g.median(), 'chl_mean': g.mean()})
    out = ours.join(TABLE1, rsuffix='_paper').reindex(TABLE1.index)
    out['n_ok'] = out['n'] == out['n_paper']
    return out


def second_derivative_qc(wave, Rrs, band=(610.0, 660.0), threshold=2e-4,
                         step=None):
    """Lange et al. (2020)-style noise screen on the Rrs second derivative.

    Lange et al. (2020, Opt. Express 28, 25682) rejected spectra whose second
    derivative fell outside +/-2e-4 (quoted as sr^-1 nm^-1) between 610 and
    660 nm, after interpolating to a 2 nm grid. This is our objective stand-in
    for Kramer's visual QC of the same band (Q&A #6).

    Parameters
    ----------
    wave : ndarray, shape (n_wave,)
        Uniform grid [nm].
    Rrs : ndarray, shape (n_samples, n_wave)
    band : (float, float), optional
        Inclusive wavelength window tested [nm].
    threshold : float, optional
        Maximum allowed ``|d2 Rrs / d lambda2|``.
    step : float, optional
        Resample to this spacing [nm] (linear interpolation) before
        differencing; ``None`` uses the native grid. Lange used 2 nm.

    Returns
    -------
    reject : ndarray of bool, shape (n_samples,)
    max_abs_d2 : ndarray, shape (n_samples,)
        ``max |d2 Rrs|`` within ``band`` for each spectrum.
    """
    wave = np.asarray(wave, dtype=float)
    Rrs = np.atleast_2d(np.asarray(Rrs, dtype=float))
    if step is not None:
        grid = np.arange(wave[0], wave[-1] + 0.5 * step, step)
        Rrs = np.vstack([np.interp(grid, wave, r) for r in Rrs])
        wave = grid
    dl = np.diff(wave)
    if not np.allclose(dl, dl[0]):
        raise ValueError('second_derivative_qc needs a uniform wavelength grid')
    dl = dl[0]
    # Centred 2nd-order finite difference (Catlett & Siegel 2018, Eq. 2).
    d2 = (Rrs[:, 2:] - 2.0 * Rrs[:, 1:-1] + Rrs[:, :-2]) / dl**2
    w_mid = wave[1:-1]
    sel = (w_mid >= band[0]) & (w_mid <= band[1])
    max_abs = np.max(np.abs(d2[:, sel]), axis=1)
    return max_abs > threshold, max_abs
