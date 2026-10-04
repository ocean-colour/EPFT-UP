"""Minimal reader for SeaBASS ``.sb`` files (Execution #9).

SeaBASS files are plain text: a ``/begin_header … /end_header`` block of
``/key=value`` lines (``/fields=``, ``/units=``, ``/missing=``,
``/below_detection_limit=``, ``/delimiter=``, …), then delimited rows.

:func:`read_sb` returns the header dict and a DataFrame with ``missing`` →
NaN, ``below_detection_limit`` → a separate boolean mask (the caller decides
whether that is 0, as Kramer does, or a censoring flag), and a ``datetime``
column built from ``date`` + ``time`` when present.

:func:`index_dir` lists the ``.sb`` files of a download directory with their
SHA-256, strips the browser-bundle prefix (``[0-9a-f]{10}_``) that SeaBASS
adds to bundled files, and keeps one copy of byte-identical duplicates.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd

BUNDLE_PREFIX = re.compile(r'^[0-9a-f]{10}_')


def read_sb(path):
    """Read one SeaBASS file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    header : dict
        ``/key=value`` pairs (keys without the slash, lower case).
    df : DataFrame
        Data with ``missing`` values as NaN; column ``datetime`` (UTC) when
        ``date`` and ``time`` fields exist.
    bdl : DataFrame of bool
        True where a value was flagged ``below_detection_limit`` (set to NaN
        in ``df``).
    """
    path = Path(path)
    header, rows, in_header = {}, [], True
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n').rstrip('\r')
            if in_header:
                if line.lower().startswith('/end_header'):
                    in_header = False
                elif line.startswith('/') and '=' in line:
                    k, v = line[1:].split('=', 1)
                    header[k.strip().lower()] = v.strip()
                continue
            if line.strip() and not line.startswith('!'):
                rows.append(line)
    if 'fields' not in header:
        raise ValueError(f'{path}: no /fields= in header')
    fields = [f.strip() for f in header['fields'].split(',')]
    delim = {'comma': ',', 'space': None, 'tab': '\t'}.get(
        header.get('delimiter', 'comma').lower(), ',')
    parsed = []
    for r in rows:
        v = r.split(delim) if delim else r.split()
        v = [x.strip() for x in v][:len(fields)]
        v += [''] * (len(fields) - len(v))
        parsed.append(v)
    df = pd.DataFrame(parsed, columns=fields)
    for c in df.columns:
        num = pd.to_numeric(df[c], errors='coerce')
        if num.notna().sum() >= max(1, 0.5 * df[c].astype(str).str.len().gt(0).sum()):
            df[c] = num
    miss = float(header.get('missing', '-9999').split('[')[0])
    bdl_val = header.get('below_detection_limit')
    bdl = pd.DataFrame(False, index=df.index, columns=df.columns)
    for c in df.columns:
        if pd.api.types.is_numeric_dtype(df[c]):
            if bdl_val is not None:
                b = np.isclose(df[c], float(bdl_val.split('[')[0]))
                bdl[c] = b
                df.loc[b, c] = np.nan
            df.loc[np.isclose(df[c], miss), c] = np.nan
    df = df.copy()
    if 'date' in df.columns and 'time' in df.columns:
        df['datetime'] = pd.to_datetime(df['date'].astype('Int64').astype(str) + ' '
                                        + df['time'].astype(str), format='%Y%m%d %H:%M:%S',
                                        utc=True, errors='coerce')
    return header, df, bdl


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def index_dir(directory, pattern='*.sb'):
    """Unique ``.sb`` files in ``directory`` (recursive), de-duplicated by content.

    Returns
    -------
    DataFrame
        Columns ``name`` (prefix stripped), ``path``, ``sha256``,
        ``n_copies``; one row per distinct name.
    """
    recs = {}
    for p in sorted(Path(directory).rglob(pattern)):
        name = BUNDLE_PREFIX.sub('', p.name)
        sha = sha256_of(p)
        if name in recs:
            if recs[name]['sha256'] != sha:
                raise ValueError(f'two different files named {name}')
            recs[name]['n_copies'] += 1
            continue
        recs[name] = {'name': name, 'path': str(p), 'sha256': sha, 'n_copies': 1}
    return pd.DataFrame(recs.values()).sort_values('name').reset_index(drop=True)
