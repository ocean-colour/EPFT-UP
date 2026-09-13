"""Validation of MOANA against the three locked targets (design doc §1).

- :func:`validate_amt24` — target (i): reproduce Lange et al. (2020)
  Tables 1–2 on AMT24, in the CTD-only configuration we have (Q&A #35): the
  published NASA model applied to our matchups, our own retrained full-fit
  (Lange "Arrangement 1"), and the 80/20 bootstrap cross-validation
  ("Arrangement 2"). Reference numbers from the paper ship in
  :data:`LANGE_TABLE1` so results diff against the target in one place.
- :func:`validate_heldout_cruises` — target (ii): the published model on the
  Brewin et al. (2023) in-situ hyperspectral Rrs for AMT23/25/28. Retrievals
  run now; *scoring* additionally needs those cruises' flow-cytometry
  deposits, which are not yet on disk (see the function's docstring).
- :func:`bitexact_pace` — target (iii)a: reproduce one PACE L4M MOANA granule
  from its L3M AOP Rrs input under both PC mappings (`nasa_compat`), which
  settles Q&A #9 empirically. Needs `earthaccess` + ``~/.netrc``.
- :func:`match_pace_to_insitu` — target (iii)b machinery: match daily L4M
  MOANA composites to any in-situ cell-count table (the Bailey & Werdell
  departure documented per Q&A #15). The SeaBASS counts themselves are a
  pending data acquisition.

Metrics (:func:`seegers_metrics`): multiplicative log₁₀-space bias and MAE
(Seegers et al. 2018) plus R² on the taxon's model scale, with clipped/
non-positive retrievals treated as **censored** — excluded from the log-space
scores but reported as a headline ``frac_unphysical`` (Q&A #13).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from epft_up.moana.algorithm import run_moana
from epft_up.moana.io import (load_brewin2023, load_fcm, load_luts,
                              load_uway_sst, surface_fcm)
from epft_up.moana.pipeline import build_training_matrix, process_cruise
from epft_up.moana.train import train_moana

#: Lange et al. (2020) Table 1 — hyperspectral, in-situ Rrs, AMT24.
#: bias/MAE multiplicative (log10-space); cv_mae is their bootstrap 80/20.
LANGE_TABLE1 = {
    'pro':   {'n': 73, 'bias': 1.08, 'mae': 1.31, 'r2': 0.82, 'cv_mae': 1.35},
    'syn':   {'n': 73, 'bias': 1.00, 'mae': 1.27, 'r2': 0.92, 'cv_mae': 1.36},
    'peuk':  {'n': 78, 'bias': 1.00, 'mae': 1.21, 'r2': 0.95, 'cv_mae': 1.26},
}

#: Lange Table 2, the CTD-only retraining sensitivity (hyperspectral): the
#: fair comparison line for our CTD-only configuration (design §5.6).
LANGE_CTD_ONLY = {'syn': {'bias': 0.84, 'mae': 1.37}}

_TAXA = ('pro', 'syn', 'peuk')


def seegers_metrics(pred, obs, scale='log10'):
    """Log-space bias/MAE (Seegers et al. 2018) with censored handling.

    Parameters
    ----------
    pred, obs : (n,) arrays — predicted and observed abundances [cells mL⁻¹].
    scale : {'log10', 'linear'} — scale for the R² statistic (the bias/MAE
        are always multiplicative log-space, as in Lange's tables).

    Returns
    -------
    dict — ``n`` usable pairs, ``bias`` (10^mean(Δlog10); 1.08 = +8 %),
    ``mae`` (10^mean|Δlog10|), ``r2`` (on ``scale``),
    ``frac_unphysical`` (fraction of finite predictions ≤ 0 — NASA's clip
    regime; censored out of the log-space scores per Q&A #13),
    ``n_censored``.
    """
    pred = np.asarray(pred, dtype=np.float64)
    obs = np.asarray(obs, dtype=np.float64)
    finite = np.isfinite(pred) & np.isfinite(obs) & (obs > 0)
    censored = finite & (pred <= 0)
    use = finite & (pred > 0)
    out = {'n': int(use.sum()),
           'n_censored': int(censored.sum()),
           'frac_unphysical': (float(censored.sum() / finite.sum())
                               if finite.any() else np.nan)}
    if use.sum() < 3:
        out.update(bias=np.nan, mae=np.nan, r2=np.nan)
        return out
    dlog = np.log10(pred[use]) - np.log10(obs[use])
    out['bias'] = float(10 ** dlog.mean())
    out['mae'] = float(10 ** np.abs(dlog).mean())
    if scale == 'linear':
        p, o = pred[use], obs[use]
    else:
        p, o = np.log10(pred[use]), np.log10(obs[use])
    ss_res = float(((o - p) ** 2).sum())
    ss_tot = float(((o - o.mean()) ** 2).sum())
    out['r2'] = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return out


def _predict_from_trained(model, rrs, sst):
    """Evaluate a :func:`train_moana` model on spectra (in-sample or held out).

    Parameters
    ----------
    model : dict — from :func:`train_moana`.
    rrs : (n, 124) array — raw Rrs on the MOANA grid.
    sst : (n,) array — SST [°C].

    Returns
    -------
    dict — (n,) prediction arrays per taxon (cells mL⁻¹; Pro raw floats).
    """
    from epft_up.moana.algorithm import standardize
    rrs_std = standardize(rrs)
    scores = (rrs_std - model['pca']['mean']) @ model['pca']['V']
    preds = {}
    for taxon in _TAXA:
        m = model['models'][taxon]
        y = np.full(rrs_std.shape[0], m['coef'][0])
        for name, b in zip(m['names'], m['coef'][1:]):
            if name == 'logSST':
                y = y + b * np.log10(np.asarray(sst, dtype=np.float64))
            else:
                y = y + b * scores[:, int(name[1:]) - 1]
        preds[taxon] = y if m['response'] == 'linear' else 10.0 ** y
    return preds


def validate_amt24(config=None, n_boot=100, seed=0, verbose=True):
    """Target (i): reproduce Lange Tables 1–2 on AMT24 (CTD-only config).

    Runs the full Level-2 chain, builds the matchup table, then scores three
    things per taxon against the flow-cytometry truth:

    - ``published`` — the operational NASA model (LUT coefficients,
      operational PC mapping) applied to our matchup Rrs + `uway_sst`;
    - ``retrained`` — our own PCA + stepwise regressions, in-sample
      (Lange Arrangement 1, "full-fit");
    - ``cv`` — bootstrap 80/20 cross-validation of the retraining
      (Arrangement 2): ``n_boot`` refits on random 80 % subsets, scored on
      the held-out 20 %.

    Parameters
    ----------
    config : dict, optional — pipeline overrides (e.g. Lange-strict
        ``{'match_window_min': 0}``).
    n_boot : int, optional — bootstrap refits (Lange's count is unstated).
    seed : int, optional — RNG seed for the bootstrap splits.
    verbose : bool, optional — print a comparison table against
        :data:`LANGE_TABLE1`.

    Returns
    -------
    dict — ``matchups`` (DataFrame), ``published`` / ``retrained`` / ``cv``
    (per-taxon metric dicts), ``model`` (the retrained model),
    ``reference`` (:data:`LANGE_TABLE1`).
    """
    stream = process_cruise(verbose=False)
    tm = build_training_matrix(stream, surface_fcm(load_fcm()),
                               load_uway_sst(), config)
    mu, rrs = tm['matchups'], tm['rrs']
    lut = load_luts()
    scales = {'pro': 'linear', 'syn': 'log10', 'peuk': 'log10'}

    # -- published operational model ---------------------------------------
    ret = run_moana(lut['wave'], rrs, sst=mu['sst'].to_numpy(), lut=lut)
    published = {t: seegers_metrics(ret[{'peuk': 'apeuk'}.get(t, t)],
                                    mu[t].to_numpy(), scales[t])
                 for t in _TAXA}

    # -- our retraining, in-sample (Arrangement 1) --------------------------
    counts = {t: mu[t].to_numpy() for t in _TAXA}
    sst = mu['sst'].to_numpy()
    model = train_moana(rrs, counts, sst=sst)
    preds = _predict_from_trained(model, rrs, sst)
    retrained = {t: seegers_metrics(preds[t], counts[t], scales[t])
                 for t in _TAXA}

    # -- bootstrap 80/20 cross-validation (Arrangement 2) -------------------
    rng = np.random.default_rng(seed)
    n = len(mu)
    cv_pred = {t: [] for t in _TAXA}
    cv_obs = {t: [] for t in _TAXA}
    for _ in range(n_boot):
        train_idx = rng.choice(n, size=int(round(0.8 * n)), replace=False)
        test_idx = np.setdiff1d(np.arange(n), train_idx)
        try:
            m_b = train_moana(rrs[train_idx],
                              {t: counts[t][train_idx] for t in _TAXA},
                              sst=sst[train_idx])
        except ValueError:
            continue        # a starved resample; skip it
        p_b = _predict_from_trained(m_b, rrs[test_idx], sst[test_idx])
        for t in _TAXA:
            cv_pred[t].append(p_b[t])
            cv_obs[t].append(counts[t][test_idx])
    cv = {t: seegers_metrics(np.concatenate(cv_pred[t]),
                             np.concatenate(cv_obs[t]), scales[t])
          for t in _TAXA}

    out = {'matchups': mu, 'published': published, 'retrained': retrained,
           'cv': cv, 'model': model, 'reference': LANGE_TABLE1}
    if verbose:
        print(f"target (i) — AMT24, CTD-only, n = {len(mu)} matchups "
              f"(Lange n = 73–78 incl. underway FCM)")
        hdr = f"{'taxon':6s} {'mode':10s} {'n':>4s} {'bias':>6s} {'MAE':>6s} " \
              f"{'R2':>6s}   Lange T1: bias/MAE/R2"
        print(hdr)
        for t in _TAXA:
            ref = LANGE_TABLE1[t]
            for mode, res in (('published', published[t]),
                              ('retrained', retrained[t]), ('cv', cv[t])):
                print(f"{t:6s} {mode:10s} {res['n']:4d} {res['bias']:6.2f} "
                      f"{res['mae']:6.2f} {res['r2']:6.2f}   "
                      f"{ref['bias']:.2f}/{ref['mae']:.2f}/{ref['r2']:.2f}")
    return out


def validate_heldout_cruises(fcm_tables=None, verbose=True):
    """Target (ii): the published model on AMT23/25/28 in-situ Rrs.

    Retrievals need only the Brewin et al. (2023) deposit (on disk).
    **Scoring needs those cruises' flow-cytometry counts, which are not yet
    downloaded** — BODC deposits 10.5285/a2104adc-e990- (AMT23), -e98e-
    (AMT25), a147c314-688b- (AMT28); browser-only, like the AMT24 ones
    (Q&A #24). Pass them via ``fcm_tables`` when available.

    Parameters
    ----------
    fcm_tables : dict, optional
        ``{cruise_number: DataFrame}`` with columns datetime/lat/lon/depth +
        pro/syn/peuk [cells mL⁻¹] (the :func:`epft_up.moana.io.load_fcm`
        layout). Stations are matched to the nearest count sample within
        ±3 h; scoring is skipped for cruises without a table.

    Returns
    -------
    dict — ``retrievals`` : DataFrame (cruise, datetime, lat, lon, sst,
    pro/syn/apeuk, flags, recon_residual); ``metrics`` : per-taxon metric
    dicts per cruise plus a ``'pooled'`` entry across all scored cruises;
    ``missing_counts`` : cruise numbers lacking a truth table;
    ``pairs`` : DataFrame (cruise, taxon, pred, obs) — one row per
    matchup, i.e. the points behind the metrics (used by the report
    figures).
    """
    brewin = load_brewin2023()
    st = brewin['stations']
    ret = run_moana(brewin['wave'], brewin['rrs'], sst=st['sst'].to_numpy())
    retrievals = st.copy()
    for k in ('pro', 'syn', 'apeuk'):
        retrievals[k] = ret[k]
    retrievals['flags'] = ret['flags']
    retrievals['recon_residual'] = ret['recon_residual']

    metrics, missing = {}, []
    pooled = {t: ([], []) for t in _TAXA}      # (pred, obs) across cruises
    pair_rows = []                             # every matchup, for plotting
    for cruise in sorted(st['cruise'].unique()):
        table = (fcm_tables or {}).get(cruise)
        if table is None:
            missing.append(int(cruise))
            continue
        sub = retrievals[retrievals['cruise'] == cruise]
        surf = table[table['depth'] <= 10.0]
        t_obs = pd.to_datetime(surf['datetime']).astype('datetime64[ns]')
        rows = {}
        for t in _TAXA:
            pred, obs = [], []
            for _, r in sub.iterrows():
                dt = (t_obs - pd.Timestamp(r['datetime'])).abs()
                j = dt.idxmin()
                if dt[j] <= pd.Timedelta(hours=3):
                    pred.append(r[{'peuk': 'apeuk'}.get(t, t)])
                    obs.append(surf.loc[j, t])
            rows[t] = seegers_metrics(
                np.array(pred), np.array(obs),
                'linear' if t == 'pro' else 'log10')
            pooled[t][0].extend(pred)
            pooled[t][1].extend(obs)
            pair_rows.extend(
                {'cruise': int(cruise), 'taxon': t,
                 'pred': float(pp), 'obs': float(oo)}
                for pp, oo in zip(pred, obs))
        metrics[int(cruise)] = rows
    if metrics:
        metrics['pooled'] = {
            t: seegers_metrics(np.array(pooled[t][0]), np.array(pooled[t][1]),
                               'linear' if t == 'pro' else 'log10')
            for t in _TAXA}
    if verbose:
        print(f"target (ii) — {len(retrievals)} stations retrieved on "
              f"cruises {sorted(st['cruise'].unique())}; "
              f"counts missing for {missing}")
        for c, rows in metrics.items():
            for t, m in rows.items():
                print(f"  {str(c):6s} {t:5s} n={m['n']:3d} "
                      f"bias {m['bias']:5.2f}  MAE {m['mae']:5.2f}  "
                      f"R2 {m['r2']:6.2f}  unphys {m['frac_unphysical']:.2f}")
    return {'retrievals': retrievals, 'metrics': metrics,
            'missing_counts': missing,
            'pairs': pd.DataFrame(pair_rows)}


# ---------------------------------------------------------------------------
# Target (iii): the operational PACE product
# ---------------------------------------------------------------------------

#: In-file variable names in PACE_OCI_L4M_MOANA (misspelled; report §7.4),
#: and the product's non-value sentinels (report Appendix A).
PACE_VARS = {'pro': 'prococcus_moana', 'syn': 'syncoccus_moana',
             'peuk': 'picoeuk_moana'}
PACE_LAND = 254
PACE_FILL = -32767


def _pace_dir():
    """Local cache directory for PACE granules (under the data tree)."""
    root = os.environ.get('OS_COLOR')
    if root is None:
        raise RuntimeError('$OS_COLOR is not set')
    d = Path(root) / 'PACE' / 'moana_validation'
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch_pace_pair(date='2025-07-01'):
    """Download one day's L3M AOP Rrs + L4M MOANA 0.1° granule pair.

    Parameters
    ----------
    date : str, optional — the day (both products daily, 0.1°).

    Returns
    -------
    (aop_path, moana_path) : Paths to the local netCDF files (cached —
    already-downloaded files are not re-fetched). Needs ``~/.netrc``.
    """
    import earthaccess
    earthaccess.login(strategy='netrc')
    out = {}
    for short_name, tag in (('PACE_OCI_L3M_AOP', 'aop'),
                            ('PACE_OCI_L4M_MOANA', 'moana')):
        found = earthaccess.search_data(short_name=short_name,
                                        temporal=(date, date), count=40)
        picks = [g for g in found
                 if '.DAY.' in g['umm']['GranuleUR']
                 and '0p1deg' in g['umm']['GranuleUR']]
        if not picks:
            raise RuntimeError(f'no daily 0.1° granule for {short_name} {date}')
        paths = earthaccess.download(picks[:1], str(_pace_dir()))
        out[tag] = Path(paths[0])
    return out['aop'], out['moana']


def bitexact_pace(date='2025-07-01', max_pixels=200_000, seed=0, verbose=True):
    """Target (iii)a: reproduce NASA's MOANA granule from its Rrs input.

    Runs our retrieval in ``nasa_compat`` mode under **both** PC mappings on
    the L3M AOP Rrs, on the MOANA regional grid, and compares pixel-for-pixel
    with the shipped L4M product. *Synechococcus* and picoeukaryotes need no
    SST, so their comparison is ancillary-free and is the discriminating test
    for the Q&A #9 mapping question; Prochlorococcus is skipped (its SST
    field, GHRSST CMC, is not reproduced here).

    **Exact equality is structurally unreachable** (found empirically,
    2026-08-16): NASA retrieves MOANA at L2 and *then* composites to L4M,
    while this test retrieves from the already-composited L3M Rrs — a
    nonlinear-operation ordering difference (and the likely reason no L3M
    MOANA product exists). The mapping question is still cleanly decidable in
    log space: the correct mapping tracks NASA to within the compositing
    scatter (median |Δlog₁₀| ≲ 0.01) while the wrong one is offset by the
    relocated coefficient (~+0.08 median on Synechococcus). The verdict
    below therefore uses the median |Δlog₁₀| of positive pixels.

    Parameters
    ----------
    date : str, optional — granule day.
    max_pixels : int, optional — cap on ocean pixels retrieved (random
        subsample for speed; the verdict needs far fewer).
    seed : int, optional — subsample seed.
    verbose : bool, optional.

    Returns
    -------
    dict — per mapping, per taxon (syn/peuk): ``n`` compared,
    ``frac_exact`` (identical int32 after our truncation), ``frac_within_1``
    (|Δ| ≤ 1 cell mL⁻¹, the truncation-boundary tolerance); plus
    ``verdict`` — which mapping reproduces NASA (or 'inconclusive').
    """
    import xarray as xr
    aop_path, moana_path = fetch_pace_pair(date)
    aop = xr.open_dataset(aop_path)
    moana = xr.open_dataset(moana_path)

    # The MOANA regional grid is an exact subset of the AOP global grid
    # (report Appendix B): align by coordinate value. 'nearest' with a
    # hundredth-of-a-cell tolerance absorbs float32 storage rounding while
    # still failing loudly on any genuine grid mismatch.
    rrs = aop['Rrs'].sel(lat=moana['lat'].values, lon=moana['lon'].values,
                         method='nearest', tolerance=1e-3)
    wave = aop['wavelength'].values.astype(np.float64)

    nasa = {t: moana[PACE_VARS[t]].values for t in ('syn', 'peuk')}
    # Real-retrieval mask: not land-254, not fill, finite, in NASA's product.
    valid = np.isfinite(nasa['syn']) & np.isfinite(nasa['peuk'])
    for t in ('syn', 'peuk'):
        valid &= (nasa[t] != PACE_LAND) & (nasa[t] != PACE_FILL) & (nasa[t] >= 0)
    iy, ix = np.nonzero(valid)
    if iy.size > max_pixels:
        pick = np.random.default_rng(seed).choice(iy.size, max_pixels,
                                                  replace=False)
        iy, ix = iy[pick], ix[pick]
    spectra = rrs.values[iy, ix, :]   # (n, n_wave)

    results = {}
    for mapping in ('operational', 'atbd'):
        ours = run_moana(wave, spectra, sst=None, pc_mapping=mapping,
                         nasa_compat=True)
        res = {}
        for t in ('syn', 'peuk'):
            mine = ours[{'peuk': 'apeuk'}.get(t, t)]
            theirs = nasa[t][iy, ix].astype(np.float64)
            ok = np.isfinite(mine)
            d = mine[ok] - theirs[ok]
            pos = ok & (mine > 0)
            dlog = (np.log10(mine[pos])
                    - np.log10(nasa[t][iy, ix][pos].astype(np.float64)))
            res[t] = {'n': int(ok.sum()),
                      'frac_exact': float((d == 0).mean()),
                      'frac_within_1': float((np.abs(d) <= 1).mean()),
                      'median_dlog': float(np.median(dlog)),
                      'mad_dlog': float(np.median(np.abs(dlog)))}
        results[mapping] = res
    aop.close(), moana.close()

    # Verdict on Synechococcus (the taxon whose disputed coefficient moves):
    # the correct mapping sits inside the L2-vs-L3M compositing scatter, the
    # wrong one carries the relocated coefficient as a systematic offset.
    op, at = (results[m]['syn']['mad_dlog'] for m in ('operational', 'atbd'))
    if min(op, at) > 0.02 or max(op, at) < 3 * min(op, at):
        verdict = 'inconclusive'
    else:
        verdict = 'operational' if op < at else 'atbd'
    results['verdict'] = verdict
    if verbose:
        for m in ('operational', 'atbd'):
            r = results[m]
            print(f"{m:12s}: syn med dlog {r['syn']['median_dlog']:+.4f} "
                  f"(MAD {r['syn']['mad_dlog']:.4f}, exact "
                  f"{r['syn']['frac_exact']:.3f})  "
                  f"peuk med dlog {r['peuk']['median_dlog']:+.4f}  "
                  f"n={r['syn']['n']}")
        print("verdict:", verdict)
    return results


def match_pace_to_insitu(insitu, dates=None, window_days=0, verbose=True):
    """Target (iii)b machinery: daily L4M MOANA vs an in-situ count table.

    Uses daily 0.1° composites rather than Bailey & Werdell's L2 5×5 boxes —
    the documented protocol departure (Q&A #15). The in-situ table is caller-
    supplied because the SeaBASS picophytoplankton counts are a pending data
    acquisition (SeaBASS account needed; Q&A #24).

    Parameters
    ----------
    insitu : pandas.DataFrame — columns datetime/lat/lon + any of
        pro/syn/peuk [cells mL⁻¹], surface samples only.
    dates : iterable of str, optional — restrict to these days (default: the
        unique days in ``insitu``).
    window_days : int, optional — also accept composites within ± this many
        days of the sample (0 = same day only).
    verbose : bool, optional.

    Returns
    -------
    dict — ``matchups`` : DataFrame (one row per sample with a valid
    same-cell product value: observed + product values per taxon, distance
    to cell centre); ``metrics`` : per-taxon :func:`seegers_metrics` with
    NASA's clipped zeros censored (Q&A #13).
    """
    import xarray as xr
    insitu = insitu.copy()
    insitu['datetime'] = pd.to_datetime(insitu['datetime'])
    if dates is None:
        dates = sorted({d.strftime('%Y-%m-%d') for d in insitu['datetime']})
    rows = []
    for date in dates:
        day = pd.Timestamp(date)
        sel = insitu[(insitu['datetime'] - day).abs()
                     <= pd.Timedelta(days=window_days, hours=23, minutes=59)]
        if not len(sel):
            continue
        try:
            _, moana_path = fetch_pace_pair(date)
        except RuntimeError:
            continue
        ds = xr.open_dataset(moana_path)
        for _, s in sel.iterrows():
            cell = ds.sel(lat=s['lat'], lon=s['lon'], method='nearest')
            row = {'datetime': s['datetime'], 'lat': s['lat'],
                   'lon': s['lon'], 'granule_date': date}
            good = False
            for t in _TAXA:
                v = float(cell[PACE_VARS[t]].values)
                bad = (v == PACE_LAND) or (v == PACE_FILL) or not np.isfinite(v)
                row[f'{t}_pace'] = np.nan if bad else v
                row[f'{t}_obs'] = s.get(t, np.nan)
                good |= not bad
            if good:
                rows.append(row)
        ds.close()
    matchups = pd.DataFrame(rows)
    metrics = {}
    for t in _TAXA:
        if len(matchups) and matchups[f'{t}_obs'].notna().any():
            metrics[t] = seegers_metrics(
                matchups[f'{t}_pace'].to_numpy(),
                matchups[f'{t}_obs'].to_numpy(),
                'linear' if t == 'pro' else 'log10')
    if verbose:
        print(f"target (iii)b — {len(matchups)} matchups over "
              f"{len(dates)} days")
    return {'matchups': matchups, 'metrics': metrics}
