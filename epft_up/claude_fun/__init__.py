"""Claude Fun — reproduce and extend MOANA with a Bayesian set-encoder network
(cf_design.md).

Public surface
--------------
run_cf       checkpoint 1: stages 0–5, training matrix, transfer set, diagnostics
train_cf     train one model on a training matrix
predict_cf   predictive distribution for any spectra (any grid)
validate_cf  checkpoints 2–4: CV, baselines, ablations, robustness, attribution, transfer
"""

import copy
import glob
import json
import os
import time

import numpy as np
import pandas as pd
import xarray as xr

from epft_up.claude_fun import algorithm, attribution, io, pipeline, train, validation
from epft_up.claude_fun.algorithm import DEFAULT_MODEL
from epft_up.claude_fun.attribution import DEFAULT_ATTR
from epft_up.claude_fun.pipeline import DEFAULT_PIPELINE
from epft_up.claude_fun.train import DEFAULT_TRAIN

#: §8.3 ablations: name → (model overrides, train overrides).
ABLATIONS = {
    'corr_diag': ({'corr_mode': 'diag'}, {}),
    'corr_input': ({'corr_mode': 'input'}, {}),
    'deterministic_head': ({'bayesian_head': False}, {}),
    'no_degradation_aug': ({}, {'augment': {'degrade': False}}),
    'sst_off': ({'use_sst': False}, {}),
    'range_414_660': ({'wave_range_nm': (414.0, 660.0)}, {}),
    'counterbin_mlp': ({'encoder': 'counterbin'}, {}),
}


def _json(obj, file):
    with open(file, 'w') as fh:
        json.dump(obj, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o))


def run_cf(pipeline_cfg=None, out_dir=None, transfer=True, verbose=True, hash_files=True):
    """Checkpoint 1: stages 0–5 on AMT24, the Lange-strict variant, the
    AMT23/25/28 transfer set and every checkpoint-1 diagnostic.

    Inputs: pipeline_cfg (default DEFAULT_PIPELINE); out_dir (default
    <AMT_ROOT>/AMT24/derived/claude_fun); transfer (bool); verbose; hash_files.
    Outputs: dict with ``matrix``, ``matrix_strict``, ``transfer``, ``stream``
    and the tables; everything is also written to out_dir (stream.nc,
    matrix.nc, matrix_strict.nc, transfer.nc, attrition.csv,
    transfer_attrition.csv, count_steps.csv, provider_crosscheck.csv,
    dphi_check.csv, diagnostics.json, provenance.json).
    """
    pc = pipeline_cfg or DEFAULT_PIPELINE
    out_dir = out_dir or pipeline.path(pc, 'derived_dir')
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    if verbose:
        print('Stages 0-3 ...', flush=True)
    stream = pipeline.process_cruise(pc, verbose=verbose)
    mat, stream, extras = pipeline.build_training_matrix(pc, stream)
    strict, _, _ = pipeline.build_training_matrix(pc, stream, window_min=0)
    files = json.loads(stream.attrs['files'])
    per_day = pd.DataFrame(json.loads(stream.attrs['per_day'])).sort_values('qc_pass')
    clear = [f for f in files if pipeline.cfio.parse_hsas_name(f)[1] in per_day['doy'].values[-3:]]
    dphi = pipeline.dphi_convention_check([io.read_hsas_day(f) for f in clear])
    steps = pipeline.count_step_table(pc)
    xcheck = pipeline.provider_crosscheck(stream, pc)
    B = io.read_brewin(pipeline.path(pc, 'brewin_file'))
    floor_fit = io.fit_sigma_floor(B['rrs'], B['unc_pct'])
    out = {'matrix': mat, 'matrix_strict': strict, 'stream': stream, 'attrition': extras['attrition'],
           'count_steps': steps, 'provider_crosscheck': xcheck, 'dphi_check': dphi,
           'floor_fit': floor_fit, 'nir_sigma_abs': extras['nir_sigma_abs'],
           'rho_bounds': extras['rho_bounds'], 'surface_counts': extras['surface_counts']}
    if transfer:
        tset, tatt = pipeline.build_transfer_set(pc)
        out.update(transfer=tset, transfer_attrition=tatt)
        pipeline.save_dataset(tset, os.path.join(out_dir, 'transfer.nc'))
        tatt.to_csv(os.path.join(out_dir, 'transfer_attrition.csv'))
    pipeline.save_dataset(stream, os.path.join(out_dir, 'stream.nc'))
    pipeline.save_dataset(mat, os.path.join(out_dir, 'matrix.nc'))
    pipeline.save_dataset(strict, os.path.join(out_dir, 'matrix_strict.nc'))
    extras['attrition'].to_csv(os.path.join(out_dir, 'attrition.csv'), index=False)
    steps.to_csv(os.path.join(out_dir, 'count_steps.csv'))
    xcheck.to_csv(os.path.join(out_dir, 'provider_crosscheck.csv'))
    dphi.to_csv(os.path.join(out_dir, 'dphi_check.csv'), index=False)
    diag = {'nir_sigma_abs': extras['nir_sigma_abs'], 'rho_bounds': extras['rho_bounds'],
            'brewin_floor_fit': floor_fit, 'surface_counts': extras['surface_counts'],
            'dphi_best': dphi.iloc[0].to_dict(), 'seconds': time.time() - t0}
    _json(diag, os.path.join(out_dir, 'diagnostics.json'))
    inputs = files + [pipeline.path(pc, pc['afc_files'][c]) for c in pc['afc_files']] + \
        [pipeline.path(pc, 'brewin_file'), pipeline.path(pc, 'jordan_file')]
    _json(pipeline.provenance(pc, files=inputs, hash_files=hash_files),
          os.path.join(out_dir, 'provenance.json'))
    if verbose:
        print_checkpoint1(out)
    return out


def print_checkpoint1(out):
    """Print the checkpoint-1 summary. Inputs: run_cf output."""
    pd.set_option('display.width', 160)
    print('\n=== Attrition (AMT24) ===')
    print(out['attrition'].to_string(index=False))
    print('\n=== Detection-limit fit vs Stage 5 table ===')
    print(out['count_steps'][['step', 'volume_ml', 'method', 'delta_loglik', 'second_step', 'power',
                              'configured_step']].round(3).to_string())
    if 'hist' in out['count_steps']:
        print('counts of the values 1..19 (rounded files):')
        for cr, h in out['count_steps']['hist'].dropna().items():
            print(f'  {cr}: {h}')
    n = out['nir_sigma_abs']
    print(f"\nNIR-residual σ_abs = {n['sigma_abs_nir']:.2e} sr^-1 (ratio to 8.2e-5: {n['ratio']:.2f}; "
          f"{'OK' if n['within_tolerance'] else 'OUTSIDE 2x: stop and report'})")
    f = out['floor_fit']
    print(f"Brewin σ-floor refit: f_rel = {f['f_rel']:.4f}, σ_abs = {f['sigma_abs']:.2e}")
    print(f"ρ_sky at bounds: {out['rho_bounds']}")
    print('\n=== Δϕ convention check (best first) ===')
    print(out['dphi_check'].round(2).to_string(index=False))
    xc = out['provider_crosscheck']
    sel = xc.loc[(xc.index >= 400) & (xc.index <= 700)].iloc[::10]
    print('\n=== Provider fixed-ρ minus fitted Rrs (400–700 nm, every 10th band) ===')
    print(sel.to_string(float_format=lambda v: f'{v:.2e}'))
    if 'transfer_attrition' in out:
        print('\n=== Transfer set attrition ===')
        print(out['transfer_attrition'].to_string())


def train_cf(matrix, rows=None, model_cfg=None, train_cfg=None, seed=None, verbose=False,
             pipeline_cfg=None):
    """Train one model. Inputs: matrix (Dataset or train.prepare dict), rows,
    configs, seed, verbose. Outputs: model dict (algorithm.save_model to keep)."""
    data = matrix if isinstance(matrix, dict) else train.prepare(matrix, pipeline_cfg)
    return train.train_cf(data, rows, model_cfg, train_cfg, seed, verbose)


def predict_cf(model, inputs, n_samples=None, seed=12345, return_draws=False):
    """Predictive distribution (§6.4). Inputs: model (or list), inputs
    (pipeline Dataset or spectra set), n_samples, seed, return_draws.
    Outputs: xarray.Dataset (algorithm.predict_cf)."""
    spec = inputs if isinstance(inputs, dict) else algorithm.spectra_from_dataset(inputs)
    return algorithm.predict_cf(model, spec, n_samples, seed, return_draws)


def _write_cv(directory, res, data):
    os.makedirs(directory, exist_ok=True)
    rep = validation.cv_report(res)
    rep['metrics'].to_csv(os.path.join(directory, 'metrics.csv'))
    rep['reliability'].to_csv(os.path.join(directory, 'reliability.csv'))
    for k in ('robustness', 'sst_reliance'):
        if k in rep:
            rep[k].to_csv(os.path.join(directory, f'{k}.csv'))
    if rep['pit']:
        np.savez(os.path.join(directory, 'pit.npz'), **rep['pit'])
    pipeline.save_dataset(validation.predictions_to_dataset(res, data),
                          os.path.join(directory, 'predictions.nc'))
    _json(res['models_meta'], os.path.join(directory, 'models_meta.json'))
    _json({'moana_pro_clipped': res.get('moana_clipped', 0)}, os.path.join(directory, 'notes.json'))
    return rep


def validate_cf(matrix, transfer=None, out_dir='cf_results',
                steps=('loso', 'latblock', 'ablations', 'known_answer', 'transfer'),
                model_cfg=None, train_cfg=None, pipeline_cfg=None, attr_cfg=None,
                ablations=None, ablation_cv='latblock', seaflow=None, n_samples=None,
                verbose=True):
    """Checkpoints 2–4.

    Steps
    -----
    loso          leave-one-station-out CV: Claude Fun, climatology, retrained
                  MOANA; SST permutation; degraded grids (checkpoint 2)
    latblock      latitude-blocked 5-fold CV, same, plus attribution maps with
                  the retrained-MOANA Jacobian and the randomisation check
    ablations     §8.3 ablations on ``ablation_cv`` folds
    known_answer  IG and occlusion on retrained MOANA vs its analytic Jacobian
    transfer      train on all AMT24, evaluate on AMT23/25/28 (+ maps); with
                  ``seaflow`` (path to a SeaFlow table) also SeaFlow+PACE

    Inputs: matrix, transfer (Datasets from run_cf), out_dir, steps, configs,
    ablations (dict like ABLATIONS), ablation_cv ('latblock' or 'loso'),
    seaflow, n_samples, verbose.
    Outputs: dict of reports; files under out_dir/<step>/.
    """
    mc = model_cfg or DEFAULT_MODEL
    tc = train_cfg or DEFAULT_TRAIN
    pc = pipeline_cfg or DEFAULT_PIPELINE
    ac = {**DEFAULT_ATTR, **(attr_cfg or {})}
    os.makedirs(out_dir, exist_ok=True)
    data = train.prepare(matrix, pc)
    out = {}
    _json(pipeline.provenance(pc, mc, tc, seeds=tc['ensemble_seeds'] or (tc['seed'],),
                              hash_files=False, attribution=ac, steps=list(steps)),
          os.path.join(out_dir, 'provenance.json'))
    if 'loso' in steps:
        r = validation.run_cv(data, train.folds_loso(data), mc, tc, pc, n_samples=n_samples,
                              verbose=verbose, label='LOSO')
        out['loso'] = _write_cv(os.path.join(out_dir, 'loso'), r, data)
        if verbose:
            print(out['loso']['metrics'].round(3).to_string())
    if 'latblock' in steps:
        calls = []

        def fn(models, d, trn, tst):
            calls.append(1)
            return attribution.fold_maps(models, d, trn, tst, ac, tc, pc, randomize=len(calls) == 1)
        r = validation.run_cv(data, train.folds_latblock(data, tc['cv_lat_blocks']), mc, tc, pc,
                              attribution_fn=fn, n_samples=n_samples, verbose=verbose, label='lat-block')
        d = os.path.join(out_dir, 'latblock')
        out['latblock'] = _write_cv(d, r, data)
        maps = attribution.combine_maps(r['attribution'])
        pipeline.save_dataset(maps, os.path.join(d, 'attribution.nc'))
        strat = attribution.stratify_by_latitude(maps, ac)
        if strat is not None:
            pipeline.save_dataset(strat, os.path.join(d, 'attribution_by_latitude.nc'))
        if verbose:
            print(out['latblock']['metrics'].round(3).to_string())
    if 'ablations' in steps:
        rows = []
        for name, (mo, to) in (ablations or ABLATIONS).items():
            m2 = algorithm.config(**{**mc, **mo})
            t2 = copy.deepcopy(tc)
            for k, v in to.items():
                if isinstance(t2[k], dict):
                    t2[k].update(v)
                else:
                    t2[k] = v
            folds = train.folds_latblock(data, tc['cv_lat_blocks']) if ablation_cv == 'latblock' \
                else train.folds_loso(data)
            r = validation.run_cv(data, folds, m2, t2, pc, robustness=False, sst_perm=False,
                                  baselines=False, n_samples=n_samples, verbose=verbose, label=name)
            s, _ = validation.score(r['cf'], data['labels'])
            rows.append(pd.concat({name: s}, names=['ablation']))
        out['ablations'] = pd.concat(rows)
        os.makedirs(os.path.join(out_dir, 'ablations'), exist_ok=True)
        out['ablations'].to_csv(os.path.join(out_dir, 'ablations', f'metrics_{ablation_cv}.csv'))
        if verbose:
            print(out['ablations'].round(3).to_string())
    if 'known_answer' in steps:
        R124 = validation.to_moana_grid(data, cfg=pc)
        sst = np.asarray(data['spec']['sst'])
        mm = validation.moana_fit(R124, data['labels'], data['ds']['abundance'].values, sst, tc['moana'])
        ka = attribution.moana_known_answer(mm, data['wave'], data['rrs'], sst, np.nanmean(R124, 0), ac, pc)
        out['known_answer'] = ka
        os.makedirs(os.path.join(out_dir, 'known_answer'), exist_ok=True)
        _json(ka, os.path.join(out_dir, 'known_answer', 'moana_known_answer.json'))
        if verbose:
            print(f"MOANA known answer: IG corr median per taxon "
                  f"{np.round(np.nanmedian(ka['ig_corr'], 0), 3)} (pass {ka['pass_ig']}); "
                  f"occlusion vs J·ΔR corr {np.round(np.nanmedian(ka['occ_corr'], 0), 3)}; "
                  f"top-3 overlap {np.round(np.nanmean(ka['occ_top3_overlap'], 0), 2)}; "
                  f"sum rules {ka['sum_JR']:.1e}, {ka['sum_J']:.1e}")
    if 'transfer' in steps and transfer is not None:
        tdata = train.prepare(transfer, pc)
        r = validation.transfer_test(data, tdata, mc, tc, pc, n_samples=n_samples, verbose=verbose)
        d = os.path.join(out_dir, 'transfer')
        os.makedirs(d, exist_ok=True)
        r['per_cruise'].to_csv(os.path.join(d, 'metrics_per_cruise.csv'))
        rep = validation.cv_report({**r, 'fold': np.zeros(tdata['n'], int)})
        for k in ('robustness', 'sst_reliance', 'reliability'):
            if k in rep:
                rep[k].to_csv(os.path.join(d, f'{k}.csv'))
        pipeline.save_dataset(validation.predictions_to_dataset(r, tdata), os.path.join(d, 'predictions.nc'))
        for i, m in enumerate(r['models']):
            algorithm.save_model(m, os.path.join(d, f'model_seed{m["meta"]["seed"]}'))
        R_tr = validation.to_moana_grid(data, cfg=pc)
        sst_tr = np.asarray(data['spec']['sst'])
        mm = validation.moana_fit(R_tr, data['labels'], data['ds']['abundance'].values, sst_tr, tc['moana'])
        baseline = {'wave': data['wave'], 'rrs': np.nanmean(data['rrs'], axis=0)}
        tm = attribution.maps_for(r['models'], tdata['spec'], tdata['labels'], baseline, mm,
                                  validation.to_moana_grid(tdata, cfg=pc), np.asarray(tdata['spec']['sst']),
                                  tdata['lat'], ac, randomize=True)
        pipeline.save_dataset(attribution.combine_maps([tm]), os.path.join(d, 'attribution.nc'))
        out['transfer'] = {'per_cruise': r['per_cruise'], **rep}
        if 'published_moana_error' in r:
            _json({'published_moana_error': r['published_moana_error']}, os.path.join(d, 'notes.json'))
        if verbose:
            print(r['per_cruise'].xs('all', level='cruise').round(3).to_string())
        if seaflow is not None:
            df = io.read_seaflow(seaflow, abundance_factor=pc['seaflow_abundance_factor'])
            sp = validation.seaflow_pace_matchups(df, pc, download_dir=os.path.join(out_dir, 'pace_l2'))
            if sp is not None:
                pipeline.save_dataset(sp, os.path.join(d, 'seaflow_pace_matchups.nc'))
                sdata = train.prepare(sp, pc)
                rs = validation.transfer_test(data, sdata, mc, tc, pc, models=r['models'],
                                              n_samples=n_samples, verbose=verbose)
                rs['per_cruise'].to_csv(os.path.join(d, 'seaflow_pace_metrics.csv'))
                out['seaflow_pace'] = rs['per_cruise']
    return out
