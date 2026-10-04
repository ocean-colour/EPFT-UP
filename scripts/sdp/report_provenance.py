"""Execution #10: provenance manifest and headline-number check for the SDP report.

``reports/SDP_Claude_Report.md`` must trace every number to a script. This
script makes that checkable:

1. **Section map.** For each report section, the scripts that produced it
   and their outputs in ``reports/figures/sdp/``. Each output gets its
   existence, SHA-256 and, for JSON outputs, the ``script`` field it records
   (which must match).
2. **Coverage.** Every ``figures/sdp/*`` file the report references exists,
   every script in ``scripts/sdp/`` is cited in the report, and every output
   file is claimed by a section.
3. **Headline claims.** The numbers quoted in §1 (Introduction) and §10
   (Conclusions) are re-read from the outputs and compared with the value
   printed in the report, to the printed precision.

Output: ``reports/figures/sdp/report_manifest.json``. The script exits
non-zero if any check fails.

Run with::

    conda run -n ocean14 python scripts/sdp/report_provenance.py
"""
import datetime
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
FIG = REPO / 'reports' / 'figures' / 'sdp'
REPORT = REPO / 'reports' / 'SDP_Claude_Report.md'

#: report section -> (scripts, outputs in reports/figures/sdp)
SECTIONS = {
    '3 Data': (['fetch_kramer2022.py', 'ingest_kramer2022.py'],
               ['ingest_summary.json', 'kramer_fig1_map.png', 'kramer_fig2a_rrs.png']),
    '4.1 GSM residual': (['reproduce_gsm.py', 'fetch_woa_ts.py'],
                         ['gsm_summary.json', 'kramer_fig4_chl.png',
                          'kramer_fig2bc_model_residual.png', 'gsm_sensitivity.png']),
    '4.2 PCR': (['reproduce_pcr.py'],
                ['pcr_summary.json', 'kramer_table2.csv', 'pcr_coefficients.png',
                 'kramer_fig6_pigments.png', 'kramer_fig3_dendrograms.png']),
    '5 Maths': (['maths_section.py'],
                ['maths_summary.json', 'maths_effective_weights.png', 'maths_prior.png',
                 'maths_snr.png', 'maths_dof.png', 'maths_filter_factors.png']),
    '6.1-6.2 Diagnostics': (['diagnostics_null_loco.py'],
                            ['diag_summary.json', 'diag_skill_table.csv',
                             'diag_skill_absolute.png', 'diag_skill_ratio.png',
                             'diag_loco_tchla.png']),
    '6.3 Source spaces': (['source_spaces.py'],
                          ['source_spaces_summary.json', 'source_spaces_summary.csv',
                           'source_spaces_table.csv', 'source_spaces_heatmap.png']),
    '7 Learned weighting': (['learned_weighting.py'],
                            ['weighting_summary.json', 'weighting_summary.csv',
                             'weighting_table.csv', 'weighting_ranks.csv',
                             'weighting_latent_spectra.png', 'weighting_target_spectra.png']),
    '8 Uncertainty': (['uncertainty.py'],
                      ['uncertainty_summary.json', 'uncertainty_summary.csv',
                       'uncertainty_scores.csv', 'uncertainty_coverage.png',
                       'uncertainty_pit.png', 'uncertainty_skill.png',
                       'uncertainty_examples.png']),
    '9 EXPORTS-NA hold-out': (['exports_na_matchups.py', 'exports_na_holdout.py'],
                              ['exports_na_summary.json', 'exports_na_scores.csv',
                               'exports_na_vs_null.csv', 'exports_na_tchla.png',
                               'exports_na_ratios.png']),
    '10 Provenance (this script)': (['report_provenance.py'], []),
}


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def jload(name):
    return json.loads((FIG / name).read_text())


def csv(name):
    return pd.read_csv(FIG / name)


def row(df, **kw):
    m = pd.Series(True, index=df.index)
    for k, v in kw.items():
        m &= df[k] == v
    if m.sum() != 1:
        raise KeyError(f'{kw}: {int(m.sum())} rows')
    return df[m].iloc[0]


def claims():
    """(id, report section, printed value, value recomputed from outputs)."""
    pcr = jload('pcr_summary.json')['table2']['n145']
    gsm = jload('gsm_summary.json')['fig4']
    mth = jload('maths_summary.json')
    diag = csv('diag_skill_table.csv').set_index('pigment')
    ss = csv('source_spaces_summary.csv').set_index('model')
    ws = csv('weighting_summary.csv').set_index('model')
    us = csv('uncertainty_summary.csv')
    usc = csv('uncertainty_scores.csv')
    ex = jload('exports_na_summary.json')
    exs = csv('exports_na_scores.csv')
    exv = csv('exports_na_vs_null.csv')
    dof = mth['part_ii']['dof']
    pn = mth['part_ii']['prediction_noise']['Tchla']

    out = [
        ('table2_within_1sd', '4.2', 13,
         sum(abs(v['dR2_in_sd']) <= 1 for v in pcr.values())),
        ('table2_Tchla_R2', '4.2', 0.73, pcr['Tchla']['R2']),
        ('fig4_gsm_slope_n144', '4.1', 0.961, gsm['GSM_without_worst']['slope']),
        ('fig4_gsm_R2_n144', '4.1', 0.864, gsm['GSM_without_worst']['R2']),
        ('weff_power_ge_0.1_Tchla', '5.1', 0.997,
         mth['part_i']['Tchla']['frac_power_f_ge_0.1_eff']),
        ('prior_hi_lo_power_ratio', '5.1', 5300,
         round(mth['implied_prior']['prior_power_ratio_hi_lo_d2'], -2)),
        ('dof_insitu_1nm', '5.2', 10.4,
         dof['pct:0.02 + floor (in-situ-like)']['1.0']['d_s_dRrs']),
        ('dof_pace_white_1nm', '5.2', 6.0, dof['PACE white']['1.0']['d_s_dRrs']),
        ('noise_Tchla_pace_white_d2', '5.2', 11.6,
         pn['PACE white']['ratio_to_pigment_sd_d2']),
        ('noise_Tchla_pace_white_direct', '5.2', 0.83,
         pn['PACE white']['ratio_to_pigment_sd_direct']),
        # derived: pixels N to average so that the δRrs''-PCR Tchla noise SD (∝ 1/√N)
        # falls below the natural SD of Tchla under PACE white noise
        ('pace_white_pixels_for_unit_noise', '10.3', 135,
         math.ceil(pn['PACE white']['ratio_to_pigment_sd_d2']**2)),
        ('null_gsm_beats_pcr_random_abs', '6.1', 11,
         int((diag['abs_linear_random_dR2'] < 0).sum())),
        ('loco_pcr_Tchla_R2', '6.2', 0.44, diag.loc['Tchla', 'abs_linear_loco_R2_PCR']),
        ('loco_nullgsm_Tchla_R2', '6.2', 0.66, diag.loc['Tchla', 'abs_linear_loco_R2_null']),
        ('loco_abs_dRrs_d2', '6.3', 0.19, ss.loc['dRrs_d2', 'abs_R2_loco']),
        ('loco_abs_dRrs', '6.3', 0.38, ss.loc['dRrs', 'abs_R2_loco']),
        ('loco_abs_M1', '6.3', 0.40, ss.loc['M1', 'abs_R2_loco']),
        ('loco_ratio_dRrs_d2', '6.3', 0.17, ss.loc['dRrs_d2', 'ratio_logR2_loco']),
        ('loco_ratio_M3_21', '6.3', 0.27, ss.loc['M3_w21', 'ratio_logR2_loco']),
        ('W_group_loco_abs', '7.2', 0.39, ws.loc['RRR_dRrs_aux_group', 'abs_R2_loco']),
        ('W_group_loco_ratio', '7.2', 0.24, ws.loc['RRR_dRrs_aux_group', 'ratio_logR2_loco']),
        ('W_group_ratio_beats_null_loco', '7.2', 6,
         ws.loc['RRR_dRrs_aux_group', 'ratio_beats_null_loco']),
        ('PCR_ratio_beats_null_loco', '7.2', 1, ws.loc['PCR_dRrs_d2', 'ratio_beats_null_loco']),
        ('W_cov68_clean', '8.2', 0.71,
         row(us, scenario='clean_1nm', model='sharedW_unaware', kind='ratio')['cov68']),
        ('W_cov95_clean', '8.2', 0.94,
         row(us, scenario='clean_1nm', model='sharedW_unaware', kind='ratio')['cov95']),
        ('PCR_ratio_R2_insitu', '8.3', 0.00,
         row(us, scenario='insitu_1nm', model='PCR_dRrs_d2', kind='ratio')['mean_R2']),
        ('Waware_Fuco_pace_white', '8.3', 0.36,
         row(usc, scenario='pace_white_1nm', model='sharedW_aware', target='ratio:Fuco')['R2']),
        ('nullOC4_Fuco_pace_white', '8.3', 0.27,
         row(usc, scenario='pace_white_1nm', model='null_OC4', target='ratio:Fuco')['R2']),
        ('exna_pcr_Tchla_rms', '9.3', 0.09,
         row(exs, kind='abs', model='PCR', target='Tchla')['rms_log']),
        ('exna_pcr_abs_beats_null', '9.3', 7, ex['vs_null']['abs:PCR'].get('beats null', 0)),
        ('exna_null_Fuco_bias', '9.4', -0.62,
         row(exs, kind='abs', model='null_GSM', target='Fuco')['bias_log']),
        ('exna_W_ratio_cov68', '9.3', 0.28, ex['coverage_sharedW']['ratio']['pooled']['cov68']),
        ('exna_W_ratio_cov95', '9.3', 0.89, ex['coverage_sharedW']['ratio']['pooled']['cov95']),
        ('exna_pcr_Fuco_dRMS_hi95', '9.4', 'negative',
         'negative' if row(exv, kind='abs', target='Fuco', model='PCR')['hi95'] < 0
         else 'non-negative'),
    ]
    res = []
    for cid, sec, printed, value in out:
        if isinstance(printed, str):
            ok = printed == value
        else:
            dec = len(str(printed).split('.')[1]) if '.' in str(printed) else 0
            ok = abs(float(value) - printed) <= 0.5 * 10**-dec + 1e-12 if dec else \
                round(float(value)) == printed
        res.append({'id': cid, 'section': sec, 'printed': printed,
                    'from_outputs': value if isinstance(value, str) else float(value),
                    'ok': bool(ok)})
    return res


def main():
    text = REPORT.read_text()
    git = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True,
                         text=True).stdout.strip()
    sections, problems = {}, []
    claimed = set()
    for sec, (scripts, outputs) in SECTIONS.items():
        srec = []
        for s in scripts:
            p = REPO / 'scripts' / 'sdp' / s
            srec.append({'script': f'scripts/sdp/{s}', 'exists': p.exists(),
                         'sha256': sha256(p) if p.exists() else None,
                         'cited_in_report': f'scripts/sdp/{s}' in text})
            if not p.exists():
                problems.append(f'missing script {s}')
        orec = []
        for o in outputs:
            claimed.add(o)
            p = FIG / o
            r = {'output': f'reports/figures/sdp/{o}', 'exists': p.exists()}
            if p.exists():
                r['sha256'] = sha256(p)
                if o.endswith('.json'):
                    r['records_script'] = json.loads(p.read_text()).get('script')
                    if r['records_script'] not in [f'scripts/sdp/{s}' for s in scripts]:
                        problems.append(f'{o} records script {r["records_script"]}')
            else:
                problems.append(f'missing output {o}')
            orec.append(r)
        sections[sec] = {'scripts': srec, 'outputs': orec}

    referenced = sorted(set(re.findall(r'figures/sdp/([A-Za-z0-9_.\-]+)', text)))
    for f in referenced:
        if not (FIG / f).exists():
            problems.append(f'report references missing {f}')
    for p in sorted((REPO / 'scripts' / 'sdp').glob('*.py')):
        if f'scripts/sdp/{p.name}' not in text:
            problems.append(f'script not cited in report: {p.name}')
    unclaimed = sorted(p.name for p in FIG.iterdir()
                       if p.is_file() and p.name not in claimed
                       and p.name != 'report_manifest.json')
    problems += [f'output not mapped to a section: {u}' for u in unclaimed]

    cl = claims()
    problems += [f'claim {c["id"]}: printed {c["printed"]}, outputs give {c["from_outputs"]}'
                 for c in cl if not c['ok']]

    manifest = {'script': 'scripts/sdp/report_provenance.py', 'git_head': git,
                'report_sha256': sha256(REPORT), 'sections': sections,
                'figures_referenced': referenced, 'claims': cl, 'problems': problems,
                'run_utc': datetime.datetime.now(datetime.timezone.utc)
                .isoformat(timespec='seconds')}
    (FIG / 'report_manifest.json').write_text(json.dumps(manifest, indent=1, default=str))
    for c in cl:
        print(f"{'ok ' if c['ok'] else 'BAD'} §{c['section']:<4} {c['id']:<34} "
              f"printed {c['printed']!s:<9} outputs {c['from_outputs']}")
    print(f'{len(referenced)} figures referenced; {len(problems)} problems')
    for p in problems:
        print('  -', p)
    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
