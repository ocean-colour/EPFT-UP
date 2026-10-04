"""Execution #11: slide deck for the SDP report.

Builds ``reports/slides/SDP_Claude_Slides.pptx`` (16:9) from:

* figures already made by the ``scripts/sdp/`` pipeline
  (``reports/figures/sdp/*.png``);
* three summary charts drawn here from the saved outputs (no new
  calculations): Table 2 ours vs paper, LOCO skill by predictor space, and
  LOCO skill by model. They go to ``reports/slides/figures/``;
* figures from the Kramer papers in ``context/papers/``, extracted from the
  PDFs at run time into ``$OS_COLOR/PANGAEA/Kramer2022/paper_figs/`` (not
  committed). They are embedded in the deck with a citation on each slide.

House rule from the prompt: **no text smaller than 20 pt.** Every text run
set here is ≥ 20 pt (checked at the end). Text *inside* figure images is not
controllable here.

A provenance sidecar ``reports/slides/SDP_Claude_Slides.json`` records the
inputs (with SHA-256) and the git HEAD.

Run with::

    conda run -n ocean14 python scripts/sdp/make_slides.py
"""
import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

import fitz  # PyMuPDF
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from PIL import Image  # noqa: E402
from pptx import Presentation  # noqa: E402
from pptx.dml.color import RGBColor  # noqa: E402
from pptx.enum.text import PP_ALIGN  # noqa: E402
from pptx.util import Inches, Pt  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FIG = REPO / 'reports' / 'figures' / 'sdp'
SLIDES = REPO / 'reports' / 'slides'
SFIG = SLIDES / 'figures'
PAPERS = REPO / 'context' / 'papers'
PFIG = Path(os.environ['OS_COLOR']) / 'PANGAEA' / 'Kramer2022' / 'paper_figs'
OUT = SLIDES / 'SDP_Claude_Slides.pptx'

MIN_PT = 20
TITLE_PT, BODY_PT, CAP_PT = 34, 22, 20
NAVY = RGBColor(0x1F, 0x3A, 0x5F)
GREY = RGBColor(0x55, 0x55, 0x55)
ACCENT = RGBColor(0xB0, 0x3A, 0x2E)
W, H = 13.333, 7.5

#: (key, pdf, page (1-based), image width, image height, citation)
PAPER_FIGS = [
    ('k22_fig1', 'kramer2022.pdf', 4, 1889, 860, 'Kramer et al. 2022, Fig. 1'),
    ('k22_fig2', 'kramer2022.pdf', 4, 2156, 1153, 'Kramer et al. 2022, Fig. 2'),
    ('k22_fig3', 'kramer2022.pdf', 5, 1778, 2114, 'Kramer et al. 2022, Fig. 3'),
    ('k22_fig5', 'kramer2022.pdf', 7, 2156, 1211, 'Kramer et al. 2022, Fig. 5'),
    ('k22_fig6', 'kramer2022.pdf', 9, 1660, 2510, 'Kramer et al. 2022, Fig. 6'),
    ('k24_fig1', 'Kramer_etal_2024.pdf', 4, 1145, 766, 'Kramer et al. 2024, Fig. 1'),
    ('k24_fig3', 'Kramer_etal_2024.pdf', 6, 1145, 682, 'Kramer et al. 2024, Fig. 3'),
    ('ks19_dendro', 'Kramer2019.pdf', 9, 2610, 1733, 'Kramer & Siegel 2019, Fig. 3'),
    ('ks19_eof', 'Kramer2019.pdf', 10, 3661, 2426, 'Kramer & Siegel 2019, Fig. 4'),
]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ----------------------------------------------------------------- figures
def extract_paper_figs():
    """Largest-matching embedded image per (pdf, page, size) → PNG."""
    PFIG.mkdir(parents=True, exist_ok=True)
    out = {}
    for key, pdf, page, w, h, cite in PAPER_FIGS:
        doc = fitz.open(PAPERS / pdf)
        hits = [im for im in doc[page - 1].get_images(full=True)
                if im[2] == w and im[3] == h]
        if not hits:
            raise SystemExit(f'{pdf} p{page}: no {w}x{h} image')
        pix = fitz.Pixmap(doc, hits[0][0])
        if pix.n - pix.alpha >= 4:  # CMYK
            pix = fitz.Pixmap(fitz.csRGB, pix)
        dest = PFIG / f'{key}.png'
        pix.save(dest)
        out[key] = {'path': dest, 'cite': cite, 'pdf': pdf, 'page': page,
                    'pdf_sha256': sha256(PAPERS / pdf)}
    return out


def chart_style():
    plt.rcParams.update({'font.size': 16, 'axes.titlesize': 18, 'axes.labelsize': 17,
                         'xtick.labelsize': 15, 'ytick.labelsize': 15,
                         'legend.fontsize': 15})


def chart_table2():
    t = json.loads((FIG / 'pcr_summary.json').read_text())['table2']['n145']
    pigs = sorted(t, key=lambda p: -t[p]['R2_paper'])
    x = np.arange(len(pigs))
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.errorbar(x - 0.12, [t[p]['R2_paper'] for p in pigs],
                yerr=[t[p]['R2_sd_paper'] for p in pigs], fmt='o', ms=9, capsize=4,
                color='0.4', label='Kramer et al. 2022, Table 2')
    ax.errorbar(x + 0.12, [t[p]['R2'] for p in pigs], yerr=[t[p]['R2_sd'] for p in pigs],
                fmt='s', ms=9, capsize=4, color='#b03a2e', label='this work (N = 145)')
    ax.set_xticks(x, pigs, rotation=45, ha='right')
    ax.set_ylabel('validation R² (mean ± SD, 100 splits)')
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc='upper right')
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    p = SFIG / 'slide_table2.png'
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def chart_spaces():
    s = pd.read_csv(FIG / 'source_spaces_summary.csv').set_index('model')
    order = [('Rrs', 'Rrs'), ('Rrs_d2', "Rrs''"), ('dRrs', 'δRrs'),
             ('dRrs_d2', "δRrs'' (paper)"), ('M1', 'M1 spline resid.'),
             ('M3_w7', 'M3 SG 7 nm'), ('M3_w11', 'M3 SG 11 nm'), ('M3_w21', 'M3 SG 21 nm')]
    lab = [o[1] for o in order]
    y = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(10, 5.2))
    a = [s.loc[o[0], 'abs_R2_loco'] for o in order]
    r = [s.loc[o[0], 'ratio_logR2_loco'] for o in order]
    d2 = ["''" in o[1] or 'SG 7' in o[1] for o in order]
    ax.barh(y + 0.2, a, 0.4, color=['#b03a2e' if d else '#1f3a5f' for d in d2],
            label='absolute pigments (R²)')
    ax.barh(y - 0.2, r, 0.4, color=['#e8a49c' if d else '#8fa8c8' for d in d2],
            label='log pigment:Tchla (log-R²)')
    ax.set_yticks(y, lab)
    ax.invert_yaxis()
    ax.set_xlabel('mean skill, leave-one-campaign-out (PCR in every space)')
    ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.grid(axis='x', alpha=0.3)
    fig.tight_layout()
    p = SFIG / 'slide_spaces.png'
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


def chart_models():
    s = pd.read_csv(FIG / 'weighting_summary.csv').set_index('model')
    order = [('PCR_dRrs_d2', "PCR δRrs'' (paper)"), ('ridge_dRrs_d2', "ridge δRrs''"),
             ('PLS_dRrs_d2', "PLS δRrs''"), ('PLS_dRrs_aux', 'PLS δRrs ⊕ GSM'),
             ('smooth_dRrs_aux_group', 'smooth, per pigment'),
             ('RRR_dRrs_aux_group', 'shared W (group CV)'),
             ('RRR_M1_aux_group', 'shared W on M1')]
    y = np.arange(len(order))
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True,
                            gridspec_kw={'width_ratios': [2, 1]})
    ax = axs[0]
    ax.barh(y + 0.2, [s.loc[o[0], 'abs_R2_loco'] for o in order], 0.4, color='#1f3a5f',
            label='absolute (R²)')
    ax.barh(y - 0.2, [s.loc[o[0], 'ratio_logR2_loco'] for o in order], 0.4,
            color='#8fa8c8', label='ratios (log-R²)')
    ax.set_yticks(y, [o[1] for o in order])
    ax.invert_yaxis()
    ax.set_xlabel('mean LOCO skill')
    ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.grid(axis='x', alpha=0.3)
    ax = axs[1]
    ax.barh(y, [s.loc[o[0], 'ratio_beats_null_loco'] for o in order], 0.6, color='#b03a2e')
    ax.set_xlabel('ratios beating\nTchla null (LOCO)')
    ax.set_xlim(0, 12)
    ax.grid(axis='x', alpha=0.3)
    fig.tight_layout()
    p = SFIG / 'slide_models.png'
    fig.savefig(p, dpi=200)
    plt.close(fig)
    return p


# ----------------------------------------------------------------- pptx helpers
def _runs(par, text, size, bold=False, color=None):
    """Plain text with **bold** spans."""
    parts = text.split('**')
    for i, t in enumerate(parts):
        if not t:
            continue
        r = par.add_run()
        r.text = t
        r.font.size = Pt(size)
        r.font.bold = bold or (i % 2 == 1)
        if color is not None:
            r.font.color.rgb = color


def add_title(slide, text, sub=None):
    tb = slide.shapes.add_textbox(Inches(0.5), Inches(0.25), Inches(W - 1), Inches(0.9))
    tf = tb.text_frame
    tf.word_wrap = True
    _runs(tf.paragraphs[0], text, TITLE_PT, bold=True, color=NAVY)
    if sub:
        p = tf.add_paragraph()
        _runs(p, sub, BODY_PT, color=GREY)
    line = slide.shapes.add_shape(1, Inches(0.5), Inches(1.15 if not sub else 1.55),
                                  Inches(W - 1), Inches(0.04))
    line.fill.solid()
    line.fill.fore_color.rgb = ACCENT
    line.line.fill.background()


def add_bullets(slide, items, x, y, w, h, size=BODY_PT):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for i, it in enumerate(items):
        level = 1 if it.startswith('- ') else 0
        txt = it[2:] if level else it
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = level
        p.space_after = Pt(8)
        _runs(p, ('– ' if level else '• ') + txt, size)
    return tb


def add_image(slide, path, x, y, w, h):
    """Fit inside the box, keep aspect, centre."""
    iw, ih = Image.open(path).size
    s = min(w / iw, h / ih)
    pw, ph = iw * s, ih * s
    slide.shapes.add_picture(str(path), Inches(x + (w - pw) / 2), Inches(y + (h - ph) / 2),
                             Inches(pw), Inches(ph))


def add_caption(slide, text, y=6.95):
    tb = slide.shapes.add_textbox(Inches(0.5), Inches(y), Inches(W - 1), Inches(0.45))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.RIGHT
    _runs(p, text, CAP_PT, color=GREY)


class Deck:
    def __init__(self):
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = Inches(W), Inches(H)
        self.blank = self.prs.slide_layouts[6]

    def new(self, title, sub=None):
        s = self.prs.slides.add_slide(self.blank)
        add_title(s, title, sub)
        return s

    def fig_text(self, title, img, bullets, cap=None, img_frac=0.6, sub=None):
        """Figure left, bullets right."""
        s = self.new(title, sub)
        top = 1.35 if not sub else 1.75
        iw = (W - 1.0) * img_frac
        add_image(s, img, 0.4, top, iw, 6.85 - top)
        add_bullets(s, bullets, 0.4 + iw + 0.25, top, W - iw - 1.05, 6.85 - top)
        if cap:
            add_caption(s, cap)
        return s

    def fig_full(self, title, img, bullets=None, cap=None, img_h=4.3):
        """Figure on top, bullets below."""
        s = self.new(title)
        add_image(s, img, 0.4, 1.3, W - 0.8, img_h)
        if bullets:
            add_bullets(s, bullets, 0.6, 1.35 + img_h, W - 1.2, 6.9 - 1.35 - img_h)
        if cap:
            add_caption(s, cap)
        return s

    def two_figs(self, title, img1, img2, bullets=None, cap=None, img_h=4.2):
        s = self.new(title)
        half = (W - 1.0) / 2
        add_image(s, img1, 0.4, 1.3, half, img_h)
        add_image(s, img2, 0.6 + half, 1.3, half, img_h)
        if bullets:
            add_bullets(s, bullets, 0.6, 1.4 + img_h, W - 1.2, 6.9 - 1.4 - img_h)
        if cap:
            add_caption(s, cap)
        return s

    def text(self, title, bullets, sub=None, size=BODY_PT):
        s = self.new(title, sub)
        add_bullets(s, bullets, 0.7, 1.45 if not sub else 1.85, W - 1.4, 5.4, size)
        return s


# ----------------------------------------------------------------- deck
def build(pf, charts):
    d = Deck()
    F = lambda n: FIG / n  # noqa: E731

    # 1 title
    s = d.prs.slides.add_slide(d.blank)
    tb = s.shapes.add_textbox(Inches(0.8), Inches(2.0), Inches(W - 1.6), Inches(3.5))
    tf = tb.text_frame
    tf.word_wrap = True
    _runs(tf.paragraphs[0], 'The 2nd-derivative (SDP) approach to phytoplankton '
          'pigments from hyperspectral Rrs', 40, bold=True, color=NAVY)
    p = tf.add_paragraph()
    _runs(p, 'Reproduction, critique and alternatives, with uncertainty', 28, color=ACCENT)
    p = tf.add_paragraph()
    p.space_before = Pt(30)
    _runs(p, 'Claude (Opus 5.5) for J. Xavier Prochaska · EPFT-UP · October 2026', 22,
          color=GREY)
    p = tf.add_paragraph()
    _runs(p, 'Report: reports/SDP_Claude_Report.md', 20, color=GREY)

    # --- WHY
    d.fig_text('Why: PFTs from hyperspectral ocean color',
               pf['k22_fig1']['path'],
               ['Phytoplankton functional types (PFTs) are what ecosystem and '
                'carbon models need',
                'Multispectral color gives chlorophyll, not composition',
                'PACE OCI now delivers global hyperspectral Rrs',
                '**SDP** (Kramer et al. 2022) is the leading derivative method; '
                'NASA ships an OCI notebook (oci_sdp)',
                'Training set: 145 global HPLC + Rrs matchups, 8 campaigns'],
               cap=pf['k22_fig1']['cite'] + ' (training matchups)', img_frac=0.55)

    # --- WHAT
    d.fig_text('What SDP does',
               pf['k22_fig2']['path'],
               ['Fit a GSM-like model to Rrs (Tchla, a_dg, b_bp)',
                'Residual **δRrs = Rrs − Rrs_mod**',
                'Second derivative δRrs″ (1 nm)',
                'PCR: z-score, PCs by inner CV, 100 random 75/25 splits',
                'Output: 13 HPLC pigment concentrations'],
               cap=pf['k22_fig2']['cite'] + ': (A) Rrs, (B) Rrs_mod, (C) δRrs', img_frac=0.6)

    d.fig_text('Why derivatives? Pigment features in δRrs',
               pf['k22_fig5']['path'],
               ['Pigment–spectrum correlations sharpen from δRrs to δRrs′ to δRrs″',
                'Lineage: Catlett & Siegel 2018 used a_ph′, a_ph″ with PCR',
                'Their argument: pigment **covariance** makes communities retrievable'],
               cap=pf['k22_fig5']['cite'], img_frac=0.6)

    d.two_figs('What HPLC itself can resolve',
               pf['ks19_dendro']['path'], pf['ks19_eof']['path'],
               ['Globally **4** pigment groups (diatoms + dinos merge); 4–6 locally',
                'Leading EOF: diatoms/dinos vs picophytoplankton (24%)'],
               cap='Kramer & Siegel 2019 (4,480 samples), Figs. 3–4')

    d.fig_text('What SDP claims (2022)',
               pf['k22_fig6']['path'],
               ['Validation R² **0.37–0.72** for 13 pigments (Table 2)',
                '1–5 nm sampling is enough; 10 nm degrades',
                'The 5 HPLC pigment groups re-emerge from modelled pigments',
                'Validated on random 75/25 splits; no chlorophyll-only baseline'],
               cap=pf['k22_fig6']['cite'], img_frac=0.5)

    d.two_figs('2024: δRrs defines optical communities',
               pf['k24_fig1']['path'], pf['k24_fig3']['path'],
               ['145 + 17 EXPORTS-NA samples; 3 δRrs communities vs 3 HPLC communities',
                '**74%** of samples assigned to matching communities; same at 5 nm'],
               cap='Kramer et al. 2024, Fig. 1 (left) and Fig. 3 (right)')

    # --- HOW
    d.text('How we tested it', [
        '**Data:** PANGAEA 937536 (145 samples), SHA-256 pinned; Kramer\'s MATLAB '
        'code and the Python port at pinned commits',
        '**Re-derived independently** in Python (epft_up/sdp/), then checked '
        'against both codes',
        '**Two tests for every model:**',
        '- skill above a 2-parameter **Tchla-only null** (pigment = a·Tchla^b)',
        '- skill under **leave-one-campaign-out** (LOCO), not only random splits',
        '**Then:** maths of the derivative, alternatives, Bayesian uncertainty, '
        'PACE-like noise, an independent campaign (EXPORTS-NA 2021)',
        'Every number traces to a script (13 scripts, 100 tests)'])

    # --- REPRODUCTION
    d.fig_text('Reproduction 1: the GSM residual',
               F('kramer_fig4_chl.png'),
               ['OC4: y = 0.873x − 0.138, R² 0.746 (paper 0.87, 0.75)',
                'GSM: slope **0.961**, R² **0.864**, as the paper, once one '
                'degenerate SABOR fit is dropped',
                'Python port δRrs agrees below the deposit\'s rounding',
                'Paper typos fixed from the code (S_dg sign, η ratio)'],
               cap='this work, §4.1', img_frac=0.6)

    d.fig_text('Reproduction 2: Table 2',
               charts['table2'],
               ['**All 13** R² within one quoted SD of the paper',
                'Kramer\'s own trained coefficients reproduce on our δRrs″',
                'Open: our Fig. 6 slopes are 0.02–0.23 lower (same R²)'],
               cap='this work, §4.2 (pcr_summary.json)', img_frac=0.66)

    d.fig_full('Reproduction 3: Fig. 6 and the pigment groups',
               F('kramer_fig6_pigments.png'),
               ['Measured ratios give the paper\'s 5 groups exactly; modelled ones only '
                'loosely, and depend on how zeros are handled'],
               cap='this work, §4.2', img_h=4.6)

    # --- MATHS
    d.fig_text('The derivative is a fixed linear operator',
               F('maths_effective_weights.png'),
               ['p̂ = Aᵀ(Dx) = (**DᵀA**)ᵀx: a PCR on δRrs″ **is** a linear model on δRrs',
                'It cannot add information; it only changes the prior',
                '**99.7%** of the effective weight power is at periods < 10 nm',
                'Absorption features are 10–40 nm wide'],
               cap='this work, §5.1', img_frac=0.6)

    d.fig_full('… and it implies an anti-smooth prior',
               F('maths_prior.png'),
               ['Derivative + z-score + shrinkage ⇒ Cov(w) ∝ DᵀS⁻²D: **≈5300×** more prior '
                'power at high than at low frequency',
                'Largest weights forced into the red, where the signal is weakest'],
               cap='this work, §5.1', img_h=4.0)

    d.fig_text('Noise: information content and amplification',
               F('maths_dof.png'),
               ['Degrees of freedom for signal: ~**10** (2% in situ), ~**6** (PACE white), '
                'not 300',
                '1 → 5 nm costs ≈ 25%; noise matters more than resolution',
                'PACE white noise: SDP Tchla noise = **11.6×** natural SD of Tchla '
                '(0.83× for weights fit on δRrs)'],
               cap='this work, §5.2', img_frac=0.58)

    # --- DIAGNOSTICS
    d.fig_text('Skill beyond chlorophyll?',
               F('diag_skill_absolute.png'),
               ['A power law in **GSM Tchla** beats SDP for **11 of 13** pigments',
                'Only DVchla and Zea (cyanobacteria) beat it, in random splits',
                'Ratios implied by SDP\'s products: R² 0.02–0.09'],
               cap='this work, §6.1', img_frac=0.6)

    d.fig_text('Leakage: hold out a whole campaign',
               F('diag_loco_tchla.png'),
               ['SDP Tchla R² **0.72 → 0.44** under LOCO; null 0.66',
                'Roughly half of SDP\'s skill is within-cruise similarity',
                'No ratio beats the null significantly out of campaign'],
               cap='this work, §6.2', img_frac=0.55)

    d.fig_text('What the residual buys — and what the derivative costs',
               charts['spaces'],
               ['Red bars: finite-difference spaces — **worst** out of campaign',
                'Same information without differentiating: ~2× better',
                'Empirical spline residual (M1) ≈ GSM residual: the NOMAD baseline '
                'is not what helps',
                'Wider Savitzky–Golay windows transfer better'],
               cap='this work, §6.3 (source_spaces_summary.csv)', img_frac=0.6)

    # --- ALTERNATIVES
    d.fig_text('Alternative: a shared, smooth, noise-aware weighting',
               charts['models'],
               ['One low-rank W for all pigments on δRrs ⊕ GSM parameters',
                'Penalty: noise variance + smoothness; rank by leave-campaign-out',
                'Changing only the estimator (ridge, PLS) does **nothing**',
                'Sharing across pigments lifts the ratios: **6** beat the null '
                'out of campaign (PCR: 1)'],
               cap='this work, §7 (weighting_summary.csv)', img_frac=0.6)

    d.two_figs('Only two compositional axes transfer',
               F('weighting_latent_spectra.png'), pf['ks19_eof']['path'],
               ['Latent 1 ≈ diatom ↔ cyanobacteria = Kramer & Siegel EOF mode 1',
                'Latent 2: green algae / haptophytes. Weight sits at 10–40 nm scales'],
               cap='this work, §7.3 (left); Kramer & Siegel 2019 Fig. 4 (right)')

    d.fig_text('Which wavelengths carry the weight',
               F('weighting_target_spectra.png'),
               ['Shared W: **no** power at periods < 10 nm (PCR: 99.7%)',
                '≈ 46% of the weight at **600–700 nm** for Fuco, Zea, DVchla, Chlc12',
                'Lobes at the Chl red peak / 683 nm fluorescence and Chl c, '
                'phycobilin bands (suggestive)'],
               cap='this work, §7.4', img_frac=0.6)

    # --- UNCERTAINTY
    d.fig_full('Calibrated intervals out of campaign',
               F('uncertainty_coverage.png'),
               ['Bayesian shared W, cross-fitted, Monte Carlo noise through the GSM: '
                '68/95% coverage **0.71/0.94**, flat PIT (report §8.2)',
                'Without input-noise propagation, coverage drops under PACE noise'],
               cap='this work, §8.2', img_h=3.9)

    d.fig_full('Degradation: in situ → PACE noise, 1 → 5 nm',
               F('uncertainty_skill.png'),
               ['SDP\'s PCR: R² ≈ **0** under any realistic noise',
                'Noise-aware shared W, PACE white: Fuco:Tchla R² **0.36** (null 0.27), '
                '±0.33 dex; **5 nm sampling costs nothing**'],
               cap='this work, §8.3', img_h=3.9)

    # --- HOLD-OUT
    d.two_figs('Independent campaign: EXPORTS-NA 2021 (17 samples)',
               F('exports_na_tchla.png'), F('exports_na_ratios.png'),
               ['SDP gets the bloom\'s diatom pigment **levels** right (nulls biased '
                'low by 0.6 dex), but over-predicts trace pigments by 0.5–1.1 dex',
                'One campaign = one draw from the LOCO spread; 68% intervals miss the '
                'campaign bias'],
               cap='this work, §9', img_h=3.9)

    # --- CONCLUSIONS
    d.text('What we found', [
        'SDP **reproduces** from public data and code',
        'The 2nd derivative is a **prior, not information**: anti-smooth, noise '
        'amplifying (≥10×)',
        'Most of SDP\'s skill is **chlorophyll** and **within-campaign** similarity',
        'Hyperspectral Rrs does carry composition: ~**two** smooth (10–40 nm) axes '
        'transfer between campaigns',
        'A shared, smooth, noise-aware W with calibrated intervals beats the null for '
        'Fuco, Zea, DVchla, Chlc12 ratios',
        'At PACE noise: SDP has no skill; the aware W keeps a modest edge'], size=24)

    d.text('Recommendations', [
        '**Always report** skill against a Tchla-only null **and** a held-out campaign',
        '**Target** log pigment:Tchla ratios; shared low-rank W on δRrs (or M1)',
        '**Add** a between-campaign variance term to the intervals',
        '**Data:** raw (unsmoothed) spectra + uncertainties; more campaigns '
        '(El Hourany 237 stations, PACE-era matchups)',
        '**Fix** the a_w smoothing mismatch in the forward model',
        '**Concentrations:** censored log-normal likelihood (Tom Jordan)'], size=24)

    d.text('What PACE OCI application would require', [
        'An OCI noise model **split into white and smooth parts** (Rrs_unc gives '
        'only the diagonal)',
        'Training spectra convolved to OCI bands, **not pre-smoothed**; same for '
        'a_w, A, B',
        'Handle fluorescence (683 nm), Raman and gas bands that δRrs keeps',
        'Validate on PACE-era matchups against both baselines, and against '
        'NASA\'s SDP product',
        'Per-pixel intervals + campaign term; stated averaging footprint',
        'SDP coefficients at OCI white noise: ≥ **135** pixels averaged just to '
        'reach one natural SD'], size=24)

    d.text('Provenance and references', [
        'Code: epft_up/sdp/; scripts: scripts/sdp/ (run order in report §10.4)',
        'Checks: scripts/sdp/report_provenance.py — 34/34 headline numbers verified',
        'This deck: scripts/sdp/make_slides.py',
        'Kramer, Siegel, Maritorena & Catlett 2022, RSE 270, 112879',
        'Kramer et al. 2024, Opt. Express 32, 34482',
        'Kramer & Siegel 2019, JGR Oceans; Catlett & Siegel 2018, JGR Oceans',
        'El Hourany & Kramer 2026 (preprint); Lange et al. 2020, Opt. Express'],
        size=20)
    return d.prs


def check_fonts(prs):
    bad = []
    for i, sl in enumerate(prs.slides, 1):
        for sh in sl.shapes:
            if not sh.has_text_frame:
                continue
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if r.font.size is None or r.font.size.pt < MIN_PT:
                        bad.append((i, r.text[:30], r.font.size))
    return bad


def main():
    SFIG.mkdir(parents=True, exist_ok=True)
    chart_style()
    pf = extract_paper_figs()
    charts = {'table2': chart_table2(), 'spaces': chart_spaces(), 'models': chart_models()}
    prs = build(pf, charts)
    bad = check_fonts(prs)
    if bad:
        raise SystemExit(f'runs below {MIN_PT} pt: {bad}')
    prs.save(OUT)
    git = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True,
                         text=True).stdout.strip()
    used = sorted({str(p.relative_to(REPO)) for p in FIG.glob('*.png')})
    prov = {'script': 'scripts/sdp/make_slides.py', 'git_head': git,
            'n_slides': len(prs.slides), 'min_font_pt': MIN_PT,
            'paper_figures': {k: {kk: str(vv) for kk, vv in v.items()} for k, v in pf.items()},
            'charts': {k: {'path': str(v.relative_to(REPO)), 'sha256': sha256(v)}
                       for k, v in charts.items()},
            'pipeline_figures_available': used,
            'output': {'file': str(OUT.relative_to(REPO)), 'sha256': sha256(OUT)},
            'run_utc': datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec='seconds')}
    OUT.with_suffix('.json').write_text(json.dumps(prov, indent=1))
    print(f'wrote {OUT} ({len(prs.slides)} slides)')


if __name__ == '__main__':
    main()
