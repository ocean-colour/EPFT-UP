"""Generate the figures for ``reports/MOANA_Claude_Report.md``.

Two groups of figures:

* The PCA-basis figure is derived purely from the vendored MOANA reference tables
  in ``epft_up/data/moana/``, so it runs anywhere the repo is checked out.
* The product figures need one PACE MOANA granule under ``$OS_COLOR/PACE/``. They
  are skipped with a message if it is absent, mirroring the repo's convention of
  letting data-dependent work opt out rather than fail.

Run with::

    python reports/scripts/moana_report_figs.py

Outputs (written to ``reports/figures/``):
    moana_pc_loadings.png     -- loadings of the principal components MOANA uses,
                                 annotated with which taxa use each one.
    moana_product_masks.png   -- maps of the three products with the land-254
                                 sentinel and the unflagged clipped zeros exposed.
    moana_clipping.png        -- how often each product is clipped to zero, and the
                                 abundance distributions of what survives.
    moana_mapping_verdict.png -- which PC mapping the shipping product uses.
    moana_mapping_consequence.png -- what the misplaced coefficients cost,
                                 operational vs as-published (section 7.1).
    moana_heldout_skill.png   -- retrieved vs observed on Lange's held-out
                                 cruises AMT23/25/28 (section 12.3).

Colours follow the house data-viz palette: a single-hue blue sequential ramp for
magnitude, the first three categorical slots for taxon identity (that subset is
validated all-pairs for colour-vision deficiency), the reserved ``critical`` red
for the clipped-zero defect, and recessive greys for chrome.
"""

import json
import os
import sys

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless-safe; no display needed
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

# Repo-relative paths, resolved from this script's location so the working
# directory does not matter.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
DATA_DIR = os.path.join(_REPO, "epft_up", "data", "moana")
FIG_DIR = os.path.join(_REPO, "reports", "figures")
# Run-as-a-script support: sys.path[0] is this file's directory, so the
# repo root (and hence ``epft_up``) is not importable without this.
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

# One PACE MOANA granule (prompt-3 original location, then the prompt-14
# earthaccess cache) plus its matching L3M AOP Rrs input for the mapping-
# verdict figure. First existing path wins.
def _first_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return paths[-1]


_PACE = os.path.join(os.environ.get("OS_COLOR", ""), "PACE")
GRANULE = _first_existing(
    os.path.join(_PACE, "PACE_OCI.20250701.L4m.DAY.MOANA.V3_2.0p1deg.nc"),
    os.path.join(_PACE, "moana_validation",
                 "PACE_OCI.20250701.L4m.DAY.MOANA.V3_2.0p1deg.nc"),
)
AOP_GRANULE = _first_existing(
    os.path.join(_PACE, "moana_validation",
                 "PACE_OCI.20250701.L3m.DAY.AOP.V3_2.0p1deg.nc"),
)

# --- palette (house data-viz reference instance, light surface) ---------------
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"
# Categorical slots 1-3: identity of the three taxa. This subset validates
# all-pairs in both light and dark modes, so it is safe for map/scatter forms.
TAXON_COLOR = {"pro": "#2a78d6", "syn": "#eb6834", "apeuk": "#1baf7a"}
# Reserved status colour, used ONLY for the clipped-zero defect -- never as a series.
CRITICAL = "#d03b3b"
LAND_GREY = "#c3c2b7"
# Single-hue sequential ramp (blue 100 -> 700) for continuous magnitude.
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_CMAP = LinearSegmentedColormap.from_list("moana_blue", BLUE_RAMP)


# Diverging ramp for SIGNED differences (blue <-> red, neutral grey midpoint),
# per the house palette's diverging pair. palette.md tabulates only the blue
# arm, so the red arm is stepped to mirror its lightness progression.
_DIV_BLUE = ["#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4",
             "#cde2fb"]
_DIV_GREY = "#f0efec"
_DIV_RED = ["#f7d4d4", "#efabab", "#e68282", "#e34948", "#c22f2f", "#9b2222",
            "#6e1616"]
DIV_CMAP = LinearSegmentedColormap.from_list(
    "moana_div", _DIV_BLUE + [_DIV_GREY] + _DIV_RED)

#: The two disputed coefficient placements (report §7.1): operational
#: picophyt.json vs ATBD v1.2. Only these two terms move.
DISPUTED_PC = {"syn": ("PC16", "U13"), "pro": ("PC7", "U17")}

# Sentinels used by the MOANA products.
LAND, FILL, I32MIN = 254, -32767, -2147483648

# Display order and labels for the three taxa.
TAXA = [
    ("pro", "prococcus_moana", "Prochlorococcus"),
    ("syn", "syncoccus_moana", "Synechococcus"),
    ("apeuk", "picoeuk_moana", "picoeukaryotes"),
]


def load_moana_tables(data_dir=DATA_DIR):
    """Load the MOANA PCA loadings and regression coefficients.

    Inputs
    ------
    data_dir : str
        Directory holding ``pca_picophyto.h5`` and ``picophyt.json``.

    Outputs
    -------
    wave : ndarray, shape (124,)
        Wavelengths of the algorithm's grid [nm], 414 to 660 in 2 nm steps.
    loadings : ndarray, shape (124, 45)
        Loading matrix ``V[lambda, i]``; columns are orthonormal.
    coefs : dict
        ``{'npc': int, 'pro': list, 'syn': list, 'apeuk': list}`` -- the raw
        zero-padded coefficient vectors exactly as OCSSW reads them.
    """
    with h5py.File(os.path.join(data_dir, "pca_picophyto.h5"), "r") as h5:
        loadings = h5["component"][:]
        wave = h5["wavelength"][:]

    with open(os.path.join(data_dir, "picophyt.json")) as fh:
        raw = json.load(fh)

    coefs = {
        "npc": int(raw["npc"]),
        "pro": [float(x) for x in raw["pro_coef"].split(",")],
        "syn": [float(x) for x in raw["syn_coef"].split(",")],
        "apeuk": [float(x) for x in raw["apeuk_coef"].split(",")],
    }
    return wave, loadings, coefs


def pcs_used_by_taxon(coefs):
    """Map each taxon to the 1-based principal components it actually uses.

    Replicates the indexing in OCSSW ``get_Cpicophyt.c``: for i in 0..npc-1,
    ``pro_coef[i+2]``, ``syn_coef[i+1]`` and ``apeuk_coef[i+1]`` each multiply
    score ``U_(i+1)``. Zero-valued slots mean the PC is unused.

    Inputs
    ------
    coefs : dict
        As returned by :func:`load_moana_tables`.

    Outputs
    -------
    dict
        ``{taxon: {pc_number: coefficient}}`` for the three taxa.
    """
    npc = coefs["npc"]
    # offset of the first PC coefficient within each vector
    offsets = {"pro": 2, "syn": 1, "apeuk": 1}
    used = {}
    for taxon, off in offsets.items():
        vec = coefs[taxon]
        used[taxon] = {
            i + 1: vec[i + off] for i in range(npc) if vec[i + off] != 0.0
        }
    return used


def plot_pc_loadings(wave, loadings, used, n_show=6, outfile=None):
    """Plot the loadings of the leading principal components.

    Inputs
    ------
    wave : ndarray
        Wavelength grid [nm].
    loadings : ndarray, shape (n_wave, n_pc)
        Loading matrix.
    used : dict
        Output of :func:`pcs_used_by_taxon`, used to annotate which taxa
        depend on each component.
    n_show : int
        How many leading components to draw.
    outfile : str or None
        Path to write the PNG to. If None, nothing is written.

    Outputs
    -------
    matplotlib.figure.Figure
        The figure, so a caller can further customise or save it.
    """
    labels = {"pro": "Pro", "syn": "Syn", "apeuk": "picoeuk"}

    fig, axes = plt.subplots(2, 3, figsize=(13, 6.5), sharex=True)
    for k, ax in enumerate(axes.ravel()[:n_show]):
        pc = k + 1
        ax.axhline(0.0, color="0.8", lw=0.8, zorder=0)
        ax.plot(wave, loadings[:, k], color="C0", lw=1.8)

        # Which taxa use this component, and with what sign?
        consumers = [
            f"{labels[t]} ({used[t][pc]:+.3g})" for t in ("pro", "syn", "apeuk")
            if pc in used[t]
        ]
        subtitle = "; ".join(consumers) if consumers else "unused by all three taxa"
        ax.set_title(f"PC{pc}\n{subtitle}", fontsize=9)
        ax.tick_params(labelsize=8)

    for ax in axes[-1]:
        ax.set_xlabel("wavelength [nm]", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("loading", fontsize=9)

    # Note: no '$...$' anywhere in the title -- matplotlib would interpret it as
    # mathtext and mangle the underscores in the filename.
    fig.suptitle(
        "MOANA PCA loadings (pca_picophyto.h5): the fixed basis onto which "
        "every standardised spectrum is projected",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    if outfile is not None:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        fig.savefig(outfile, dpi=150)
        print(f"wrote {outfile}")
    return fig


def _style_axes(ax):
    """Apply the recessive chrome the house style calls for.

    Inputs
    ------
    ax : matplotlib.axes.Axes
        Axes to restyle in place.

    Outputs
    -------
    None
    """
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelsize=8)


def load_granule(path=GRANULE):
    """Read the three MOANA fields and their coordinates from a PACE granule.

    Values are returned exactly as stored (int32 with sentinels intact), because
    the sentinels are the subject of two of the figures.

    Inputs
    ------
    path : str
        Path to a ``*.L4m.*.MOANA.*.nc`` file.

    Outputs
    -------
    fields : dict
        ``{short_key: {'data': ndarray, 'valid_max': int}}`` per taxon.
    lat, lon : ndarray
        Coordinate vectors [degrees].
    """
    import xarray as xr  # imported here so the basis figure needs no xarray

    ds = xr.open_dataset(path, mask_and_scale=False)
    fields = {}
    for key, var, _label in TAXA:
        fields[key] = {
            "data": np.asarray(ds[var].values),
            "valid_max": int(ds[var].attrs["valid_max"]),
        }
    lat = np.asarray(ds["lat"].values)
    lon = np.asarray(ds["lon"].values)
    return fields, lat, lon


def summarize_product(fields):
    """Compute the ocean-only statistics the report quotes.

    Inputs
    ------
    fields : dict
        Output of :func:`load_granule`.

    Outputs
    -------
    dict
        ``{key: {...}}`` with pixel counts, the percentage of real retrievals that
        are exactly zero, the percentage exceeding ``valid_max``, and the positive
        values themselves (for the distribution panel).
    """
    stats = {}
    for key, _var, label in TAXA:
        d = fields[key]["data"]
        vmax = fields[key]["valid_max"]
        # "Real" = anything that is not land, not fill, and not the int32 underflow.
        ocean = (d != LAND) & (d != FILL) & (d != I32MIN)
        vals = d[ocean].astype("int64")
        stats[key] = {
            "label": label,
            "n_ocean": vals.size,
            "n_land": int((d == LAND).sum()),
            "n_fill": int((d == FILL).sum()),
            "n_i32min": int((d == I32MIN).sum()),
            "pct_zero": 100.0 * (vals == 0).mean(),
            "pct_over": 100.0 * (vals > vmax).mean(),
            "valid_max": vmax,
            "positive": vals[vals > 0],
        }
    return stats


def plot_product_masks(fields, lat, lon, outfile=None):
    """Map each product with the land sentinel and the clipped zeros exposed.

    The point of the figure: land is stored as 254 cells/mL, *inside* the declared
    valid range, and pixels clipped to zero carry no flag. Both are drawn in their
    own colours so the reader can see how much of the field they occupy.

    Inputs
    ------
    fields, lat, lon
        As returned by :func:`load_granule`.
    outfile : str or None
        Path to write the PNG to.

    Outputs
    -------
    matplotlib.figure.Figure
    """
    # lat runs north -> south in the file, which is what origin='upper' expects.
    assert lat[0] > lat[-1], "expected descending latitudes"
    extent = [lon.min() - 0.05, lon.max() + 0.05, lat.min() - 0.05, lat.max() + 0.05]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.6), facecolor=SURFACE)
    for ax, (key, _var, label) in zip(axes, TAXA):
        d = fields[key]["data"]
        ocean = (d != LAND) & (d != FILL) & (d != I32MIN)

        # Layer 1: magnitude, log10 of the positive retrievals only.
        pos = ocean & (d > 0)
        mag = np.full(d.shape, np.nan)
        mag[pos] = np.log10(d[pos].astype("float64"))
        im = ax.imshow(mag, extent=extent, origin="upper", cmap=SEQ_CMAP,
                       interpolation="nearest")

        # Layer 2: land, recessive grey.
        ax.imshow(np.where(d == LAND, 1.0, np.nan), extent=extent, origin="upper",
                  cmap=LinearSegmentedColormap.from_list("l", [LAND_GREY, LAND_GREY]),
                  interpolation="nearest", vmin=0, vmax=1)

        # Layer 3: the defect -- retrievals clipped to exactly zero, drawn last so
        # they are never hidden underneath anything else.
        ax.imshow(np.where(ocean & (d == 0), 1.0, np.nan), extent=extent,
                  origin="upper",
                  cmap=LinearSegmentedColormap.from_list("c", [CRITICAL, CRITICAL]),
                  interpolation="nearest", vmin=0, vmax=1)

        pct = 100.0 * (d[ocean] == 0).mean()
        ax.set_title(f"{label}\n{pct:.1f}% of retrievals clipped to zero",
                     fontsize=9.5, color=INK_PRIMARY)
        ax.set_facecolor(SURFACE)
        _style_axes(ax)
        # No axis labels: the degree ticks are self-explanatory on a map, and a
        # "longitude" label collides with the horizontal colorbar below it.

        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.06, fraction=0.045)
        cb.set_label("log$_{10}$(cells mL$^{-1}$)", fontsize=8, color=INK_MUTED)
        cb.ax.tick_params(labelsize=7, colors=INK_MUTED)
        cb.outline.set_edgecolor(BASELINE)

    # (no y-axis label either, for the same reason as the x-axis)

    # Legend names every non-data colour, so nothing is carried by colour alone.
    fig.legend(
        handles=[
            Patch(facecolor=CRITICAL, label="clipped to 0 — unflagged"),
            Patch(facecolor=LAND_GREY, label="land, stored as 254 cells mL$^{-1}$"),
            Patch(facecolor=SURFACE, edgecolor=BASELINE, label="no retrieval (fill)"),
        ],
        loc="lower center", ncol=3, frameon=False, fontsize=9,
        bbox_to_anchor=(0.5, -0.01),
    )
    fig.suptitle(
        "PACE MOANA, 2025-07-01 daily 0.1°: the two traps a naive reader walks into",
        fontsize=11, color=INK_PRIMARY,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))

    if outfile is not None:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        fig.savefig(outfile, dpi=150, facecolor=SURFACE)
        print(f"wrote {outfile}")
    return fig


def plot_clipping(stats, outfile=None):
    """Quantify the clipping, and show the distribution of what survives it.

    Inputs
    ------
    stats : dict
        Output of :func:`summarize_product`.
    outfile : str or None
        Path to write the PNG to.

    Outputs
    -------
    matplotlib.figure.Figure
    """
    fig, (ax_bar, ax_hist) = plt.subplots(1, 2, figsize=(12.5, 4.6), facecolor=SURFACE)

    # --- left: how often each product is exactly zero -------------------------
    keys = [k for k, _v, _l in TAXA]
    ypos = np.arange(len(keys))[::-1]  # first taxon at the top
    for y, key in zip(ypos, keys):
        s = stats[key]
        ax_bar.barh(y, s["pct_zero"], height=0.55, color=TAXON_COLOR[key],
                    zorder=3)
        # Direct labels: required for the aqua slot, and good practice for all.
        ax_bar.text(s["pct_zero"] + 0.4, y, f"{s['pct_zero']:.1f}%",
                    va="center", fontsize=9, color=INK_PRIMARY)
    ax_bar.set_yticks(ypos)
    ax_bar.set_yticklabels([stats[k]["label"] for k in keys], fontsize=9,
                           color=INK_PRIMARY)
    ax_bar.set_xlabel("share of real ocean retrievals equal to exactly 0 [%]",
                      fontsize=9, color=INK_MUTED)
    ax_bar.set_xlim(0, max(s["pct_zero"] for s in stats.values()) * 1.25 + 1)
    ax_bar.xaxis.grid(True, color=GRIDLINE, lw=0.8, zorder=0)
    ax_bar.set_axisbelow(True)
    ax_bar.set_facecolor(SURFACE)
    _style_axes(ax_bar)
    ax_bar.set_title("Silent clipping", fontsize=9, color=INK_PRIMARY)
    # The valid_max exceedances are 1-3 orders of magnitude smaller than the clipping
    # rates, so a second bar series would be invisible; annotate instead. Placed in
    # the empty lower-right of the panel rather than the title, which overflowed.
    note = "\n".join(
        f"{stats[k]['label']}: {stats[k]['pct_over']:.2f}% exceed valid_max"
        for k in keys
    )
    ax_bar.text(0.98, 0.06, note, transform=ax_bar.transAxes, ha="right", va="bottom",
                fontsize=8, color=INK_MUTED, linespacing=1.5)

    # --- right: distribution of the surviving positive values ------------------
    for key, _var, label in TAXA:
        vals = stats[key]["positive"]
        if vals.size == 0:
            continue
        ax_hist.hist(np.log10(vals), bins=60, histtype="step", lw=2.0,
                     color=TAXON_COLOR[key], label=label, zorder=3)
    ax_hist.set_xlabel("log$_{10}$(cells mL$^{-1}$)", fontsize=9, color=INK_MUTED)
    ax_hist.set_ylabel("pixels", fontsize=9, color=INK_MUTED)
    ax_hist.yaxis.grid(True, color=GRIDLINE, lw=0.8, zorder=0)
    ax_hist.set_axisbelow(True)
    ax_hist.set_facecolor(SURFACE)
    ax_hist.legend(frameon=False, fontsize=9, labelcolor=INK_PRIMARY)
    ax_hist.set_title("Abundance of what survives clipping", fontsize=9,
                      color=INK_PRIMARY)
    _style_axes(ax_hist)

    fig.tight_layout()
    if outfile is not None:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        fig.savefig(outfile, dpi=150, facecolor=SURFACE)
        print(f"wrote {outfile}")
    return fig


def plot_mapping_verdict(n_pixels=60_000, seed=0, outfile=None):
    """Our Synechococcus retrieval vs NASA's, under both PC mappings (§7.1).

    The figure that settles the mapping question: two log-log density panels
    against the shipped product — the operational mapping sits on the 1:1
    line to within the L2-vs-L3M compositing scatter; the ATBD mapping is
    displaced by its relocated coefficient. Synechococcus needs no SST, so
    the comparison is ancillary-free.

    Needs both cached granules (AOP Rrs + MOANA) under ``$OS_COLOR/PACE``.
    """
    import xarray as xr
    from epft_up.moana import run_moana

    aop = xr.open_dataset(AOP_GRANULE)
    moana = xr.open_dataset(GRANULE)
    rrs = aop["Rrs"].sel(lat=moana["lat"].values, lon=moana["lon"].values,
                         method="nearest", tolerance=1e-3)
    wave = aop["wavelength"].values.astype(float)
    nasa = moana["syncoccus_moana"].values
    valid = np.isfinite(nasa) & (nasa != LAND) & (nasa != FILL) & (nasa > 0)
    iy, ix = np.nonzero(valid)
    pick = np.random.default_rng(seed).choice(iy.size, min(n_pixels, iy.size),
                                              replace=False)
    iy, ix = iy[pick], ix[pick]
    spectra = rrs.values[iy, ix, :]
    theirs = np.log10(nasa[iy, ix].astype(float))
    aop.close(), moana.close()

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6), facecolor=SURFACE,
                             sharex=True, sharey=True)
    lims = (0, 5)
    for ax, mapping, label in zip(
            axes, ("operational", "atbd"),
            ("operational mapping (picophyt.json: PC16)",
             "ATBD mapping (U13)")):
        out = run_moana(wave, spectra, sst=None, pc_mapping=mapping)
        ok = np.isfinite(out["syn"]) & (out["syn"] > 0)
        mine = np.log10(out["syn"][ok])
        dlog = mine - theirs[ok]
        # Density of points: single-hue sequential ramp on a log count scale.
        ax.hexbin(theirs[ok], mine, gridsize=70, cmap=SEQ_CMAP,
                  bins="log", extent=(*lims, *lims), linewidths=0)
        ax.plot(lims, lims, color=BASELINE, lw=1.2, ls=(0, (4, 3)), zorder=3)
        ax.text(1.62, 1.28, "1:1", color=INK_MUTED, ha="left", fontsize=9,
                rotation=45, rotation_mode="anchor")
        ax.text(0.04, 0.96,
                f"median $\\Delta$log$_{{10}}$ = {np.median(dlog):+.3f}\n"
                f"MAD = {np.median(np.abs(dlog)):.3f}",
                transform=ax.transAxes, va="top", fontsize=10,
                color=INK_PRIMARY)
        ax.set_title(label, fontsize=11, color=INK_PRIMARY)
        ax.set_xlabel("NASA product  log$_{10}$ Syn [cells mL$^{-1}$]",
                      color=INK_PRIMARY)
        ax.set_xlim(lims), ax.set_ylim(lims)
        ax.set_aspect("equal")
        _style_axes(ax)
    axes[0].set_ylabel("our retrieval  log$_{10}$ Syn [cells mL$^{-1}$]",
                       color=INK_PRIMARY)
    fig.suptitle("Which PC mapping does the shipping product use? "
                 "(2025-07-01 daily granule pair)", fontsize=12,
                 color=INK_PRIMARY)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=150, facecolor=SURFACE)
        print(f"wrote {outfile}")
    plt.close(fig)


def _retrieve_both_mappings(n_pixels=None, seed=0):
    """Run our retrieval under both PC mappings on the PACE granule's Rrs.

    Both runs use the *same* constant SST, which is legitimate because the two
    mappings share the intercept and the ``log10(SST)`` term: those cancel in
    the Prochlorococcus difference, so ``pro_op - pro_atbd`` is exact and
    SST-independent even though absolute Pro is not recoverable without the
    real GHRSST field. ``nasa_compat=False`` keeps raw floats, so neither the
    negative-Pro clamp nor the int32 truncation contaminates the comparison.

    Inputs
    ------
    n_pixels : int or None
        Random subsample of valid ocean pixels; None uses all of them.
    seed : int
        Subsample seed.

    Outputs
    -------
    dict with ``out`` ({mapping: retrieval dict}), ``iy``/``ix`` (pixel
    indices into the granule grid), ``shape``, ``lat``, ``lon``.
    """
    import xarray as xr
    from epft_up.moana import run_moana

    aop = xr.open_dataset(AOP_GRANULE)
    moana = xr.open_dataset(GRANULE)
    rrs = aop["Rrs"].sel(lat=moana["lat"].values, lon=moana["lon"].values,
                         method="nearest", tolerance=1e-3)
    wave = aop["wavelength"].values.astype(float)
    nasa = moana["syncoccus_moana"].values
    nasa_pro = moana["prococcus_moana"].values
    lat = moana["lat"].values
    lon = moana["lon"].values

    # Real retrievals only: not land-254, not fill, positive in NASA's product.
    valid = np.isfinite(nasa) & (nasa != LAND) & (nasa != FILL) & (nasa > 0)
    iy, ix = np.nonzero(valid)
    if n_pixels is not None and iy.size > n_pixels:
        pick = np.random.default_rng(seed).choice(iy.size, n_pixels,
                                                  replace=False)
        iy, ix = iy[pick], ix[pick]
    spectra = rrs.values[iy, ix, :]
    shape = nasa.shape
    pro_ref = nasa_pro[iy, ix].astype(float)   # NASA Pro, for scaling
    aop.close(), moana.close()

    out = {m: run_moana(wave, spectra, sst=20.0, pc_mapping=m,
                        nasa_compat=False)
           for m in ("operational", "atbd")}
    return {"out": out, "iy": iy, "ix": ix, "shape": shape,
            "lat": lat, "lon": lon, "pro_ref": pro_ref}


def plot_mapping_consequence(data=None, outfile=None):
    """How much do the misplaced coefficients change the product? (§7.1)

    Distinct from ``plot_mapping_verdict``, which asks *which* mapping NASA
    ships. This asks what the misplacement costs: our own retrieval under the
    operational (in-service) mapping against the same retrieval under the
    ATBD (as-published) mapping, on the one PACE granule.

    Inputs
    ------
    data : dict, optional — from :func:`_retrieve_both_mappings`; computed if
        omitted.
    outfile : str or None — PNG path.

    Outputs
    -------
    matplotlib.figure.Figure
    """
    d = data if data is not None else _retrieve_both_mappings()
    op, at = d["out"]["operational"], d["out"]["atbd"]
    iy, ix, shape = d["iy"], d["ix"], d["shape"]
    lat, lon = d["lat"], d["lon"]
    extent = [lon.min() - 0.05, lon.max() + 0.05,
              lat.min() - 0.05, lat.max() + 0.05]

    # Synechococcus: both mappings are SST-free, so the ratio is exact.
    ok_s = (np.isfinite(op["syn"]) & np.isfinite(at["syn"])
            & (op["syn"] > 0) & (at["syn"] > 0))
    dlog_s = np.log10(op["syn"][ok_s]) - np.log10(at["syn"][ok_s])
    # Prochlorococcus: the difference is exact (intercept + SST cancel).
    ok_p = np.isfinite(op["pro"]) & np.isfinite(at["pro"])
    dpro = op["pro"][ok_p] - at["pro"][ok_p]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), facecolor=SURFACE)

    # (a) Syn, operational vs ATBD, log-log density
    ax = axes[0]
    x = np.log10(at["syn"][ok_s])
    y = np.log10(op["syn"][ok_s])
    # Robust limits floored at 0: our raw (unclipped) retrieval has a thin
    # tail below 1 cell/mL that would otherwise squash the populated range.
    lims = (max(0.0, float(np.nanpercentile(np.r_[x, y], 1.0))),
            float(np.nanpercentile(np.r_[x, y], 99.9)))
    ax.hexbin(x, y, gridsize=70, cmap=SEQ_CMAP, bins="log",
              extent=(*lims, *lims), linewidths=0)
    ax.plot(lims, lims, color=BASELINE, lw=1.2, ls=(0, (4, 3)), zorder=3)
    med = float(np.median(dlog_s))
    ax.text(0.04, 0.96,
            f"median $\\Delta$log$_{{10}}$ = {med:+.3f}\n"
            f"i.e. operational is {10**med:.2f}$\\times$ the ATBD value\n"
            f"n = {ok_s.sum():,}",
            transform=ax.transAxes, va="top", fontsize=9.5, color=INK_PRIMARY)
    ax.set_xlabel("ATBD mapping (as published, U13)  log$_{10}$ Syn",
                  fontsize=9.5, color=INK_PRIMARY)
    ax.set_ylabel("operational mapping (in service, PC16)  log$_{10}$ Syn",
                  fontsize=9.5, color=INK_PRIMARY)
    ax.set_xlim(lims), ax.set_ylim(lims), ax.set_aspect("equal")
    ax.set_title("$\\it{Synechococcus}$: the two mappings disagree\n"
                 "systematically, not randomly", fontsize=10.5,
                 color=INK_PRIMARY)
    _style_axes(ax)

    # (b) map of the Syn ratio
    ax = axes[1]
    field = np.full(shape, np.nan)
    field[iy[ok_s], ix[ok_s]] = dlog_s
    v = float(np.nanpercentile(np.abs(dlog_s), 95))
    im = ax.imshow(field, extent=extent, origin="upper", cmap=DIV_CMAP,
                   vmin=-v, vmax=v, interpolation="nearest")
    ax.set_title("$\\Delta$log$_{10}$ $\\it{Synechococcus}$\n"
                 "(operational $-$ ATBD)", fontsize=10.5, color=INK_PRIMARY)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.06,
                      fraction=0.045)
    cb.set_label("$\\Delta$log$_{10}$ cells mL$^{-1}$", fontsize=8.5,
                 color=INK_MUTED)
    cb.ax.tick_params(labelsize=7, colors=INK_MUTED)
    cb.outline.set_edgecolor(BASELINE)
    _style_axes(ax)

    # (c) map of the Pro difference (SST-free, therefore exact)
    ax = axes[2]
    field = np.full(shape, np.nan)
    field[iy[ok_p], ix[ok_p]] = dpro / 1e3     # 10^3 cells/mL: compact ticks
    v = float(np.nanpercentile(np.abs(dpro / 1e3), 95))
    im = ax.imshow(field, extent=extent, origin="upper", cmap=DIV_CMAP,
                   vmin=-v, vmax=v, interpolation="nearest")
    ax.set_title("$\\Delta$ $\\it{Prochlorococcus}$\n"
                 "(operational $-$ ATBD; SST cancels)", fontsize=10.5,
                 color=INK_PRIMARY)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.06,
                      fraction=0.045)
    cb.set_label("$\\Delta$  $10^3$ cells mL$^{-1}$", fontsize=8.5,
                 color=INK_MUTED)
    cb.ax.tick_params(labelsize=7, colors=INK_MUTED)
    cb.outline.set_edgecolor(BASELINE)
    # Scale to NASA's shipped Pro, which supplies the SST term we cannot
    # reconstruct: this turns an abstract offset into a relative error.
    ref = d["pro_ref"][ok_p]
    frac = np.abs(dpro[ref > 0]) / ref[ref > 0]
    ax.text(0.03, 0.03,
            f"median |$\\Delta$| = {np.median(np.abs(dpro))/1e3:,.0f}"
            f"$\\times10^3$ cells mL$^{{-1}}$\n"
            f"= {100*np.median(frac):.0f}% of NASA's own value",
            transform=ax.transAxes, va="bottom", fontsize=9,
            color=INK_PRIMARY)
    _style_axes(ax)

    fig.suptitle("What the misplaced coefficients cost: our retrieval under "
                 "the operational vs the as-published mapping "
                 "(2025-07-01 granule)", fontsize=12, color=INK_PRIMARY)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    if outfile:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        fig.savefig(outfile, dpi=150, facecolor=SURFACE)
        print(f"wrote {outfile}")
    plt.close(fig)
    return fig


#: Cruise identity for the held-out figure (categorical slots 1-3; the taxon
#: is the facet there, so colour is free to carry the cruise).
CRUISE_COLOR = {23: "#2a78d6", 25: "#eb6834", 28: "#1baf7a"}

#: Lange et al. (2020) Table 3, held-out cruises via Aqua-MODIS -- the only
#: previously published held-out numbers, for reference in the same panels.
LANGE_HELDOUT = {"pro": (1.75, 2.26, 0.54), "syn": (0.93, 2.20, 0.40),
                 "peuk": (1.05, 1.53, 0.60)}


def heldout_skill_data():
    """Run target (ii) and return its matchups and metrics (report §12.3).

    Needs the Brewin et al. (2023) in-situ hyperspectral Rrs and the AMT23/25/28
    BODC flow-cytometry deposits, both under ``$OS_COLOR/AMT/``.

    Inputs
    ------
    None.

    Outputs
    -------
    dict — ``pairs`` : DataFrame (cruise, taxon, pred, obs) of every matchup;
    ``metrics`` : per-cruise and pooled metric dicts.
    """
    from epft_up.moana.io import load_fcm
    from epft_up.moana.validation import validate_heldout_cruises

    tables = {c: load_fcm(cruise=c) for c in (23, 25, 28)}
    res = validate_heldout_cruises(fcm_tables=tables, verbose=False)
    return {"pairs": res["pairs"], "metrics": res["metrics"]}


def plot_heldout_skill(data=None, outfile=None):
    """Retrieved vs observed on Lange's held-out cruises (report §12.3).

    One panel per taxon, points coloured by cruise, with 1:1 and factor-of-3
    guides. This is the figure behind the section's table: picoeukaryotes fall
    on the 1:1 line, while Prochlorococcus sits high and Synechococcus fans
    out -- the transferability failure the table reports as numbers.

    Inputs
    ------
    data : dict, optional — from :func:`heldout_skill_data`; computed if omitted.
    outfile : str or None — PNG path.

    Outputs
    -------
    matplotlib.figure.Figure
    """
    d = data if data is not None else heldout_skill_data()
    pairs, metrics = d["pairs"], d["metrics"]

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.6), facecolor=SURFACE)
    order = [("pro", "Prochlorococcus"), ("syn", "Synechococcus"),
             ("peuk", "picoeukaryotes")]

    for ax, (taxon, label) in zip(axes, order):
        sub = pairs[pairs["taxon"] == taxon]
        sub = sub[np.isfinite(sub["pred"]) & np.isfinite(sub["obs"])
                  & (sub["pred"] > 0) & (sub["obs"] > 0)]
        # Robust limits: a single wild retrieval (Syn has one at ~1e-5)
        # would otherwise compress the populated range to invisibility.
        both = np.r_[sub["pred"].to_numpy(), sub["obs"].to_numpy()]
        lo = float(np.percentile(both, 1.0)) * 0.5
        hi = float(np.percentile(both, 99.0)) * 2.0
        n_out = int(((both < lo) | (both > hi)).sum())

        # Guides first, so the data sit on top of them.
        ax.plot([lo, hi], [lo, hi], color=BASELINE, lw=1.2, ls=(0, (4, 3)),
                zorder=1)
        for f in (3.0, 1 / 3.0):
            ax.plot([lo, hi], [lo * f, hi * f], color=GRIDLINE, lw=1.0,
                    ls=(0, (1, 2)), zorder=1)

        for cruise, colour in CRUISE_COLOR.items():
            c = sub[sub["cruise"] == cruise]
            ax.scatter(c["obs"], c["pred"], s=42, facecolor=colour,
                       edgecolor=SURFACE, linewidth=0.8, alpha=0.9,
                       label=f"AMT{cruise}", zorder=3)

        m = metrics["pooled"][taxon]
        lb, lm, lr = LANGE_HELDOUT[taxon]
        r2_scale = "linear" if taxon == "pro" else "log$_{10}$"
        ax.text(0.04, 0.96,
                f"pooled  n = {m['n']}\n"
                f"bias {m['bias']:.2f}   MAE {m['mae']:.2f}\n"
                f"R$^2$ {m['r2']:+.2f} ({r2_scale})",
                transform=ax.transAxes, va="top", fontsize=9.5,
                color=INK_PRIMARY)
        ax.text(0.96, 0.06,
                f"Lange+2020 held-out\n(MODIS): {lb:.2f} / {lm:.2f} / {lr:+.2f}",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=8,
                color=INK_MUTED)

        if n_out:
            ax.text(0.04, 0.09,
                    f"{n_out} point{'s' if n_out > 1 else ''} beyond axes",
                    transform=ax.transAxes, fontsize=8, color=INK_MUTED)
        ax.set_xscale("log"), ax.set_yscale("log")
        ax.set_xlim(lo, hi), ax.set_ylim(lo, hi), ax.set_aspect("equal")
        ax.set_title(f"$\\it{{{label}}}$", fontsize=11, color=INK_PRIMARY)
        ax.set_xlabel("observed (flow cytometry)  cells mL$^{-1}$",
                      fontsize=9.5, color=INK_PRIMARY)
        ax.grid(True, which="major", color=GRIDLINE, lw=0.7, zorder=0)
        ax.set_axisbelow(True)
        _style_axes(ax)

    axes[0].set_ylabel("retrieved (MOANA, in-situ hyperspectral Rrs)  "
                       "cells mL$^{-1}$", fontsize=9.5, color=INK_PRIMARY)
    handles, labels = axes[0].get_legend_handles_labels()
    # Legend under the title: the bottom-centre slot collides with the
    # middle panel's x-axis label.
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               fontsize=9.5, labelcolor=INK_PRIMARY,
               bbox_to_anchor=(0.5, 0.945))
    fig.suptitle("Held-out skill: the published MOANA coefficients on "
                 "AMT23/25/28 in-situ hyperspectral Rrs  "
                 "(dashed 1:1, dotted $\\pm$3$\\times$)",
                 fontsize=12, color=INK_PRIMARY)
    # Explicit margins, not tight_layout: the equal-aspect log panels fight
    # its rect and the x-axis labels end up clipped off the canvas.
    fig.subplots_adjust(left=0.062, right=0.995, top=0.845,
                        bottom=0.155, wspace=0.22)
    if outfile:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        fig.savefig(outfile, dpi=150, facecolor=SURFACE)
        print(f"wrote {outfile}")
    plt.close(fig)
    return fig


def main():
    """Build every report figure and print a short summary of the tables."""
    wave, loadings, coefs = load_moana_tables()
    used = pcs_used_by_taxon(coefs)

    # Sanity checks worth printing: these are the invariants the report relies on.
    gram = loadings.T @ loadings
    offdiag = np.abs(gram - np.diag(np.diag(gram))).max()
    print(f"loadings shape      : {loadings.shape}")
    print(f"wavelengths         : {wave.min()}-{wave.max()} nm, step {np.diff(wave)[0]}")
    print(f"max |offdiag(V^T V)|: {offdiag:.2e}  (columns orthonormal)")
    print(f"||Rrs'||_2 for any standardised spectrum: {np.sqrt(len(wave) - 1):.3f}")
    for taxon in ("pro", "syn", "apeuk"):
        print(f"  {taxon:6s} uses PCs {sorted(used[taxon])}")

    plot_pc_loadings(
        wave, loadings, used,
        outfile=os.path.join(FIG_DIR, "moana_pc_loadings.png"),
    )

    # --- data-dependent figures: skip cleanly if the granule is not present ----
    if not os.path.isfile(GRANULE):
        print(f"\nskipping product figures: granule not found at {GRANULE}")
        return
    print(f"\nreading {os.path.basename(GRANULE)}")
    fields, lat, lon = load_granule()
    stats = summarize_product(fields)
    for key, _var, label in TAXA:
        s = stats[key]
        print(f"  {label:16s} n_ocean={s['n_ocean']:>8,d}  zero={s['pct_zero']:5.1f}%  "
              f"over valid_max={s['pct_over']:5.2f}%")
    plot_product_masks(fields, lat, lon,
                       outfile=os.path.join(FIG_DIR, "moana_product_masks.png"))
    plot_clipping(stats, outfile=os.path.join(FIG_DIR, "moana_clipping.png"))

    # --- the mapping-verdict figure additionally needs the AOP Rrs granule -----
    if not os.path.isfile(AOP_GRANULE):
        print(f"\nskipping mapping-verdict figure: no AOP granule at {AOP_GRANULE}")
        return
    plot_mapping_verdict(
        outfile=os.path.join(FIG_DIR, "moana_mapping_verdict.png"))

    # Same granule pair, but our retrieval under each mapping against the
    # other -- what the misplacement costs, rather than which NASA ships.
    both = _retrieve_both_mappings()
    plot_mapping_consequence(
        data=both,
        outfile=os.path.join(FIG_DIR, "moana_mapping_consequence.png"))

    # --- held-out skill needs the AMT trees, not the PACE granules -------
    try:
        hd = heldout_skill_data()
    except Exception as exc:            # missing AMT data is not an error
        print(f"\nskipping held-out skill figure: {exc}")
        return
    plot_heldout_skill(
        data=hd, outfile=os.path.join(FIG_DIR, "moana_heldout_skill.png"))


if __name__ == "__main__":
    main()
