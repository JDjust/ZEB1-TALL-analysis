"""
Shared publication figure style for the ZEB1/T-ALL reappraisal manuscript.

Goal: Nature/Cell-like main-figure aesthetics that stay fully editable in
Adobe Illustrator (vector PDF, real text objects, TrueType 42 fonts).

All figure scripts import from here so the whole set is visually consistent.
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
TOTAL_DIR = os.path.dirname(CODE_DIR)
ROOT = os.path.dirname(TOTAL_DIR)
MODULES = os.path.join(ROOT, "modules")
FIG_DIR = os.path.join(TOTAL_DIR, "figures")
PNG_DIR = os.path.join(FIG_DIR, "png")
SUPP_DIR = os.path.join(FIG_DIR, "Supplementary")
SUPP_PNG_DIR = os.path.join(SUPP_DIR, "png")
DATA_DIR = os.path.join(TOTAL_DIR, "data")
MS_DIR = os.path.join(TOTAL_DIR, "manuscript")
for d in (FIG_DIR, PNG_DIR, SUPP_DIR, SUPP_PNG_DIR, DATA_DIR, MS_DIR):
    os.makedirs(d, exist_ok=True)


def tbl(module, name):
    """Use live module export, or its bundled figure snapshot when absent."""
    source = os.path.join(MODULES, module, "tables", name)
    if os.path.isfile(source):
        return source
    bundled = os.path.join(TOTAL_DIR, "data", "module_tables", module, name)
    if os.path.isfile(bundled):
        return bundled
    raise FileNotFoundError(f"Missing module table: {source}; bundled fallback: {bundled}")


# ----------------------------------------------------------------------------
# Colour system (colour-blind safe, restrained)
# ----------------------------------------------------------------------------
C = {
    "tall":   "#C1272D",   # T-ALL / focus red
    "ball":   "#0072B2",   # B-ALL / comparator blue
    "aml":    "#E69F00",
    "normal": "#7F7F7F",
    "grey":   "#9AA0A6",
    "grid":   "#E6E8EB",
    "ink":    "#1A1A1A",
    "up":     "#C1272D",
    "down":   "#0072B2",
    "ns":     "#C7CBD1",
    "green":  "#2CA02C",
    "purple": "#6A3D9A",
    "teal":   "#1B9E77",
    "gold":   "#D4A017",
}
# categorical palette for subtypes / lineages
PAL = ["#C1272D", "#0072B2", "#E69F00", "#009E73", "#6A3D9A",
       "#CC79A7", "#56B4E9", "#D55E00", "#999999", "#117733"]

# lineage palette (Fig 1 immune atlas etc.)
LIN = {"T": "#C1272D", "B": "#0072B2", "NK": "#009E73",
       "Myeloid": "#E69F00", "Progenitor": "#6A3D9A", "Other": "#B0B0B0"}


def set_style():
    """Global rcParams. Call once at the top of every figure script."""
    for fam in ("Arial", "Helvetica", "DejaVu Sans"):
        if any(fam in f.name for f in font_manager.fontManager.ttflist):
            base = fam
            break
    else:
        base = "DejaVu Sans"
    mpl.rcParams.update({
        "pdf.fonttype": 42,          # editable text in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": base,
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7.5,
        "axes.titlepad": 4,
        "axes.labelpad": 2.5,
        "axes.linewidth": 0.7,
        "axes.edgecolor": "#3A3A3A",
        "axes.axisbelow": True,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "xtick.color": "#3A3A3A",
        "ytick.color": "#3A3A3A",
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "legend.handlelength": 1.1,
        "legend.handletextpad": 0.4,
        "legend.columnspacing": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        "lines.solid_capstyle": "round",
        "mathtext.default": "regular",
    })


MM = 1 / 25.4  # millimetres -> inches


def new_fig(w_mm, h_mm, constrained=True):
    """Figure sized in millimetres. Constrained layout is default; pass
    constrained=False for nested GridSpecs with long y-labels."""
    fig = plt.figure(figsize=(w_mm * MM, h_mm * MM),
                     layout="constrained" if constrained else None)
    if constrained:
        fig.get_layout_engine().set(w_pad=2.2 / 72, h_pad=2.2 / 72,
                                    wspace=0.06, hspace=0.08)
    return fig


def add_grid(ax, axis="y"):
    """Subtle background grid (top-journal look)."""
    ax.grid(axis=axis, color=C["grid"], lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def repel(ax, texts, arrow=True, **kw):
    """Non-overlapping text labels via adjustText (fallback: leave as-is)."""
    try:
        from adjustText import adjust_text
        ap = dict(arrowprops=dict(arrowstyle="-", color="#9AA0A6", lw=0.4)) if arrow else {}
        adjust_text(texts, ax=ax, expand=(1.15, 1.35),
                    min_arrow_len=2, **ap, **kw)
    except Exception:
        pass


def conf_ellipse(ax, x, y, n_std=1.96, fc="none", ec="#333", lw=0.9, alpha=1.0, ls="-"):
    """Draw an n_std confidence ellipse for a 2D point cloud."""
    import numpy as np
    from matplotlib.patches import Ellipse
    from matplotlib.transforms import Affine2D
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    pear = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    rx, ry = np.sqrt(1 + pear), np.sqrt(1 - pear)
    ell = Ellipse((0, 0), width=2 * rx, height=2 * ry, facecolor=fc,
                  edgecolor=ec, lw=lw, alpha=alpha, ls=ls, zorder=1)
    sx, sy = np.sqrt(cov[0, 0]) * n_std, np.sqrt(cov[1, 1]) * n_std
    tr = Affine2D().rotate_deg(45).scale(sx, sy).translate(x.mean(), y.mean())
    ell.set_transform(tr + ax.transData)
    ax.add_patch(ell)


def panel_label(ax, letter, x=-0.16, y=1.06, size=11):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right", family="sans-serif")


def plabel(ax, letter, dx=-0.5, dy=3):
    """Panel letter placed in figure-relative points (constrained-layout safe)."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx * 12, dy), textcoords="offset points",
                fontsize=11, fontweight="bold", va="bottom", ha="right")


def text_right(ax, x, y, s, pad=0.06, **kw):
    """Place a numeric label to the RIGHT of a point (avoids title collisions)."""
    kw.setdefault("va", "center")
    kw.setdefault("ha", "left")
    kw.setdefault("fontsize", 5.6)
    ax.text(x + pad, y, s, **kw)


def text_above_safe(ax, x, y, s, ypad=0.18, **kw):
    """Small in-panel annotation that stays inside data coords."""
    kw.setdefault("ha", "center")
    kw.setdefault("va", "bottom")
    kw.setdefault("fontsize", 5.5)
    ax.text(x, y + ypad, s, **kw)


def stars(p):
    if p is None:
        return "ns"
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def style_ax(ax):
    ax.tick_params(length=2.5, width=0.7)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#333333")
    return ax


def save(fig, stem, pdf_dir=None, png_dir=None):
    """Save vector PDF (main) + high-res PNG (Word embedding / preview)."""
    pdf_dir = FIG_DIR if pdf_dir is None else pdf_dir
    png_dir = PNG_DIR if png_dir is None else png_dir
    os.makedirs(pdf_dir, exist_ok=True)
    os.makedirs(png_dir, exist_ok=True)
    pdf = os.path.join(pdf_dir, stem + ".pdf")
    png = os.path.join(png_dir, stem + ".png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=400)
    plt.close(fig)
    print("saved", os.path.relpath(pdf, TOTAL_DIR), "+ png")


def save_supp(fig, stem):
    """Save a supplementary figure into figures/Supplementary/."""
    save(fig, stem, pdf_dir=SUPP_DIR, png_dir=SUPP_PNG_DIR)
