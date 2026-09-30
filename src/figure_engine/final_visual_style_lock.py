"""Single visual style source for the frozen ZEB1/T-ALL figure release.

This module contains display constants only. It does not read or transform
biological data, fit models, select genes, or define inferential thresholds.
"""

import matplotlib
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

INK = '#252525'
ZEB2 = '#3C5488'
ZEB1 = '#E64B35'
BCL11B = '#3C5488'
ETP = '#00A087'
TLX3 = '#E64B35'
GREY = '#A9A9A9'
BG_POINT = '#D9D9D9'
REF = '#BDBDBD'
GRID = '#ECECEC'
NA = '#E6E6E6'
WHITE = '#FFFFFF'
MIDPOINT = '#F7F7F7'

SIGNED_CMAP = LinearSegmentedColormap.from_list(
    'frozen_zeb_signed', [ZEB2, MIDPOINT, ZEB1], N=256)
SIGNED_CMAP.set_bad(NA)
SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list(
    'frozen_zeb_sequential', ['#F3F3F3', ZEB2], N=256)
SEQUENTIAL_CMAP.set_bad(NA)
BINARY_CMAP = ListedColormap([NA, ZEB2], name='frozen_binary')
BINARY_CMAP.set_bad(NA)


def apply_style():
    """Set final-size typography and vector/PDF export defaults."""
    matplotlib.rcParams.update({
        'font.family': 'Arial', 'font.size': 7.2,
        'axes.titlesize': 7.8, 'axes.labelsize': 7.3,
        'xtick.labelsize': 6.7, 'ytick.labelsize': 6.7,
        'legend.fontsize': 6.5, 'axes.linewidth': .55,
        'pdf.fonttype': 42, 'ps.fonttype': 42,
        'figure.facecolor': WHITE, 'savefig.facecolor': WHITE,
        'axes.spines.top': False, 'axes.spines.right': False,
    })


def footer_colorbar(fig, mappable, label, x, width, ticks, y=.018):
    """Small figure-footer key; it does not move or cover a data panel."""
    cax = fig.add_axes([x, y, width, .0055])
    cb = fig.colorbar(mappable, cax=cax, orientation='horizontal', ticks=ticks)
    cb.ax.tick_params(labelsize=6.5, length=1, pad=0)
    cax.set_title(label, fontsize=6.5, pad=1, color=INK)
    return cb


def lock_legacy_module(module):
    """Supply one palette to imported final plotting helpers at runtime."""
    for key, value in {
        'Z1': ZEB1, 'Z2': ZEB2, 'ETP': ETP, 'BCL': BCL11B,
        'TLX': TLX3, 'GREY': GREY, 'LIGHT': REF, 'DARK': INK,
        'SOFT': MIDPOINT, 'CMAP': SIGNED_CMAP,
    }.items():
        if hasattr(module, key):
            setattr(module, key, value)
    if hasattr(module, 'panel'):
        module.panel = panel
    if hasattr(module, 'clean_grid'):
        module.clean_grid = clean_grid
    if hasattr(module, 'color_subtype'):
        module.color_subtype = color_subtype


def color_subtype(value):
    label = str(value).lower()
    if 'bcl11b' in label:
        return BCL11B
    if 'etp' in label:
        return ETP
    if 'tlx3' in label:
        return TLX3
    return GREY


def panel(ax, letter, title, subtitle=None):
    """Common panel lettering for imported final production figures."""
    if letter:
        ax.text(-.02, 1.10, letter, transform=ax.transAxes, fontsize=10,
                weight='bold', color=INK, va='bottom', clip_on=False)
    ax.text(.07 if letter else .0, 1.10, title, transform=ax.transAxes,
            fontsize=7.8, weight='semibold', color=INK,
            va='bottom', clip_on=False)
    if subtitle:
        ax.text(.07 if letter else .0, 1.02, subtitle,
                transform=ax.transAxes, fontsize=6.5,
                color=GREY, va='bottom', clip_on=False)
    ax.tick_params(length=2.2, width=.45, pad=2)
    ax.spines['left'].set_color(REF)
    ax.spines['bottom'].set_color(REF)


def clean_grid(ax, axis='y'):
    ax.grid(axis=axis, color=GRID, linewidth=.38, zorder=0)
    ax.set_axisbelow(True)


def polish_text(fig):
    """Enforce the 6.5 pt minimum after all artists have been created."""
    from matplotlib.text import Text
    for artist in fig.findobj(match=Text):
        if artist.get_text().strip():
            artist.set_fontfamily('Arial')
            if artist.get_fontsize() < 6.5:
                artist.set_fontsize(6.5)
