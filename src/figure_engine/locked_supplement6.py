"""Final visual-only redraw of the frozen seven-figure manuscript set.

All patient, cell and gene values come from locked project outputs. This script
does not fit biological models, select genes, test new hypotheses or change
analysis thresholds. The existing atlas UMAP and measured counts are displayed.
"""
from pathlib import Path
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

import final_visual_style_lock as style
import editorial_rebuild_figures as base
style.apply_style()
style.lock_legacy_module(base)
import editorial_upgrade_347 as prior
import editorial_upgrade_5 as prior5
style.lock_legacy_module(prior)
style.lock_legacy_module(prior5)

ROOT = Path(os.environ.get('ZEB1_DATA_ROOT',str(Path(__file__).resolve().parents[2]/'data/processed'))) / 'total/rebuild_2026'
DATA = ROOT / "data"
DEEP = DATA / "deepen_2026"
SOURCE = DATA / "source_data_rebuilt"
VISUAL = DATA / "editorial_visual_sources"
OUT = Path(__file__).resolve().parents[2] / 'results/figures'
SUPP_OUT = OUT.parent / "supplementary"
BLUE, ORANGE, TEAL = style.ZEB2, style.ZEB1, style.ETP
DARK, GREY, LIGHT = style.INK, style.GREY, style.REF
CMAP = style.SIGNED_CMAP


def read(path):
    return pd.read_csv(path, sep="\t", low_memory=False)

def heading(ax, letter, title):
    ax.text(0, 1.035, letter, transform=ax.transAxes, fontsize=10,
            weight="bold", color=DARK, va="bottom", clip_on=False)
    ax.text(.12, 1.035, title, transform=ax.transAxes, fontsize=7.8,
            weight="semibold", color=DARK, va="bottom", clip_on=False)
    ax.tick_params(length=2, width=.45, pad=2)
    ax.spines[["top", "right"]].set_visible(False)

def clean(ax, axis="y"):
    ax.grid(axis=axis, color=style.GRID, lw=.40, zorder=0)
    ax.set_axisbelow(True)

def sanitize(fig):
    replacements = {"尾R": "βR", "蟻": "ρ", "脳": "×", "路": "·",
                    "鈭扗P": "→DP", "鈭扵": "−T", "鈭扙TP": "−ETP",
                    "鈭?": "−", "鈮?8": "≥18", "鈮?1": "≥21",
                    "P枚l枚nen": "Pölönen", "鈥?": "—"}
    for artist in fig.findobj(match=matplotlib.text.Text):
        val = artist.get_text()
        for old, new in replacements.items():
            val = val.replace(old, new)
        artist.set_text(val)
        if val.strip() and artist.get_fontsize() < 6.5:
            artist.set_fontsize(6.5)
    assert all("鈭" not in t.get_text() and "蟻" not in t.get_text()
               for t in fig.findobj(match=matplotlib.text.Text))

def figS6():
    # The cohort inventory and full 94-gene malignant audit live in Supplement.
    fig=plt.figure(figsize=(180/25.4,175/25.4))
    g=GridSpec(3,2,figure=fig,height_ratios=[.53,1,1],width_ratios=[1,1])
    a=fig.add_subplot(g[0,:]);b=fig.add_subplot(g[1:,0])
    c=fig.add_subplot(g[1,1]);d=fig.add_subplot(g[2,1])
    prior5.panel=lambda ax,letter,title,subtitle=None: heading(ax,"B" if letter=="E" else letter,title)
    prior5._validation_matrix(a)
    prior5._malignant_gene_heatmap(b)
    frozen=read(VISUAL/"F7G_frozen100_context_effects.tsv")
    assert "Discovery βR" in frozen and "GSE146901 non−ETP" in frozen and "Malignant ρ(E,B)" in frozen
    for ax, effect, letter, title in [(c,"GSE146901 non−ETP","C","GSE146901 bulk effect"),
                                     (d,"Malignant ρ(E,B)","D","Malignant-cell effect")]:
        sub=frozen[["Discovery βR",effect,"side"]].dropna()
        if letter == "C":
            assert len(sub) == 84
        else:
            assert len(sub) == 94
        for side,color in [("ZEB1-side50",ORANGE),("ZEB2-side50",BLUE)]:
            q=sub[sub.side.eq(side)]
            ax.scatter(q.iloc[:,0],q.iloc[:,1],s=11,color=color,alpha=.7,
                       edgecolor="white",lw=.2,rasterized=True)
        ax.axvline(0,color=LIGHT,lw=.6);ax.axhline(0,color=LIGHT,lw=.6)
        ax.set_xlabel("Discovery βR");ax.set_ylabel("External effect")
        if letter == "C":
            ax.set_yscale("symlog", linthresh=1)
            ax.set_yticks([-100,-10,0,10,100], ["−100","−10","0","10","100"])
        clean(ax);heading(ax,letter,title)
    fig.subplots_adjust(left=.12,right=.975,top=.95,bottom=.075,wspace=.44,hspace=.39)
    gene_key=matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(-1,1),cmap=CMAP)
    style.footer_colorbar(fig,gene_key,'B: columnwise display-scaled effect',
                          .43,.24,[-1,0,1])
    sanitize(fig)
    style.polish_text(fig)
    fig.savefig(SUPP_OUT/"FigureS6.pdf")
    fig.savefig(SUPP_OUT/"FigureS6.png",dpi=300)
    plt.close(fig)
    print("Figure S6: PDF + PNG",flush=True)
