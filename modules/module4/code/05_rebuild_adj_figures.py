#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

FIG = Path(r"d:\_bioinformation\modules\module4\figures")
TAB = Path(r"d:\_bioinformation\modules\module4\tables")
plt.rcParams.update({
    "font.family": "Arial", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42,
})

def save(fig, stem, w, h):
    fig.set_size_inches(w, h)
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

mrg = pd.read_csv(TAB / "M4_r2_subtype_adjusted_ZEB1.tsv", sep="\t")
d = mrg[mrg["gene"] != "ZEB1"].copy()
keys = ["LMO2", "LYL1", "ZEB2", "GATA3", "TCF7", "BCL11B", "MEF2C", "CD34", "IL7R"]
lab = d[d["gene"].isin(keys)]
fig, ax = plt.subplots()
ax.axhline(0, color="0.8", lw=0.6)
ax.axvline(0, color="0.8", lw=0.6)
ax.plot([-20, 20], [-20, 20], ls="--", color="0.4", lw=0.8)
ax.scatter(d["t_unadj"], d["t_adj"], s=4, c="0.45", alpha=0.35, linewidths=0)
ax.scatter(lab["t_unadj"], lab["t_adj"], s=22, c="#E64B35", zorder=3)
for _, r in lab.iterrows():
    ax.annotate(r["gene"], (r["t_unadj"], r["t_adj"]), fontsize=8, xytext=(3, 3),
                textcoords="offset points")
ax.set_xlim(-12, 18)
ax.set_ylim(-12, 18)
ax.set_xlabel("t (ZEB1 only)")
ax.set_ylabel("t (ZEB1 + subtype + age)")
ax.set_title("Unadjusted vs subtype-adjusted ZEB1 t-statistic")
ax.text(0.02, 0.98, "Spearman rho = 0.99 (ZEB1 self-correlation excluded)\nLMO2 loses FDR after subtype; LYL1 / T-lineage genes remain",
        transform=ax.transAxes, va="top", fontsize=8, color="0.25")
save(fig, "M4_r2_4.1_unadj_vs_adj_scatter", 6.6, 6.2)

kg = pd.read_csv(TAB / "M4_r2_keygenes_adjusted.tsv", sep="\t")
kg = kg[kg["gene"] != "ZEB1"].sort_values("logFC_adj")
y = np.arange(len(kg))
fig, ax = plt.subplots()
ax.axvline(0, color="0.55", ls="--", lw=0.8)
cols = ["#3C5488" if v < 0 else "#E64B35" for v in kg["logFC_adj"]]
ax.hlines(y, 0, kg["logFC_adj"], color="#B07D3A", lw=1.6)
ax.scatter(kg["logFC_adj"], y, c=cols, s=42, edgecolors="black", linewidths=0.4, zorder=3)
ax.set_yticks(y)
ax.set_yticklabels(kg["gene"])
ax.set_xlabel("logFC per SD ZEB1 (subtype-adjusted)")
ax.set_title("Pre-specified genes after subtype + age")
save(fig, "M4_r2_4.1_keygene_forest", 6.8, 5.8)
print("rebuilt scatter and forest")
