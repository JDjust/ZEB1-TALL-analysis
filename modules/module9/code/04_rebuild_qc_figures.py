#!/usr/bin/env python3
"""Rebuild Module 9 paper figures without the failed GSE110632 polyA KD."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(r"d:\_bioinformation\modules\module9")
FIG = ROOT / "figures"
TAB = ROOT / "tables"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save(fig, stem, w, h):
    fig.set_size_inches(w, h)
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def forest_110635():
    d = pd.read_csv(TAB / "M9_keygene_effects.tsv", sep="\t")
    d = d[d["experiment"] == "gse110635_combined"].copy()
    order = ["ZEB1", "TLX1", "GATA3", "TCF7", "BCL11B", "IL7R", "DNTT",
             "RAG1", "LEF1", "LYL1", "RUNX1", "MYC", "NOTCH1"]
    d = d[d["gene"].isin(order)].copy()
    d["gene"] = pd.Categorical(d["gene"], categories=order[::-1], ordered=True)
    d = d.sort_values("gene")
    d["lo"] = d["log2FoldChange"] - 1.96 * d["lfcSE"]
    d["hi"] = d["log2FoldChange"] + 1.96 * d["lfcSE"]
    y = np.arange(len(d))
    fig, ax = plt.subplots()
    ax.axvline(0, color="0.55", ls="--", lw=0.8)
    ax.errorbar(d["log2FoldChange"], y, xerr=[d["log2FoldChange"] - d["lo"], d["hi"] - d["log2FoldChange"]],
                fmt="o", color="#3C5488", ecolor="#3C5488", capsize=2, ms=5)
    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"])
    ax.set_xlabel("log2 fold change (KD vs NTC)")
    ax.set_title("GSE110635 TLX1 KD, pre-specified genes")
    ax.set_subtitle = None
    fig.text(0.12, 0.01, "GSE110632 polyA omitted: TLX1 did not decrease. Whiskers 95% CI.",
             fontsize=8, color="0.25")
    save(fig, "M9_9.8_zeb1_bidirectional_forest", 6.4, 5.4)


def evidence_heatmap():
    d = pd.read_csv(TAB / "M9_9.24_evidence_matrix.tsv", sep="\t")
    d = d[~d["experiment"].str.contains("110632")].copy()
    rename = {
        "gse110635_combined": "GSE110635\ncombined",
        "gse110635_s1": "GSE110635\nsiTLX1-1",
        "gse110635_s2": "GSE110635\nsiTLX1-2",
    }
    d["exp"] = d["experiment"].map(rename)
    genes = ["ZEB1", "TLX1", "GATA3", "TCF7", "BCL11B", "IL7R", "DNTT",
             "RAG1", "LEF1", "LYL1", "RUNX1", "MYC", "NOTCH1"]
    cols = ["GSE110635\ncombined", "GSE110635\nsiTLX1-1", "GSE110635\nsiTLX1-2"]
    mat = (d.pivot_table(index="gene", columns="exp", values="log2FoldChange")
             .reindex(index=genes, columns=cols))
    fig, ax = plt.subplots()
    vmax = 1.8
    im = ax.imshow(mat.values, cmap="RdBu_r", aspect="auto",
                   norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if abs(v) > 0.9 else "black")
    ax.set_title("TLX1 KD evidence matrix (log2FC)\nGSE110635 only; polyA KD discarded")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03, label="log2FC")
    save(fig, "M9_9.24_evidence_heatmap", 5.8, 6.4)
    out = d.copy()
    out.to_csv(TAB / "M9_9.24_evidence_matrix.tsv", sep="\t", index=False)


if __name__ == "__main__":
    forest_110635()
    evidence_heatmap()
    print("rebuilt 9.8 and 9.24")
