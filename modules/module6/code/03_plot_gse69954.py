#!/usr/bin/env python3
"""Module 6 r2 figures for GSE69954 ZEB1 450k promoter methylation."""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(r"d:\_bioinformation\modules\module6")
FIG, TAB = ROOT / "figures", ROOT / "tables"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
PAL = {"CIMP-": "#4DBBD5", "CIMP+": "#E64B35"}


def save(fig, stem, w, h):
    fig.set_size_inches(w, h)
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def box_points(ax, df, xcol, ycol, order):
    rng = np.random.default_rng(1)
    for i, g in enumerate(order):
        v = df.loc[df[xcol] == g, ycol].dropna().to_numpy()
        if not len(v):
            continue
        ax.boxplot(v, positions=[i], widths=0.45, showfliers=False,
                   patch_artist=True,
                   boxprops=dict(facecolor=PAL[g], alpha=0.85, edgecolor="black", linewidth=0.7),
                   medianprops=dict(color="black", linewidth=1.1),
                   whiskerprops=dict(color="black", linewidth=0.7),
                   capprops=dict(color="black", linewidth=0.7))
        x = i + rng.uniform(-0.12, 0.12, size=len(v))
        ax.scatter(x, v, s=18, c="0.15", alpha=0.7, zorder=3, linewidths=0)


def main():
    samp = pd.read_csv(TAB / "M6_6.17_ZEB1_sample_region_mean.tsv", sep="\t")
    stats = pd.read_csv(TAB / "M6_6.17_ZEB1_probe_CIMP_stats.tsv", sep="\t")
    summ = pd.read_csv(TAB / "M6_6.17_ZEB1_promoter_beta_summary.tsv", sep="\t").iloc[0]
    long = pd.read_csv(TAB / "M6_6.17_ZEB1_promoter_beta.tsv", sep="\t")
    order = ["CIMP-", "CIMP+"]

    prom = samp[samp.region == "promoter"].copy()
    fig, ax = plt.subplots()
    box_points(ax, prom, "cimp", "mean_beta", order)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(order)
    ax.set_ylabel("Mean ZEB1 promoter beta")
    ax.set_title("6.17  ZEB1 promoter methylation (GSE69954)")
    ax.set_ylim(0, max(0.35, float(prom.mean_beta.max()) + 0.05))
    ax.text(0.02, 0.98,
            f"n={int(summ.n_CIMP_neg)} vs {int(summ.n_CIMP_pos)}; "
            f"delta={float(summ.promoter_delta):.3f}; "
            f"g={float(summ.promoter_hedges_g):.2f}; "
            f"P={float(summ.promoter_p):.2g}\n"
            f"35 TSS/UTR probes; promoter remains hypomethylated in both CIMP groups.",
            transform=ax.transAxes, va="top", ha="left", fontsize=8, color="0.25")
    save(fig, "M6_6.17_ZEB1_promoter_beta", 5.6, 5.4)

    body = samp[samp.region == "gene_body"].copy()
    fig, ax = plt.subplots()
    box_points(ax, body, "cimp", "mean_beta", order)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(order)
    ax.set_ylabel("Mean ZEB1 gene-body beta")
    ax.set_title("6.17  ZEB1 gene-body methylation")
    save(fig, "M6_6.17_ZEB1_genebody_beta", 5.4, 5.2)

    st = stats[stats.region == "promoter"].dropna(subset=["mean_all"]).copy()
    st = st.sort_values("pos")
    fig, ax = plt.subplots()
    y = np.arange(len(st))
    ax.scatter(st.mean_CIMP_neg, y, s=22, c=PAL["CIMP-"], label="CIMP-", zorder=3)
    ax.scatter(st.mean_CIMP_pos, y, s=22, c=PAL["CIMP+"], label="CIMP+", zorder=3)
    for _, r in st.iterrows():
        yi = list(st.probe).index(r.probe)
        ax.plot([r.mean_CIMP_neg, r.mean_CIMP_pos], [yi, yi], color="0.7", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(st.probe, fontsize=6)
    ax.set_xlabel("Mean beta")
    ax.axvline(0.2, color="0.6", ls="--", lw=0.7)
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("6.17  Per-CpG promoter beta, CIMP- vs CIMP+")
    save(fig, "M6_6.17_ZEB1_promoter_probe_means", 6.4, 7.2)

    # heatmap of promoter probes (samples ordered by CIMP then mean)
    d = long[(long.region == "promoter") & long.beta.notna()].copy()
    meta = d[["geo", "cimp"]].drop_duplicates()
    mat = d.pivot_table(index="probe", columns="geo", values="beta")
    pos = stats.set_index("probe")["pos"]
    mat = mat.loc[mat.index.intersection(pos.index)]
    mat = mat.loc[pos.loc[mat.index].sort_values().index]
    col_order = (meta.set_index("geo").loc[mat.columns]
                 .assign(_k=lambda x: x.cimp.map({"CIMP-": 0, "CIMP+": 1}))
                 .sort_values(["_k"]).index)
    mat = mat.loc[:, col_order]
    fig, ax = plt.subplots()
    im = ax.imshow(mat.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.4,
                   interpolation="nearest")
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6)
    n_neg = int((meta.set_index("geo").loc[col_order, "cimp"] == "CIMP-").sum())
    ax.axvline(n_neg - 0.5, color="white", lw=1.2)
    ax.set_xticks([])
    ax.set_xlabel(f"CIMP- n={n_neg}  |  CIMP+ n={mat.shape[1]-n_neg}")
    ax.set_title("6.17  ZEB1 promoter CpG beta (0-0.4 scale)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="beta")
    save(fig, "M6_6.17_ZEB1_promoter_heatmap", 8.2, 6.8)

    print("figures written", list(FIG.glob("M6_6.17_*")))


if __name__ == "__main__":
    main()
