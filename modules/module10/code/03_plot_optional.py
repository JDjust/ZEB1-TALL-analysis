#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(r"d:\_bioinformation\modules\module10")
FIG, TAB = ROOT / "figures", ROOT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
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


def main():
    st = pd.read_csv(TAB / "M10_r2_GSE287751_DN1_KO_vs_ctrl.tsv", sep="\t")
    order = ["Lmo2", "Zeb1", "Zeb2", "Tal1", "Lyl1", "Spi1", "Mef2c", "Cd34",
             "Tcf7", "Gata3", "Bcl11b", "Il7r", "Notch1"]
    d = st[st.gene.isin(order)].copy()
    d["gene"] = pd.Categorical(d["gene"], categories=order, ordered=True)
    d = d.sort_values("gene")
    fig, ax = plt.subplots()
    colors = ["#E64B35" if v < 0 else "#3C5488" for v in d.log2FC_KO_vs_ctrl]
    ax.barh(d.gene.astype(str), d.log2FC_KO_vs_ctrl, color=colors, edgecolor="black", linewidth=0.4)
    ax.axvline(0, color="0.3", lw=0.7)
    ax.set_xlabel("log2 FC (Lmo2 KO / control)")
    ax.set_title("GSE287751 DN1, Lmo2 KO vs control")
    ax.text(0.02, 0.02, "n=2 vs 2. Lmo2 dropped; Zeb1 did not.", transform=ax.transAxes, fontsize=8, color="0.3")
    save(fig, "M10_r2_GSE287751_DN1_lmo2KO_log2FC", 6.2, 5.6)

    long = pd.read_csv(TAB / "M10_r2_GSE287751_key_long.tsv", sep="\t")
    sub = long[long.symbol.isin(["Lmo2", "Zeb1"]) & long.stage.eq("DN1")].copy()
    fig, axes = plt.subplots(1, 2)
    pal = {"control": "#4DBBD5", "Lmo2_KO": "#E64B35"}
    for ax, gene in zip(axes, ["Lmo2", "Zeb1"]):
        g = sub[sub.symbol == gene]
        for i, geno in enumerate(["control", "Lmo2_KO"]):
            v = g.loc[g.geno == geno, "TPM"].to_numpy()
            ax.scatter(np.full(len(v), i) + np.array([-0.05, 0.05])[:len(v)], v,
                       s=40, c=pal[geno], zorder=3, edgecolors="black", linewidths=0.4)
            ax.hlines(v.mean(), i - 0.2, i + 0.2, color="black", lw=1.2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["control", "Lmo2 KO"])
        ax.set_title(gene)
        ax.set_ylabel("TPM")
    fig.suptitle("GSE287751 DN1 day-5 OP9-DLL1, n=2", fontsize=11)
    save(fig, "M10_r2_GSE287751_DN1_Lmo2_Zeb1_TPM", 7.2, 4.6)

    lib = pd.read_csv(TAB / "M10_r2_GSE260948_library_means.tsv", sep="\t")
    fig, ax = plt.subplots()
    x = np.arange(len(lib))
    w = 0.25
    ax.bar(x - w, lib.g_Zeb1, w, color="#3C5488", edgecolor="black", lw=0.4, label="Zeb1")
    ax.bar(x, lib.g_Tal1, w, color="#E64B35", edgecolor="black", lw=0.4, label="Tal1")
    ax.bar(x + w, lib.g_Lmo2, w, color="#4DBBD5", edgecolor="black", lw=0.4, label="Lmo2")
    ax.set_xticks(x)
    ax.set_xticklabels(lib.state, rotation=15)
    ax.set_ylabel("Mean log1p (library)")
    ax.legend(frameon=False)
    ax.set_title("GSE260948 Tal1/Lmo2 thymus scRNA")
    ax.text(0.02, 0.98, "n=1 WT, 1 preleukemic, 2 T-ALL libraries. Descriptive only.",
            transform=ax.transAxes, va="top", fontsize=8, color="0.3")
    save(fig, "M10_r2_GSE260948_library_means", 6.8, 5.0)
    print("ok")


if __name__ == "__main__":
    main()
