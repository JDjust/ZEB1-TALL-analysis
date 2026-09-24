"""
Figure 4 | Selected T-lineage associations, genome-wide breadth, cohort
heterogeneity, and pathway/signature sensitivity.

Panels
  A  Selected gene correlations with ZEB1 (TARGET, n=265)
  B  Six key-gene correlations in three independent bulk cohorts
  C  Genome-wide unadjusted versus adjusted association breadth
  D  Subtype-adjusted volcano with key T-lineage / immature genes labelled
  E  Functional signatures do not reproduce across cohorts (sign flips)
  F  Hallmark GSEA is null (no pathway passes FDR)
  G  LMO2–ZEB1 conditional associations in four cohorts
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import _style as S

S.set_style()
C = S.C

key = pd.read_csv(S.tbl("module4", "M4_keygene_stats.tsv"), sep="\t")
adj = pd.read_csv(S.tbl("module4", "M4_r2_keygenes_adjusted.tsv"), sep="\t")
sigT = pd.read_csv(S.tbl("module4", "M4_4.9_4.19_signature_vs_ZEB1.tsv"), sep="\t")
sigP = pd.read_csv(S.tbl("module4", "M4_4.32_pharmacotype_signature_vs_ZEB1.tsv"), sep="\t")
gw = pd.read_csv(S.tbl("module4", "M4_r2_subtype_adjusted_ZEB1.tsv"), sep="\t")
hall = pd.read_csv(S.tbl("module4", "M4_4.6_hallmark_fgsea.tsv"), sep="\t")
conc = pd.read_csv(S.tbl("module4", "M4_r2_TARGET_vs_Pharmacotype_rank.tsv"), sep="\t")

fig = S.new_fig(183, 155, constrained=False)
axd = fig.subplot_mosaic(
    [["A", "B", "C"],
     ["D", "E", "F"],
     ["G", "G", "G"]],
    gridspec_kw=dict(height_ratios=[1.05, 1.05, 0.62], left=.11, right=.96,
                     top=.94, bottom=.08, wspace=.78, hspace=.68),
)

# ---- A: key gene rho -------------------------------------------------------
axA = axd["A"]
k = key[key.gene != "ZEB1"].copy().sort_values("rho_zeb1")
y = np.arange(len(k))
cols = [C["up"] if (r > 0 and f < 0.05) else (C["down"] if (r < 0 and f < 0.05) else C["ns"])
        for r, f in zip(k["rho_zeb1"], k["fdr_rho"])]
S.add_grid(axA, "x")
axA.barh(y, k["rho_zeb1"], color=cols, height=0.72, edgecolor="white", linewidth=0.3, zorder=3)
axA.axvline(0, color="#333", lw=0.6, zorder=2)
axA.set_yticks(y); axA.set_yticklabels(k["gene"], fontsize=6)
axA.set_xlabel("Spearman ρ with ZEB1")
axA.set_title("Selected genes (TARGET)", loc="left", fontsize=7)
axA.set_xlim(-0.33, 0.5)
S.style_ax(axA); S.plabel(axA, "A")

# ---- B: key-gene replication and sensitivity ------------------------------
axB = axd["B"]
b = pd.read_csv(os.path.join(S.DATA_DIR, "SourceData_Fig4B_crosscohort.tsv"), sep="\t")
genes = ["GATA3", "TCF7", "BCL11B", "IL7R", "LMO2", "ZEB2"]
cohorts = ["TARGET-ALL-P2", "StJude_Pharmacotype", "GSE272023"]
tmat = b.pivot(index="gene", columns="cohort", values="rho").loc[genes, cohorts]
from matplotlib.colors import TwoSlopeNorm
hm = axB.imshow(tmat.to_numpy(), cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), aspect="auto")
axB.set_xticks(range(3), ["TARGET", "St Jude", "NOPHO"], fontsize=5.5)
axB.set_yticks(range(len(genes)), genes, fontsize=6)
axB.tick_params(axis="both", length=0)
axB.set_xticks(np.arange(-.5, 3, 1), minor=True)
axB.set_yticks(np.arange(-.5, len(genes), 1), minor=True)
axB.grid(which="minor", color="white", linewidth=1.2)
axB.tick_params(which="minor", bottom=False, left=False)
for yi, gene in enumerate(genes):
    for xi, cohort in enumerate(cohorts):
        val = float(tmat.loc[gene, cohort])
        axB.text(xi, yi, f"{val:+.2f}", ha="center", va="center", fontsize=5.2,
                 color="white" if abs(val) > .7 else C["ink"])
axB.set_title("Three patient cohorts", loc="left", fontsize=7.5)
cbB = fig.colorbar(hm, ax=axB, fraction=0.048, pad=0.02)
cbB.set_label("Spearman ρ", fontsize=5.6)
cbB.ax.tick_params(labelsize=5, length=2)
S.plabel(axB, "B")

# ---- C: genome-wide hexbin -------------------------------------------------
axC = axd["C"]
hb = axC.hexbin(gw["t_unadj"], gw["t_adj"], gridsize=48, cmap="rocket_r" if False else "magma_r",
                mincnt=1, linewidths=0)
lim = 26
axC.plot([-lim, lim], [-lim, lim], color=C["tall"], lw=0.9, ls="--", zorder=3)
axC.set_xlim(-lim, lim); axC.set_ylim(-lim, lim)
axC.set_xlabel("t (unadjusted)")
axC.set_ylabel("t (subtype + age adj.)")
axC.set_title("Global association breadth", loc="left", fontsize=7)
n_sig = int(gw["fdr_adj"].lt(.05).sum())
n_pos = int(gw["logFC_adj"].gt(0).sum())
axC.text(0.05, 0.96, f"{n_sig:,}/{len(gw):,} FDR<0.05\n"
         f"{n_pos:,}/{len(gw):,} coefficients positive",
         transform=axC.transAxes, fontsize=5.5, va="top",
         bbox=dict(facecolor="white", edgecolor="none", alpha=.88, pad=1.5))
cb = fig.colorbar(hb, ax=axC, fraction=0.046, pad=0.02)
cb.set_label("genes", fontsize=6); cb.ax.tick_params(labelsize=5.5)
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: volcano of subtype-adjusted associations ---------------------------
axD = axd["D"]
v = gw.copy()
v["mlog"] = -np.log10(v["fdr_adj"].clip(lower=1e-300))
S.add_grid(axD, "both")
axD.scatter(v["logFC_adj"], v["mlog"], s=3, c=C["ns"], alpha=0.35, linewidths=0, zorder=2)
hi = v[v["gene"].isin(["GATA3", "TCF7", "BCL11B", "LMO2"])]
cols_h = [C["up"] if fc > 0 else C["down"] for fc in hi["logFC_adj"]]
axD.scatter(hi["logFC_adj"], hi["mlog"], s=22, c=cols_h, edgecolor="#222", lw=0.4, zorder=4)
offsets = {"TCF7": (7, 4), "BCL11B": (-27, -10), "GATA3": (8, -5),
           "LMO2": (-27, 6)}
for _, row in hi.iterrows():
    dx, dy = offsets[row["gene"]]
    axD.annotate(row["gene"], (row["logFC_adj"], min(row["mlog"], 41.5)),
                 xytext=(dx, dy), textcoords="offset points", fontsize=5.3,
                 arrowprops=dict(arrowstyle="-", lw=.4, color="#777"))
axD.axhline(-np.log10(0.05), color="#999", lw=0.6, ls=":", zorder=1)
axD.axvline(0, color="#333", lw=0.5, zorder=1)
axD.set_ylim(0, 42)
axD.set_xlabel("logFC per ZEB1 SD\n(subtype + age adjusted)")
axD.set_ylabel(r"$-$log10 FDR")
axD.set_title("Subtype-adjusted effects", loc="left", fontsize=7)
S.style_ax(axD); S.plabel(axD, "D")

# ---- E: signature scatter across cohorts -----------------------------------
axE = axd["E"]
m = sigT.merge(sigP, on="signature", suffixes=("_T", "_P"))
short = {"Tcell_differentiation": "T-diff", "DNA_repair": "DNA repair",
         "PI3K_AKT_mTOR": "PI3K", "IL7_JAK_STAT": "JAK-STAT",
         "Apoptosis": "Apoptosis", "TGFb": "TGF-β", "NOTCH": "NOTCH",
         "MYC": "MYC", "Stemness": "Stemness", "Proliferation": "Prolif.",
         "ETP": "ETP"}
programs = m.set_index("signature")[["rho_T", "rho_P"]]
imE = axE.imshow(programs, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axE.set_yticks(range(len(programs)),[short.get(x,x) for x in programs.index],fontsize=5.2)
axE.set_xticks([0,1],["TARGET","St Jude"],fontsize=5.5)
for yi,row in enumerate(programs.to_numpy()):
    for xi,value in enumerate(row):
        axE.text(xi,yi,f"{value:+.2f}",ha="center",va="center",fontsize=5.0)
axE.tick_params(length=0)
axE.set_title("Cross-cohort programs",loc="left",fontsize=7)
cbE=fig.colorbar(imE,ax=axE,fraction=.04,pad=.04)
cbE.set_label("Spearman ρ",fontsize=5.5); cbE.ax.tick_params(labelsize=5)
S.plabel(axE,"E")

# ---- F: Hallmark GSEA null -------------------------------------------------
axF = axd["F"]
h = hall.reindex(hall["NES"].abs().sort_values(ascending=False).index).head(8).sort_values("NES")

def _short(p):
    t = p.replace("HALLMARK_", "").replace("_", " ").title()
    t = t.replace("Reactive Oxygen Species Pathway", "ROS")
    t = t.replace("Tnfa Signaling Via Nfkb", "TNFA/NFKB")
    t = t.replace("Pi3K Akt Mtor Signaling", "PI3K/mTOR")
    t = t.replace("Pancreas Beta Cells", "Pancreas β")
    return (t[:22] + "…") if len(t) > 23 else t

y = np.arange(len(h))
S.add_grid(axF, "x")
axF.hlines(y, 0, h["NES"], color="#C7CBD1", lw=1.1, zorder=2)
axF.scatter(h["NES"], y, s=24, color=C["ns"], edgecolor="#333", linewidth=0.4, zorder=3)
axF.axvline(0, color="#333", lw=0.6, zorder=2)
axF.set_yticks(y); axF.set_yticklabels([_short(p) for p in h["pathway"]], fontsize=5.6)
axF.set_xlabel("NES")
axF.set_title("Hallmark enrichment", loc="left", fontsize=7)
axF.text(0.96, 0.08, "0 / 50 FDR<0.05", transform=axF.transAxes,
         fontsize=6, ha="right", va="bottom")
axF.set_xlim(-2.1, 2.1)
S.style_ax(axF); S.plabel(axF, "F")

# ---- G: LMO2-ZEB1 false edge ----------------------------------------------
axG = axd["G"]
edge = pd.read_csv(os.path.join(S.DATA_DIR, "SourceData_FigS4_LMO2_ZEB1_false_edge.tsv"), sep="\t")
labmap = {"TARGET-ALL-P2": "TARGET", "GSE62156": "GSE62156", "GSE26713": "GSE26713",
          "GSE110636": "GSE110636"}
order = ["TARGET-ALL-P2", "GSE62156", "GSE26713", "GSE110636"]
e = edge.set_index("cohort").reindex([c for c in order if c in set(edge.cohort)]).reset_index()
y = np.arange(len(e))
S.add_grid(axG, "x")
axG.axvline(0, color="#333", lw=0.6, zorder=2)
axG.hlines(y, e["rho_raw"], e["rho_adj"], color="#C4C8CE", lw=1.4, zorder=2)
axG.scatter(e["rho_raw"], y, s=28, color=C["grey"], zorder=4, label="raw Spearman", edgecolor="white", lw=0.3)
adjcol = [C["tall"] if p >= 0.05 else C["ink"] for p in e["p_adj"]]
axG.scatter(e["rho_adj"], y, s=36, color=adjcol, zorder=4, label="subtype residual", edgecolor="white", lw=0.3)
axG.set_yticks(y)
axG.set_yticklabels([f"{labmap.get(c, c)} (n={int(n)})" for c, n in zip(e["cohort"], e["n"])], fontsize=6.2)
axG.set_xlabel("Spearman ρ (LMO2 vs ZEB1)")
axG.set_xlim(-0.95, 0.55)
axG.set_title("LMO2–ZEB1 association varies after subtype adjustment",
              loc="left", fontsize=7.5)
axG.legend(loc="lower right", fontsize=5.6)
for yi in y:
    r = e.iloc[yi]
    S.text_right(axG, max(float(r["rho_adj"]), float(r["rho_raw"])), yi,
                 f"P={r['p_adj']:.2g}", pad=0.03, fontsize=5.4,
                 color=C["tall"] if r["p_adj"] >= 0.05 else C["ink"])
S.style_ax(axG); S.plabel(axG, "G")

S.save(fig, "Figure4_transcriptome")
