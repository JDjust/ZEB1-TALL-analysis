"""
Figure 1 | ZEB1 is a T-lineage gene, stably elevated in T-ALL relative to other
leukemias and normal marrow, and inherited from normal T-cell identity.

Panels
  A  Disease-spectrum ranking of ZEB1 across MILE leukemias (GSE13159)
  B  Effect sizes: T-ALL vs each comparator (Hedges g, 95% CI)
  C  Three independent cohorts + random-effects meta (T-ALL vs B-ALL)
  D  ZEB1 vs ZEB2 vs LMO2 (T-ALL vs B-ALL) diverging effect sizes
  E  Single-gene discrimination of T-ALL (ROC AUC, 95% CI)
  F  ZEB1 across 29 healthy immune cell types (GSE107011): a T/lymphoid gene
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.colors import to_rgba
import _style as S

S.set_style()
C = S.C

rank = pd.read_csv(S.tbl("module1", "M1_1.6_ZEB1_disease_ranking.tsv"), sep="\t")
pair = pd.read_csv(S.tbl("module1", "M1_1.2_1.5_pairwise_vs_TALL.tsv"), sep="\t")
coh = pd.read_csv(S.tbl("module1", "M1_1.13_cohort_effect_sizes.tsv"), sep="\t")
meta = pd.read_csv(S.tbl("module1", "M1_1.14_random_effects_meta.tsv"), sep="\t")
genes = pd.read_csv(S.tbl("module1", "M1_1.7_ZEB1_ZEB2_LMO2_effects.tsv"), sep="\t")
roc = pd.read_csv(S.tbl("module1", "M1_1.9_multiROC_three_genes.tsv"), sep="\t")
atlas = pd.read_csv(S.tbl("module1", "M1_r2_1.1_celltype_ZEB1_rank.tsv"), sep="\t")

fig = S.new_fig(183, 215, constrained=False)
gs = fig.add_gridspec(3, 1, height_ratios=[1.08, 1.0, 1.28],
                      hspace=0.60, left=0.11, right=0.965, top=0.955, bottom=0.14)
g0 = gs[0].subgridspec(1, 2, wspace=0.55, width_ratios=[1.0, 1.15])
g1 = gs[1].subgridspec(1, 3, wspace=0.55)
axA = fig.add_subplot(g0[0, 0])
axB = fig.add_subplot(g0[0, 1])
axC = fig.add_subplot(g1[0, 0])
axD = fig.add_subplot(g1[0, 1])
axE = fig.add_subplot(g1[0, 2])
axF = fig.add_subplot(gs[2])
axd = {"A": axA, "B": axB, "C": axC, "D": axD, "E": axE, "F": axF}

# ---- A: disease ranking ----------------------------------------------------
axA = axd["A"]
r = rank.sort_values("mean_z")
cols = [C["tall"] if d == "T-ALL" else C["grey"] for d in r["disease"]]
S.add_grid(axA, "x")
axA.barh(r["disease"], r["mean_z"], color=cols, edgecolor="white", linewidth=0.4, height=0.72, zorder=3)
axA.axvline(0, color="#333", lw=0.6, zorder=2)
axA.set_xlabel("ZEB1 standardized expression (mean z)")
axA.set_title("Disease-spectrum ranking (MILE, n=2096)", loc="left", fontsize=7.5)
axA.set_xlim(-0.78, 1.2)
for yi, (d, v) in enumerate(zip(r["disease"], r["mean_z"])):
    axA.text(v + (0.03 if v >= 0 else -0.03), yi, f"{v:+.2f}", va="center",
             ha="left" if v >= 0 else "right", fontsize=5.8)
S.style_ax(axA); S.plabel(axA, "A")

# ---- B: forest vs comparators ---------------------------------------------
axB = axd["B"]
b = pair[pair.gene == "ZEB1"].copy()
b["lab"] = b["group_b"]
b = b.iloc[::-1].reset_index(drop=True)
y = np.arange(len(b))
S.add_grid(axB, "x")
axB.errorbar(b["hedges_g"], y, xerr=[b["hedges_g"] - b["ci95_lo"], b["ci95_hi"] - b["hedges_g"]],
             fmt="o", color=C["tall"], ecolor=C["tall"], ms=5, lw=1.2, capsize=2.2, zorder=3)
axB.axvline(0, color="#333", lw=0.6, ls="--", zorder=2)
axB.set_yticks(y); axB.set_yticklabels(b["lab"])
axB.set_xlabel("Hedges g (T-ALL higher →)")
axB.set_title("T-ALL vs comparator (GSE13159)", loc="left", fontsize=7.5)
axB.set_xlim(-0.2, 3.55)
for yi, g, hi in zip(y, b["hedges_g"], b["ci95_hi"]):
    S.text_right(axB, hi, yi, f"g={g:.2f}  n={int(b['n_b'].iloc[yi])}", pad=0.08, fontsize=5.4)
S.style_ax(axB); S.plabel(axB, "B")

# ---- C: three cohorts + meta ----------------------------------------------
axC = axd["C"]
c = coh[coh.gene == "ZEB1"].copy()
labels = list(c["cohort"]) + ["Random-effects meta"]
g = list(c["hedges_g"]) + [meta["hedges_g"].iloc[0]]
lo = list(c["ci95_lo"]) + [meta["ci95_lo"].iloc[0]]
hi = list(c["ci95_hi"]) + [meta["ci95_hi"].iloc[0]]
ns = list(c["n_a"]) + [""]
y = np.arange(len(labels))[::-1]
cols = [C["ball"]] * len(c) + [C["tall"]]
S.add_grid(axC, "x")
for yi, gi, loi, hii, col in zip(y, g, lo, hi, cols):
    axC.plot([loi, hii], [yi, yi], color=col, lw=1.4, zorder=3)
    axC.plot(gi, yi, "D" if col == C["tall"] else "o", color=col,
             ms=7 if col == C["tall"] else 5, zorder=4)
axC.axvline(0, color="#333", lw=0.6, ls="--", zorder=2)
axC.set_yticks(y)
axC.set_yticklabels([f"{l}" + (f" (T={n})" if n != "" else "") for l, n in zip(labels, ns)], fontsize=6)
axC.set_xlabel("Hedges g, ZEB1 (T-ALL vs B-ALL)")
axC.set_title("Independent replication + meta", loc="left", fontsize=7.5)
axC.set_xlim(0, 2.85)
mg = meta["hedges_g"].iloc[0]
S.text_right(axC, meta["ci95_hi"].iloc[0], y[-1],
             f"pooled g={mg:.2f}, I²={meta['I2'].iloc[0]:.0f}%", pad=0.06, fontsize=5.6)
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: ZEB1/ZEB2/LMO2 -----------------------------------------------------
axD = axd["D"]
d = genes[genes.contrast == "T-ALL vs B-ALL"].set_index("gene").loc[["ZEB1", "ZEB2", "LMO2"]]
y = np.arange(len(d))[::-1]
cols = [C["up"] if v > 0 else C["down"] for v in d["hedges_g"]]
S.add_grid(axD, "x")
axD.barh(y, d["hedges_g"], color=cols, height=0.62, edgecolor="white", linewidth=0.4, zorder=3)
axD.errorbar(d["hedges_g"], y, xerr=[d["hedges_g"] - d["ci95_lo"], d["ci95_hi"] - d["hedges_g"]],
             fmt="none", ecolor="#333", lw=0.8, capsize=2.2, zorder=4)
axD.axvline(0, color="#333", lw=0.6, zorder=2)
axD.set_yticks(y); axD.set_yticklabels(d.index)
axD.set_xlabel("Hedges g (T-ALL vs B-ALL)")
axD.set_title("ZEB axis", loc="left", fontsize=7.5)
S.style_ax(axD); S.plabel(axD, "D")

# ---- E: ROC AUC ------------------------------------------------------------
axE = axd["E"]
rr = roc.set_index("gene").loc[["ZEB2", "ZEB1", "LMO2"]]
y = np.arange(len(rr))
S.add_grid(axE, "x")
axE.errorbar(rr["auc"], y, xerr=[rr["auc"] - rr["ci_lo"], rr["ci_hi"] - rr["auc"]],
             fmt="o", color=C["purple"], ms=5, lw=1.2, capsize=2.2, zorder=3)
axE.axvline(0.5, color="#999", lw=0.7, ls=":", zorder=2)
axE.set_yticks(y)
axE.set_yticklabels([f"{gname} ({dr})" for gname, dr in zip(rr.index, rr["direction"])], fontsize=6)
axE.set_xlabel("AUC (discriminate T-ALL)")
axE.set_xlim(0.42, 1.05)
axE.set_title("Single-gene ROC", loc="left", fontsize=7.5)
for yi, a, hi in zip(y, rr["auc"], rr["ci_hi"]):
    S.text_right(axE, hi, yi, f"{a:.2f}", pad=0.015)
S.style_ax(axE); S.plabel(axE, "E")

# ---- F: immune atlas -------------------------------------------------------
axF = axd["F"]
scores = pd.read_csv(S.tbl("module1", "M1_r2_GSE107011_scores.tsv"), sep="\t")
marker_genes = ["ZEB1", "ZEB2", "LMO2", "LYL1", "TCF7", "GATA3", "BCL11B",
                "MEF2C", "IL7R", "CD3D", "CD19", "NKG7", "CD14", "CD34"]
# One donor-celltype mean per unit; genes selected by role, not their ZEB1 correlation.
donor = scores.groupby(["lineage", "cell_raw", "donor"], observed=True)[marker_genes].mean()
means = donor.groupby(level=["lineage", "cell_raw"]).mean()
z = means.sub(means.mean(axis=0), axis=1).div(means.std(axis=0, ddof=0), axis=1)
assert np.isfinite(z.to_numpy()).all()
ordered = []
blocks = []
for lineage in ["T", "B", "NK", "Myeloid", "Progenitor", "Other"]:
    labels = [idx for idx in z.index if idx[0] == lineage]
    if not labels:
        continue
    if len(labels)>1:
        tree = linkage(z.loc[labels].to_numpy(), method="average", metric="euclidean", optimal_ordering=True)
        labels = [labels[i] for i in leaves_list(tree)]
    blocks.append((lineage, len(ordered), len(labels)))
    ordered.extend(labels)
matrix = z.loc[ordered, marker_genes].T
means.loc[ordered].to_csv(S.DATA_DIR + "/SourceData_Fig1F_donor_mean_expression.tsv", sep="\t")
matrix.to_csv(S.DATA_DIR + "/SourceData_Fig1F_display_zscores.tsv", sep="\t")
donor.reset_index().to_csv(S.DATA_DIR + "/SourceData_Fig1F_donor_celltype_expression.tsv", sep="\t", index=False)
limit = 2.5  # Display saturation only; exported values retain their full magnitude.
imF = axF.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit, interpolation="nearest")
axF.set_xticks(range(len(ordered)), [idx[1].replace("_", " ") for idx in ordered], rotation=90, fontsize=5.3)
axF.set_yticks(range(len(marker_genes)), marker_genes, fontsize=5.5)
axF.get_yticklabels()[0].set_color(C["tall"])
axF.get_yticklabels()[0].set_fontweight("bold")
axF.add_patch(Rectangle((-.5,-.5),len(ordered),1,fill=False,edgecolor=C["tall"],lw=.85))
strip = axF.inset_axes([0, 1.025, 1, .045])
strip.imshow(np.array([[to_rgba(S.LIN.get(idx[0],C['grey'])) for idx in ordered]]),aspect="auto")
strip.set_axis_off()
for lineage,start,n in blocks:
    short = {"Progenitor":"P", "Other":"O"}.get(lineage,lineage)
    axF.text((start+(n-1)/2+.5)/len(ordered),1.085,short,transform=axF.transAxes,
             ha="center",va="bottom",fontsize=5.3,color=S.LIN.get(lineage,C['grey']))
    if start: axF.axvline(start-.5,color="white",lw=1.2)
axF.set_title("Healthy immune atlas: lineage programs and ZEB1 context",loc="left",fontsize=7.5,pad=27)
cax = axF.inset_axes([1.012, 0, .014, 1])
cbF = fig.colorbar(imF,cax=cax,extend="both")
cbF.set_label("Within-gene z score",fontsize=5.5)
cbF.ax.tick_params(labelsize=5)
S.plabel(axF, "F")

S.save(fig, "Figure1_disease_specificity")
