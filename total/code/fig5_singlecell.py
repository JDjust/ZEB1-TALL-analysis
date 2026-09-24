"""
Figure 5 | Patient-level single-cell reanalysis and cell-composition sensitivity.

Panels
  A  Diagnosis composition (GSE227122)
  B  Discovery patient-level associations (n=10, all FDR > 0.8)
  C  Author-annotated malignant-cell raw-UMI pseudobulk (GSE248287, n=15)
  D  Mapping to normal thymus
  E  ZEB1 vs CD8SP-like mapping (descriptive)
  F  Diagnosis-to-EOI blast clearance
  G  Mixed-cell estimates versus malignant-only estimates; no meta-analysis
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import _style as S

S.set_style()
C = S.C

comp = pd.read_csv(S.tbl("module5", "M5_5.13_composition.tsv"), sep="\t")
a10 = pd.read_csv(S.tbl("module5", "M5_5.41_patient_ZEB1_associations.tsv"), sep="\t")
a15_old = pd.read_csv(S.tbl("module5", "M5_5.42_ZEB1_patient_assoc.tsv"), sep="\t")
a15 = pd.read_csv(os.path.join(S.DATA_DIR, "gse248287_malignant_job1699",
                                "zeb1_patient_associations.tsv"), sep="\t")
mp = pd.read_csv(S.tbl("module5", "M5_r2_map_patient.tsv"), sep="\t")
a15.to_csv(os.path.join(S.DATA_DIR, "SourceData_Fig5C_malignant_assoc.tsv"), sep="\t", index=False)

fig = S.new_fig(183, 155)
axd = fig.subplot_mosaic(
    [["A", "B", "C"],
     ["D", "E", "F"],
     ["G", "G", "G"]],
    gridspec_kw=dict(height_ratios=[1.05, 1.0, 0.65]),
)

# ---- A: composition --------------------------------------------------------
axA = axd["A"]
dx = comp[comp.timepoint == "Dx"].copy()
piv = dx.pivot_table(index="patient", columns="lineage", values="frac", fill_value=0)
order_p = [f"T{i}" for i in range(1, 11) if f"T{i}" in piv.index]
piv = piv.reindex(order_p)
lin_order = ["T", "B", "NK", "Myeloid", "Ery", "Unknown"]
lin_col = {"T": C["tall"], "B": C["ball"], "NK": C["green"], "Myeloid": C["aml"],
           "Ery": C["purple"], "Unknown": C["grey"]}
bottom = np.zeros(len(piv)); x = np.arange(len(piv))
for l in lin_order:
    if l in piv.columns:
        axA.bar(x, piv[l].values, bottom=bottom, color=lin_col[l], width=0.82,
                edgecolor="white", linewidth=0.3, label=l, zorder=3)
        bottom += piv[l].values
axA.set_xticks(x); axA.set_xticklabels(order_p, rotation=45, fontsize=5.6, ha="right")
axA.set_ylabel("Cell fraction"); axA.set_ylim(0, 1)
axA.set_title("Diagnosis cell composition", loc="left", fontsize=7,pad=25)
axA.legend(ncol=3, fontsize=5.0, loc="lower left",bbox_to_anchor=(0,1),frameon=False,
           handletextpad=.3,columnspacing=.7,borderpad=.1)
S.style_ax(axA); S.plabel(axA, "A",dy=25)

# ---- B: n=10 null forest ---------------------------------------------------
axB = axd["B"]
b = a10.dropna(subset=["rho"]).copy()
b["absr"] = b["rho"].abs()
b = b.sort_values("absr", ascending=False).head(10).sort_values("rho")
b["lab"] = b["feature"].str.replace("sig_", "", regex=False).str.replace("g_", "", regex=False)
y = np.arange(len(b))
S.add_grid(axB, "x")
axB.scatter(b["rho"], y, s=26, color=C["grey"], edgecolor="#333", linewidth=0.4, zorder=3)
axB.axvline(0, color="#333", lw=0.6, ls="--", zorder=2)
axB.set_yticks(y); axB.set_yticklabels(b["lab"], fontsize=6)
axB.set_xlabel("Spearman ρ with ZEB1")
axB.set_title("Patient associations (n=10)", loc="left", fontsize=7)
axB.set_xlim(-0.85, 0.85)
axB.text(0.97, 0.06, "all FDR > 0.8", transform=axB.transAxes, fontsize=6,
         ha="right", va="bottom")
S.style_ax(axB); S.plabel(axB, "B")

# ---- C: n=15 malignant-cell analysis --------------------------------------
axC = axd["C"]
c = a15.copy()
c["lab"] = c["feature"]
c = c.sort_values("rho")
y = np.arange(len(c))
cols = [C["up"] if (r > 0 and q < 0.05) else (C["down"] if (r < 0 and q < 0.05) else C["ns"])
        for r, q in zip(c["rho"], c["fdr_bh"])]
S.add_grid(axC, "x")
axC.hlines(y, 0, c["rho"], color=cols, lw=.8, zorder=2)
axC.scatter(c["rho"], y, c=cols, s=24, zorder=3)
axC.axvline(0, color="#333", lw=0.6, zorder=2)
axC.set_yticks(y); axC.set_yticklabels(c["lab"], fontsize=6)
axC.set_xlabel("Spearman ρ with ZEB1")
axC.set_title("Malignant-cell estimates (n=15)", loc="left", fontsize=7)
axC.set_xlim(-1.18, 1.18)
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: patient-by-stage mapping heatmap ----------------------------------
axD = axd["D"]
stage_cols = ["p_DN1", "p_DN2", "p_DN3", "p_ISP", "p_DP_CD3neg", "p_DP_CD3pos", "p_CD4SP", "p_CD8SP"]
mpp = mp.set_index("patient").reindex(order_p)
stage_mat = mpp[stage_cols].to_numpy(dtype=float).T
assert stage_mat.shape == (8, len(order_p)) and np.isfinite(stage_mat).all()
imD = axD.imshow(stage_mat, cmap="YlGnBu", vmin=0, vmax=max(.65, float(stage_mat.max())),
                 aspect="auto", interpolation="nearest")
axD.set_xticks(np.arange(len(order_p)), order_p, fontsize=5.5)
axD.set_yticks(np.arange(len(stage_cols)),
               ["DN1", "DN2", "DN3", "ISP", "DP CD3-", "DP CD3+", "CD4SP", "CD8SP"],
               fontsize=5.5)
axD.set_xticks(np.arange(-.5, len(order_p), 1), minor=True)
axD.set_yticks(np.arange(-.5, len(stage_cols), 1), minor=True)
axD.grid(which="minor", color="white", linewidth=.7)
axD.tick_params(axis="both", which="both", length=0)
axD.set_title("Patient-to-thymus mapping", loc="left", fontsize=7)
cbD = fig.colorbar(imD, ax=axD, fraction=.035, pad=.02)
cbD.set_label("Probability", fontsize=5.6)
cbD.ax.tick_params(labelsize=5, length=2)
S.plabel(axD, "D")

# ---- E: ZEB1 vs CD8SP ------------------------------------------------------
axE = axd["E"]
S.add_grid(axE, "both")
axE.scatter(mpp["g_ZEB1"], mpp["p_CD8SP"], s=38, color=C["tall"], edgecolor="#333",
            linewidth=0.5, zorder=4)
S.conf_ellipse(axE, mpp["g_ZEB1"], mpp["p_CD8SP"], n_std=1.96, ec=C["tall"], lw=0.8, ls="--")
z = np.polyfit(mpp["g_ZEB1"], mpp["p_CD8SP"], 1)
xs = np.linspace(mpp["g_ZEB1"].min(), mpp["g_ZEB1"].max(), 20)
axE.plot(xs, np.polyval(z, xs), color="#333", lw=0.9, ls="--", zorder=3)
axE.set_xlabel("ZEB1 mean (blasts)")
axE.set_ylabel("P(CD8SP-like)")
axE.set_title("ZEB1 vs maturity mapping", loc="left", fontsize=7.5)
axE.text(0.04, 0.07, r"ρ=−0.75, P=0.013" "\n(FDR ns, n=10)",
         transform=axE.transAxes, fontsize=6)
S.style_ax(axE); S.plabel(axE, "E")

# ---- F: Dx -> EOI ----------------------------------------------------------
axF = axd["F"]
paired = comp[comp.lineage == "T"].pivot_table(index="patient", columns="timepoint",
                                               values="frac", fill_value=np.nan)
paired = paired.dropna(subset=["Dx", "EOI"]) if "EOI" in paired.columns else paired
S.add_grid(axF, "y")
for p in paired.index:
    axF.plot([0, 1], [paired.loc[p, "Dx"], paired.loc[p, "EOI"]], "-o",
             color=C["ball"], ms=4.5, lw=1.1, alpha=0.85, zorder=3)
axF.set_xticks([0, 1]); axF.set_xticklabels(["Diagnosis", "EOI"])
axF.set_ylabel("Blast (T) fraction")
axF.set_title("Blast clearance at EOI", loc="left", fontsize=7.5)
axF.set_xlim(-0.3, 1.3); axF.set_ylim(0, 1.08)
axF.text(0.5, -0.23, "Paired ZEB1: P=0.31; n=5",
         transform=axF.transAxes, fontsize=5.5, ha="center")
S.style_ax(axF); S.plabel(axF, "F")

# ---- G: mixed-cell versus malignant-only estimates ------------------------
axG = axd["G"]
old = a15_old.copy()
old["feature"] = old["feature"].str.replace("sig_", "", regex=False).str.replace("g_", "", regex=False)
m = old.merge(a15[["feature", "rho", "fdr_bh"]], on="feature", suffixes=("_mixed", "_malignant"))
m = m.sort_values("rho_mixed").reset_index(drop=True)
m.to_csv(os.path.join(S.DATA_DIR, "SourceData_Fig5G_composition_sensitivity.tsv"), sep="\t", index=False)
y = np.arange(len(m))
S.add_grid(axG, "x")
axG.axvline(0, color="#333", lw=0.6, zorder=2)
axG.hlines(y, m["rho_mixed"], m["rho_malignant"], color="#B7BEC8", lw=1.2, zorder=2)
axG.scatter(m["rho_mixed"], y, s=30, c=C["ball"], edgecolor="#333", lw=0.35, zorder=4,
            label="all cell types (previous)")
axG.scatter(m["rho_malignant"], y, s=30, c=C["tall"], edgecolor="#333", lw=0.35, zorder=4,
            label="author malignant cells, raw-UMI pseudobulk")
axG.set_yticks(y)
axG.set_yticklabels(m["feature"], fontsize=6.2)
axG.set_xlabel("Spearman ρ with patient ZEB1 (GSE248287, n=15)")
axG.set_xlim(-1.05, 1.05)
axG.set_title("Cell composition changes the apparent validation signal", loc="left", fontsize=7.5)
axG.legend(loc="lower left",bbox_to_anchor=(0,1.0),ncol=2,fontsize=5.2,
           handletextpad=.3,columnspacing=1)
axG.set_title("Cell-composition sensitivity",loc="left",fontsize=7,pad=18)
S.style_ax(axG); S.plabel(axG, "G",dy=18)

S.save(fig, "Figure5_singlecell")
