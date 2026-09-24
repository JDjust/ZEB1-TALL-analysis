"""
Figure 7 | Drug-response associations, genetic dependency and clinical MRD.
All displayed estimates are reused from existing exports.

Panels
  A  17 analyzable drugs from an 18-drug ex-vivo panel
  B  MRD
  C  ZEB1 CRISPR gene effect heterogeneity
  D  Selected target gene effects
  E  T-versus-B mean gene-effect differences with uncertainty
  F  GDSC1/GDSC2 MTX in T-ALL, with shared models identified
  G  Paired patient MRD transitions and incremental-model AUC
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.ticker import ScalarFormatter, NullFormatter
import _style as S

S.set_style()
C = S.C

drug = pd.read_csv(S.tbl("module7", "M7_7.7_7.8_drug_spearman.tsv"), sep="\t")
mrd = pd.read_csv(S.tbl("module7", "M7_7.18_ZEB1_MRD.tsv"), sep="\t")
dep = pd.read_csv(os.path.join(S.DATA_DIR, "SourceData_Fig7C_ZEB1_dependency.tsv"), sep="\t")
tgt = pd.read_csv(os.path.join(S.DATA_DIR, "SourceData_Fig7DE_drug_target_GE.tsv"), sep="\t")
assert dep["ZEB1_GE"].notna().sum() == 8 and tgt["n_TALL"].eq(8).all()

fig = S.new_fig(183, 190, constrained=False)
positions = {'A':[.16,.735,.46,.22], 'B':[.80,.735,.16,.22],
             'C':[.16,.415,.20,.215], 'D':[.48,.415,.19,.215],
             'E':[.80,.415,.16,.215], 'F':[.16,.085,.36,.225],
             'G':[.655,.075,.32,.24]}
axd = {key:fig.add_axes(rect) for key,rect in positions.items()}

# ---- A: 18-drug lollipop ---------------------------------------------------
axA = axd["A"]
z = drug[drug.gene == "ZEB1_log"].sort_values("rho")
assert len(z) == 17 and z["n"].min() == 14 and z["n"].max() == 94
y = np.arange(len(z))
cols = [C["tall"] if f < 0.05 else C["grey"] for f in z["fdr"]]
S.add_grid(axA, "x")
axA.hlines(y, 0, z["rho"], color="#C7CBD1", lw=1.0, zorder=2)
axA.scatter(z["rho"], y, s=20, color=cols,
            edgecolor="#333", linewidth=0.4, zorder=3)
axA.axvline(0, color="#333", lw=0.6, zorder=2)
axA.set_yticks(y); axA.set_yticklabels(z["drug_short"], fontsize=5.5)
axA.set_xlabel("Spearman ρ (ZEB1 vs LC50)")
axA.set_title("17 analyzable drugs; MTX absent", loc="left", fontsize=7.5)
for yi,n in zip(y,z['n']):
    axA.text(1.025,yi,str(int(n)),transform=axA.get_yaxis_transform(),fontsize=5.2,va='center')
axA.text(1.025,1.025,'n',transform=axA.transAxes,fontsize=5.5)
axA.text(.01,-.25,"0/17 at FDR < 0.05; sample size shown per drug",transform=axA.transAxes,fontsize=5.4)
S.style_ax(axA); S.plabel(axA, "A")

# ---- B: MRD ----------------------------------------------------------------
axB = axd["B"]
x = np.arange(len(mrd))
cols = [C["down"] if p < 0.05 else C["grey"] for p in mrd["p"]]
S.add_grid(axB, "y")
axB.vlines(x, 0, mrd["rho"], color="#C5CAD1", lw=1.4, zorder=2)
axB.scatter(x, mrd["rho"], color=cols, s=43, edgecolor="#333",
            linewidth=.4, zorder=3)
axB.axhline(0, color="#333", lw=0.6, zorder=2)
axB.set_xticks(x); axB.set_xticklabels(["Mid-induction\nMRD", "End-induction\nMRD"], fontsize=6.2)
axB.set_ylabel("Spearman ρ with ZEB1")
axB.set_title("Residual disease", loc="left", fontsize=7.5)
for xi, r, p, n in zip(x, mrd["rho"], mrd["p"], mrd["n"]):
    axB.text(xi, r - 0.015, f"{S.stars(p)}\nn={n}", ha="center", va="top", fontsize=5.8)
axB.set_ylim(-0.36, 0.08)
S.style_ax(axB); S.plabel(axB, "B")

# ---- C: ZEB1 self dependency -----------------------------------------------
axC = axd["C"]
d = dep.dropna(subset=["ZEB1_GE"]).copy()
d["ZEB1_GE"] = d["ZEB1_GE"].astype(float)
d = d.sort_values("ZEB1_GE")
y = np.arange(len(d))
cols = [C["tall"] if v < -0.5 else C["grey"] for v in d["ZEB1_GE"]]
S.add_grid(axC, "x")
axC.hlines(y, d["ZEB1_GE"], 0, color="#C5CAD1", lw=1.1, zorder=2)
axC.scatter(d["ZEB1_GE"], y, color=cols, s=31, edgecolor="#333",
            linewidth=.35, zorder=3)
axC.axvline(-0.5, color=C["down"], lw=0.9, ls="--", zorder=2)
axC.axvline(0, color="#333", lw=0.6, zorder=2)
axC.set_yticks(y); axC.set_yticklabels(d["CellLineName"], fontsize=5.6)
axC.set_xlabel("ZEB1 CRISPR gene effect")
axC.set_title("ZEB1 dependence varies (n=8)", loc="left", fontsize=7.5)
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: drug target essentiality -------------------------------------------
axD = axd["D"]
targets = ["BCL2L1", "PSMB5", "CDK6", "MCL1", "HDAC3", "JAK1", "BCL2"]
t = tgt.set_index("gene").reindex(targets)
y = np.arange(len(t))[::-1]
S.add_grid(axD, "x")
for yi, (_, row) in zip(y, t.iterrows()):
    axD.plot([row.mean_TALL, row.mean_BALL], [yi, yi], color="#C5CAD1",
             lw=1.4, zorder=2)
axD.scatter(t["mean_BALL"], y, color=C["ball"], s=23, edgecolor="white",
            linewidth=.3, zorder=3, label="B-ALL (n=13)")
axD.scatter(t["mean_TALL"], y, color=C["purple"], s=29, edgecolor="white",
            linewidth=.3, zorder=4, label="T-ALL (n=8)")
axD.axvline(-0.5, color=C["down"], lw=0.9, ls="--", zorder=2)
axD.axvline(0, color="#333", lw=0.6, zorder=2)
axD.set_yticks(y); axD.set_yticklabels(t.index, fontsize=6.2)
axD.set_xlabel("Mean CRISPR gene effect")
axD.set_title("Selected target effects", loc="left", fontsize=7.0)
axD.legend(loc="upper center",bbox_to_anchor=(.5,1.01),fontsize=4.8,frameon=False,ncol=2,
           handletextpad=.2,columnspacing=.4)
axD.set_ylim(-.5,7.2)
S.style_ax(axD); S.plabel(axD, "D")

# ---- E: lineage specificity ------------------------------------------------
axE = axd["E"]
tfs = ["GATA3", "BCL11B", "NOTCH1", "ZEB1", "LMO2", "ZEB2"]
e = tgt.set_index("gene").reindex(tfs)
contrasts = pd.read_csv(os.path.join(S.DATA_DIR, "SourceData_Fig7E_lineage_contrasts.tsv"), sep="\t")
e = e.join(contrasts.set_index("gene")[["bootstrap_ci95_lo", "bootstrap_ci95_hi", "fdr_bh_32_genes"]])
y = np.arange(len(e))[::-1]
cols = [C["tall"] if q < 0.05 else C["grey"] for q in e["fdr_bh_32_genes"]]
S.add_grid(axE, "x")
for yi, (_, row), col in zip(y, e.iterrows(), cols):
    axE.plot([row.bootstrap_ci95_lo, row.bootstrap_ci95_hi], [yi, yi],
             color=col, lw=1.5, zorder=3)
    axE.plot(row.delta_T_minus_B, yi, "o", color=col, ms=4.5, zorder=4)
axE.axvline(0, color="#333", lw=0.6, zorder=2)
axE.set_yticks(y); axE.set_yticklabels(e.index, fontsize=6.2)
axE.set_xlabel("Δ gene effect (T − B)")
axE.set_xlim(-1.22, 0.88)
axE.set_title("T–B differences", loc="left", fontsize=7.0)
axE.text(.02,-.25,"95% bootstrap CI\nBH FDR across 32 genes",transform=axE.transAxes,fontsize=4.8,va='top')
S.style_ax(axE); S.plabel(axE, "E")

# ---- F: direct MTX signal in the two overlapping GDSC screens ------------
axF = axd["F"]
mtx = pd.read_csv(os.path.join(S.DATA_DIR, "F7_mtx_points.tsv"), sep="\t")
sens = pd.read_csv(os.path.join(S.DATA_DIR, "validation", "mtx_screen_sensitivity.tsv"), sep="\t")
overlap = pd.read_csv(os.path.join(S.DATA_DIR, "validation", "mtx_model_overlap.tsv"), sep="\t")
S.add_grid(axF, "both")
shared = overlap.dropna(subset=["GDSC1", "GDSC2"])
assert len(shared) == 8
for row in shared.itertuples(index=False):
    axF.plot([row.ZEB1, row.ZEB1], [row.GDSC1, row.GDSC2],
             color="#B8BCC2", lw=0.7, alpha=0.75, zorder=2)
for screen, color, mark in [("GDSC1", C["ball"], "s"), ("GDSC2", C["tall"], "o")]:
    t = mtx[(mtx.source == screen) & (mtx.lineage == "T-ALL")].dropna(subset=["ZEB1", "value"])
    r = sens.set_index("screen").loc[screen]
    axF.scatter(t.ZEB1, t.value, s=34, color=color, marker=mark,
                edgecolor="#222", lw=0.35, zorder=4,
                label=f"{screen}: n={int(r.n_models)}, ρ={r.rho:.2f}")
axF.set_xlabel("ZEB1 RNA (log2 TPM+1)")
axF.set_ylabel("MTX LN_IC50")
axF.set_title("Methotrexate response | T-ALL models", loc="left", fontsize=7.0)
axF.text(0.03, 0.05, "8 shared ModelIDs; exploratory", transform=axF.transAxes,
         fontsize=5.2, va="bottom")
axF.legend(loc="upper left", fontsize=5.4)
S.style_ax(axF); S.plabel(axF, "F")

# ---- G: paired patient MRD transitions (provisional threshold) ------------
axG = axd["G"]
clinical_dir = os.path.join(S.DATA_DIR, "clinical_mrd_persistence_20260923")
pred = pd.read_csv(os.path.join(clinical_dir, "paired_patients_predictions.tsv"), sep="\t")
flow = pd.read_csv(os.path.join(S.DATA_DIR, "validation", "mrd_transitions.tsv"), sep="\t")
with open(os.path.join(clinical_dir, "summary.json"), encoding="utf-8") as fh:
    clinical = json.load(fh)
pred.to_csv(os.path.join(S.DATA_DIR, "SourceData_Fig7G_MRD_persistence.tsv"), sep="\t", index=False)
counts = {(r.day15, r.later): int(r.n_patients) for r in flow.itertuples(index=False)}
assert counts == {("below_threshold", "below_threshold"): 26,
                  ("below_threshold", "positive"): 0,
                  ("positive", "below_threshold"): 48,
                  ("positive", "positive"): 18}
assert clinical["n_eligible"] == 66 and clinical["n_late_positive"] == 18

def ribbon(y0lo, y0hi, y1lo, y1hi, color):
    x0, x1 = .20, .80
    verts = [(x0, y0lo), (.43, y0lo), (.57, y1lo), (x1, y1lo),
             (x1, y1hi), (.57, y1hi), (.43, y0hi), (x0, y0hi), (x0, y0lo)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    axG.add_patch(PathPatch(Path(verts, codes), fc=color, ec="none", alpha=.42,
                            transform=axG.transAxes, zorder=2))

bottom, scale = .14, .69 / 92
coord = lambda n: bottom + n * scale
# Right negative group is ordered by source to expose both clinical routes.
ribbon(coord(0), coord(26), coord(48), coord(74), C["ball"])
ribbon(coord(26), coord(74), coord(0), coord(48), C["teal"])
ribbon(coord(74), coord(92), coord(74), coord(92), C["tall"])
for x, spans in [(.16, [(0, 26, C["ball"]), (26, 92, C["tall"])]),
                 (.80, [(0, 74, C["teal"]), (74, 92, C["tall"])])]:
    for lo, hi, col in spans:
        axG.add_patch(Rectangle((x, coord(lo)), .04, (hi-lo)*scale,
                                transform=axG.transAxes, fc=col, ec="white", lw=.3, zorder=4))
axG.text(.16, .94, "Mid-induction", transform=axG.transAxes, ha="left", fontsize=6.2)
axG.text(.84, .94, "End-induction", transform=axG.transAxes, ha="right", fontsize=6.2)
axG.text(.13, coord(13), "<0.01%\n26", transform=axG.transAxes, ha="right", va="center", fontsize=5.0)
axG.text(.13, coord(59), "≥0.01%\n66", transform=axG.transAxes, ha="right", va="center", fontsize=5.0)
axG.text(.87, coord(37), "<0.01%\n74", transform=axG.transAxes, ha="left", va="center", fontsize=5.0)
axG.text(.87, coord(83), "≥0.01%\n18", transform=axG.transAxes, ha="left", va="center", fontsize=5.0)
axG.text(.5, -.02, "66 mid-induction-positive patients\nLOO AUC: MRD alone → MRD + ZEB1\n"
         f"{clinical['loo_auc_baseline']:.2f}→{clinical['loo_auc_plus_zeb1']:.2f}",
         transform=axG.transAxes, ha="center", va="top",fontsize=5.3)
axG.set_xlim(0, 1); axG.set_ylim(0, 1)
axG.set_axis_off()
axG.set_title("MRD trajectories | 92 patients", loc="left", fontsize=7.0)
S.style_ax(axG); S.plabel(axG, "G")

S.save(fig, "Figure7_therapy")
