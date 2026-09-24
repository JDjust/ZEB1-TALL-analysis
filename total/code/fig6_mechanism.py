"""
Figure 6 | Chromatin context and transcriptional responses to perturbation.
Perturbations differ in target, developmental background and measurement;
these panels do not establish a universal regulatory mechanism.

Panels
  A  H3K27ac HiChIP ZEB1-locus contact rate (GSE243915)
  B  Promoter CpG methylation in CIMP+ and CIMP− (GSE69954)
  C  LMO2/TAL1/LDB1/GATA2 occupancy heatmap (GSE154675; no input)
  D  TLX1 knockdown raises (not lowers) ZEB1 (GSE110635)
  E  ZEB1/Zeb1 perturbation forest: Lmo2 ON/OFF, KO and TLX1 KD
  F  CollecTRI ULM TF activity after TLX1 KD
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import _style as S

S.set_style()
C = S.C

occ = pd.read_csv(S.tbl("module6", "M6_r2_GSE154675_occupancy.tsv"), sep="\t")
audited_occ = pd.read_csv(os.path.join(S.DATA_DIR,"validation/chromatin_windows/window_summaries.tsv"),sep="\t")
occ = audited_occ[audited_occ.window.eq("legacy_12kb")].rename(columns={"mean_covered":"hg38_prom"})
occ.to_csv(os.path.join(S.DATA_DIR,"SourceData_Fig6C_exact_window_signal.tsv"),sep="\t",index=False)
kd = pd.read_csv(S.tbl("module9", "M9_keygene_effects.tsv"), sep="\t")
probes = pd.read_csv(S.tbl("module6", "M6_6.17_ZEB1_probe_CIMP_stats.tsv"), sep="\t")
pert = pd.read_csv(os.path.join(S.DATA_DIR, "SourceData_Fig6E_perturbation_forest.tsv"), sep="\t")

fig = S.new_fig(183, 178, constrained=False)
# Each row has its own label budget; the perturbation labels no longer push
# all upper panels to the right. Alphabetical reading order is A-B-C / D-E / F.
positions = {'A':[.085,.735,.225,.215], 'B':[.425,.735,.225,.215],
             'C':[.765,.735,.205,.215], 'D':[.085,.405,.225,.22],
             'E':[.565,.405,.405,.22], 'F':[.13,.075,.84,.22]}
axd = {key:fig.add_axes(rect) for key,rect in positions.items()}

# ---- A: HiChIP contact -----------------------------------------------------
axA = axd["A"]
hic = {"non-ETP\nT-ALL": [486, 443, 487, 479], "ETP": [195, 213, 224],
       "CD34+": [174, 167], "Thymus": [193]}
groups = list(hic.keys()); gcols = [C["tall"], C["ball"], C["aml"], C["grey"]]
x = np.arange(len(groups)); means = [np.mean(v) for v in hic.values()]
S.add_grid(axA, "y")
for xi, mean, col in zip(x, means, gcols):
    axA.plot([xi-.22, xi+.22], [mean, mean], color=col, lw=2, zorder=3)
rng = np.random.default_rng(0)
for xi, v in zip(x, hic.values()):
    axA.scatter(xi + rng.uniform(-0.08, 0.08, len(v)), v, s=16, color="#222", zorder=4)
axA.set_xticks(x); axA.set_xticklabels(groups, fontsize=5.5, rotation=25, ha="right")
axA.set_ylabel("ZEB1-locus contacts / million")
axA.set_title("H3K27ac HiChIP (GSE243915)", loc="left", fontsize=7.5)
axA.plot([0, 1], [520, 520], color="#333", lw=0.7)
axA.text(0.5, 528, "P=0.057", ha="center", fontsize=6)
axA.set_ylim(0, 575)
S.style_ax(axA); S.plabel(axA, "A")

# ---- B: promoter CpGs stay hypomethylated ---------------------------------
axB = axd["B"]
pr = probes[(probes["region"] == "promoter")].dropna(subset=["mean_CIMP_pos", "mean_CIMP_neg"]).copy()
S.add_grid(axB, "both")
axB.plot([0, 0.55], [0, 0.55], color="#C4C8CE", lw=0.7, ls="--", zorder=1)
axB.axhline(0.5, color="#999", lw=0.6, ls=":", zorder=1)
axB.axvline(0.5, color="#999", lw=0.6, ls=":", zorder=1)
sig = pr["q_bh"] < 0.05 if "q_bh" in pr.columns else pd.Series(False, index=pr.index)
axB.scatter(pr.loc[~sig, "mean_CIMP_neg"], pr.loc[~sig, "mean_CIMP_pos"],
            s=16, color=C["teal"], edgecolor="#333", lw=0.3, zorder=3, label="ns")
if sig.any():
    axB.scatter(pr.loc[sig, "mean_CIMP_neg"], pr.loc[sig, "mean_CIMP_pos"],
                s=22, color=C["tall"], edgecolor="#333", lw=0.35, zorder=4, label="FDR<0.05")
axB.set_xlim(0, 0.55); axB.set_ylim(0, 0.55)
axB.set_xlabel("Mean β, CIMP− (n=25)")
axB.set_ylabel("Mean β, CIMP+ (n=40)")
axB.set_title("Promoter CpG methylation", loc="left", fontsize=7.5)
axB.text(0.04, 0.96, "GSE69954  |  mean β≈0.10\ng=0.35, P=0.23",
         transform=axB.transAxes, fontsize=5.6, va="top")
axB.legend(loc="lower right", fontsize=5.4)
S.style_ax(axB); S.plabel(axB, "B")

# ---- C: occupancy heatmap --------------------------------------------------
axC = axd["C"]
abs_ = ["LMO2", "TAL1", "LDB1", "GATA2"]
lines = ["ARR", "DU528", "HSB2", "CCRFCEM"]
M = np.array([[float(occ.loc[(occ.line == ln) & (occ.antibody == ab), "hg38_prom"].iloc[0])
               if ((occ.line == ln) & (occ.antibody == ab)).any() else np.nan
               for ab in abs_] for ln in lines])
im = axC.imshow(M, cmap="YlOrRd", aspect="auto", vmin=0)
axC.set_xticks(range(len(abs_))); axC.set_xticklabels(abs_, fontsize=6.2)
axC.set_yticks(range(len(lines))); axC.set_yticklabels(lines, fontsize=6.2)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axC.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.4,
                 color="white" if M[i, j] > 1.6 else "#222")
axC.set_title("Shared-region ChIP signal", loc="left", fontsize=7.2)
cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.03)
cb.set_label("Mean covered-base signal", fontsize=5.8)
cb.ax.tick_params(labelsize=5.2)
axC.text(0.0, -0.22, "ZEB1 / ZEB1-AS1 overlap; no input", transform=axC.transAxes,
         fontsize=5.4, style="italic", color="#555")
S.plabel(axC, "C")

# ---- D: TLX1 KD ------------------------------------------------------------
axD = axd["D"]
d = kd[kd.experiment == "gse110635_combined"].set_index("gene").reindex(
    ["TLX1", "GATA3", "ZEB1", "DNTT", "IL7R"])
y = np.arange(len(d))[::-1]
cols = [C["down"] if v < 0 else C["up"] for v in d["log2FoldChange"]]
S.add_grid(axD, "x")
axD.hlines(y, 0, d["log2FoldChange"], color=cols, lw=1, zorder=2)
axD.scatter(d["log2FoldChange"], y, c=cols, s=28, zorder=3)
axD.axvline(0, color="#333", lw=0.6, zorder=2)
axD.set_yticks(y); axD.set_yticklabels(d.index, fontsize=6.4)
axD.set_xlabel("log2FC (KD vs control)")
axD.set_title("TLX1 knockdown (GSE110635)", loc="left", fontsize=7.5)
for yi, v, p in zip(y, d["log2FoldChange"], d["padj"]):
    axD.annotate(S.stars(p), (v, yi), xytext=(0,5), textcoords="offset points",
                 va="bottom", ha="center", fontsize=5.6)
axD.set_xlim(-1.6, 2.4)
axD.set_ylim(-.45,4.65)
S.style_ax(axD); S.plabel(axD, "D")

# ---- F: decoupler TF activity ----------------------------------------------
axF = axd["F"]
tf = pd.read_csv(os.path.join(S.DATA_DIR, "F6_decoupler_ulm_focus.tsv"), sep="\t").sort_values("delta_kd_minus_ctrl")
audit_dir = os.path.join(S.DATA_DIR, "validation/regulon_coverage")
sample_meta = pd.read_csv(os.path.join(audit_dir, "samples.tsv"), sep="\t")
sample_meta["group_order"] = sample_meta.group.map({"control":0,"siTLX1_1":1,"siTLX1_2":2})
sample_meta = sample_meta.sort_values(["group_order", "rep"])
activities = pd.read_csv(os.path.join(audit_dir, "sample_ulm_activities.tsv"), sep="\t", index_col=0)
matrix = activities.loc[sample_meta.gsm, tf.tf].T
ctrl = sample_meta.loc[sample_meta.group.eq("control"), "gsm"]
matrix = matrix.sub(matrix[ctrl].mean(axis=1), axis=0)
matrix.to_csv(os.path.join(S.DATA_DIR,"SourceData_Fig6F_sample_activity.tsv"),sep="\t")
limit = float(np.abs(matrix.to_numpy()).max())
imF = axF.imshow(matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
axF.set_yticks(range(len(matrix)), matrix.index, fontsize=5.5)
axF.set_xticks(range(len(sample_meta)), [f"{['Ctrl','si1','si2'][int(r.group_order)]}-{r.rep}" for r in sample_meta.itertuples()], fontsize=5.4)
for boundary in [2.5,5.5]: axF.axvline(boundary,color="white",lw=1.3)
axF.set_title("Sample-level TF activity | centered on control mean",loc="left",fontsize=7.5,pad=18)
for x0,x1,label,col in [(-.5,2.5,'Control',C['grey']),(2.5,5.5,'TLX1 siRNA 1',C['ball']),(5.5,8.5,'TLX1 siRNA 2',C['purple'])]:
    axF.plot([x0,x1],[-.9,-.9],lw=3,color=col,clip_on=False)
    axF.text((x0+x1)/2,-1.2,label,ha='center',va='bottom',fontsize=5.5,clip_on=False)
cbF = fig.colorbar(imF,ax=axF,fraction=.025,pad=.02)
cbF.set_label("ULM activity difference",fontsize=5.5)
cbF.ax.tick_params(labelsize=5)
S.plabel(axF, "F", dy=18)

# ---- E: perturbation forest (D:/ZEB1 + module9) ----------------------------
axE = axd["E"]
p = pert.copy()
p["log2FC"] = pd.to_numeric(p["log2FC"], errors="coerce")
p["ci_lo"] = pd.to_numeric(p["ci_lo"], errors="coerce")
p["ci_hi"] = pd.to_numeric(p["ci_hi"], errors="coerce")
p["fdr_num"] = pd.to_numeric(p["fdr"], errors="coerce")
p["lab"] = ["TLX1 KD · ALL-SIL", "Lmo2-TG · young DN2", "Lmo2-TG · young DN3a",
            "Lmo2-TG · old DN3a", "Lmo2-TG · DN thymocytes", "Lmo2 KO · DN1"]
assert len(p)==6, "Revisit display labels if the source comparisons change"
y = np.arange(len(p))[::-1]
S.add_grid(axE, "x")
axE.axvline(0, color="#333", lw=0.6, zorder=2)
for yi, row in zip(y, p.itertuples(index=False)):
    if yi % 2 == 0:
        axE.axhspan(yi-.45,yi+.45,color="#F5F6F8",zorder=0)
    fdr = row.fdr_num
    sig = pd.notna(fdr) and fdr < 0.05
    col = (C["up"] if row.log2FC > 0 else C["down"]) if sig else C["grey"]
    if pd.notna(row.ci_lo) and pd.notna(row.ci_hi):
        axE.plot([row.ci_lo, row.ci_hi], [yi, yi], color=col, lw=1.5, zorder=3)
    axE.plot(row.log2FC, yi, "o", color=col, ms=6.2, zorder=4,
             markeredgecolor="#333", markeredgewidth=0.35)
    note = "NA" if pd.isna(fdr) else (f"{fdr:.2g}" if fdr >= 0.001 else f"{fdr:.1e}")
    axE.text(1.03,yi,note,transform=axE.get_yaxis_transform(),fontsize=5.2,
             va="center",color=C["ink"] if sig else "#777")
axE.text(1.03,1.025,"FDR",transform=axE.transAxes,fontsize=5.2,ha="left")
axE.set_yticks(y)
axE.set_yticklabels(p["lab"], fontsize=5.8)
axE.set_xlabel("ZEB1 / Zeb1 log2FC vs control",fontsize=6.5)
axE.set_xlim(-0.7, 1.35)
axE.set_title("ZEB1 response across perturbation contexts",
              loc="left", fontsize=7.5)
S.style_ax(axE); S.plabel(axE, "E")

S.save(fig, "Figure6_mechanism")
