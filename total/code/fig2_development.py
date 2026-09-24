"""
Figure 2 | Developmental-stage modulation and sensitivity to maturity-score definition.

Panels
  A  Bulk sorted thymocyte trajectories (GSE142522)
  B  FACS scRNA thymus (GSE195812): ZEB1 by stage
  C  OLS disease coefficient under two maturity definitions
  D  Cohort-specific ZEB1-maturity correlations
  E  T-ALL maps to a bimodal nearest normal stage
  F  Stage-resolved expression and existing stage-order correlations
  G  Independent thymus atlas, donor-stage means
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import _style as S

S.set_style()
C = S.C

traj = pd.read_csv(S.tbl("module2", "M2_stage_keygene_rho.tsv"), sep="\t")
ols = pd.read_csv(S.tbl("module2", "M2_2.16_2.17_regressions.tsv"), sep="\t")
rep = pd.read_csv(S.tbl("module2", "M2_2.18_replication.tsv"), sep="\t")
near = pd.read_csv(S.tbl("module2", "M2_2.11_nearest_stage.tsv"), sep="\t")
sc = pd.read_csv(S.tbl("module2", "M2_r2_GSE195812_stage_gene_means.tsv"), sep="\t")
scr = pd.read_csv(S.tbl("module2", "M2_r2_GSE195812_gene_vs_stage.tsv"), sep="\t")

fig = S.new_fig(183, 155)
axd = fig.subplot_mosaic(
    [["A", "B", "C"],
     ["D", "E", "F"],
     ["G", "G", "G"]],
    gridspec_kw=dict(height_ratios=[1.05, 1.05, 0.75]),
)

# ---- A: bulk trajectories --------------------------------------------------
axA = axd["A"]
stages = ["cd34", "isp", "ec", "lc", "sp4", "sp8"]
slab = ["CD34", "ISP", "EC", "LC", "SP4", "SP8"]
t = traj.set_index("gene")
show = {"ZEB1": C["tall"], "ZEB2": C["ball"], "LMO2": C["aml"], "IL7R": C["green"]}
S.add_grid(axA, "y")
for g, col in show.items():
    ys = t.loc[g, stages].values.astype(float)
    axA.plot(range(6), ys, "-o", color=col, ms=3.5, lw=1.5, zorder=3)
    if g == "IL7R":
        axA.text(3.12, ys[3] + 0.12, g, color=col, fontsize=6.2, va="bottom")
    elif g == "ZEB1":
        axA.text(5.08, ys[-1] + 0.12, g, color=col, fontsize=6.2, va="bottom")
    elif g == "ZEB2":
        axA.text(5.08, ys[-1] - 0.18, g, color=col, fontsize=6.2, va="top")
    else:
        axA.text(5.08, ys[-1] + 0.12, g, color=col, fontsize=6.2, va="bottom")
axA.set_xticks(range(6)); axA.set_xticklabels(slab)
axA.set_ylabel("log2 expression")
axA.set_title("Sorted thymocytes (GSE142522)", loc="left", fontsize=7.5)
axA.set_xlim(-0.15, 5.85)
axA.set_ylim(5.3, 14.1)
axA.text(0.03, 0.97, r"ZEB1 range 1.34, ρ=−0.43 (ns)", transform=axA.transAxes,
         fontsize=5.6, va="top", color=C["tall"])
S.style_ax(axA); S.plabel(axA, "A")

# ---- B: scRNA by stage -----------------------------------------------------
axB = axd["B"]
order = ["DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP"]
blab = ["DN1", "DN2", "DN3", "ISP", "DPneg", "DPpos", "CD4SP", "CD8SP"]
scm = sc.set_index("facs_stage").reindex(order)
vals = scm["g_ZEB1"].values.astype(float)
S.add_grid(axB, "y")
axB.plot(range(len(order)), vals, "-o", color=C["tall"], ms=4, lw=1.5, zorder=3)
axB.fill_between(range(len(order)), vals, color=C["tall"], alpha=0.12, zorder=1)
axB.set_xticks(range(len(order))); axB.set_xticklabels(blab, rotation=45, fontsize=5.8, ha="right")
axB.set_ylabel("ZEB1 mean (log1p)")
axB.set_ylim(-0.02, 0.49)
axB.set_title("FACS scRNA thymus (GSE195812)", loc="left", fontsize=7.5)
axB.text(0.98, 0.96, "ρ=−0.60 (P=0.12)\npeak at DN3 / immature DP", transform=axB.transAxes,
         fontsize=5.8, va="top", ha="right")
S.style_ax(axB); S.plabel(axB, "B")

# ---- C: OLS forest ---------------------------------------------------------
axC = axd["C"]
rows = [
    ("ZEB1 ~ T-ALL", "Unadjusted"),
    ("ZEB1 ~ T-ALL + maturity", "Literature score"),
    ("ZEB1 ~ T-ALL + empirical_maturity", "Empirical score"),
]
labs, est, lo, hi, sig, source_rows = [], [], [], [], [], []
for model, lab in rows:
    r = ols[(ols.model == model) & (ols.term == "is_tall_n")].iloc[0]
    labs.append(lab); est.append(r["coef"]); lo.append(r["ci95_lo"]); hi.append(r["ci95_hi"])
    sig.append(r["p"] < 0.05)
    source_rows.append({"model": model, "maturity_definition": lab, "n_total_libraries": int(r["n"]),
                        "n_normal_stage_libraries": 6, "coef_disease": r["coef"],
                        "ci95_lo": r["ci95_lo"], "ci95_hi": r["ci95_hi"], "p": r["p"]})
pd.DataFrame(source_rows).to_csv(os.path.join(S.DATA_DIR, "SourceData_Fig2C_maturity_sensitivity.tsv"),
                                 sep="\t", index=False)
y = np.arange(len(labs))[::-1]
S.add_grid(axC, "x")
for yi, e, l, h, s in zip(y, est, lo, hi, sig):
    col = C["tall"] if s else C["grey"]
    axC.plot([l, h], [yi, yi], color=col, lw=1.4, zorder=3)
    axC.plot(e, yi, "o", color=col, ms=5.5, zorder=4)
axC.axvline(0, color="#333", lw=0.6, ls="--", zorder=2)
axC.set_yticks(y); axC.set_yticklabels(labs, fontsize=6)
axC.set_xlabel("Coefficient (95% CI)")
axC.set_title("Disease coefficient depends on score", loc="left", fontsize=7.2)
axC.set_xlim(-1.58, 0.62)
axC.text(0.98, 0.72, "41 T-ALL + 6 normal stages\nnot donor replicates",
         transform=axC.transAxes, fontsize=4.9, ha="right", va="top", color="#555")
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: replication --------------------------------------------------------
axD = axd["D"]
rep2 = rep.copy()
adult = pd.read_csv(os.path.join(S.DATA_DIR, "validation/maturity_figure/sample_scores.tsv"), sep="\t")
adult = adult[adult.is_tall.eq(True)].dropna(subset=["ZEB1", "maturity_literature"])
adult_stat = spearmanr(adult.ZEB1, adult.maturity_literature)
extra = pd.DataFrame({"cohort": ["GSE142522 adult T-ALL"], "rho": [adult_stat.statistic], "p": [adult_stat.pvalue], "n": [len(adult)]})
rep2["cohort"] = rep2["cohort"].replace({"TARGET": "TARGET (n=265)", "Pharmacotype": "Pharmacotype (n=116)"})
rep2 = pd.concat([extra, rep2], ignore_index=True)
rep2.to_csv(os.path.join(S.DATA_DIR,"SourceData_Fig2D_maturity_correlations.tsv"),sep="\t",index=False)
cols = [C["tall"] if p < 0.05 else C["grey"] for p in rep2["p"]]
x = np.arange(len(rep2))
S.add_grid(axD, "y")
axD.vlines(x, 0, rep2["rho"], color=cols, lw=.8, zorder=2)
axD.scatter(x, rep2["rho"], c=cols, s=34, zorder=3)
axD.axhline(0, color="#333", lw=0.6, zorder=2)
axD.set_xticks(x); axD.set_xticklabels(["Adult T-ALL", "TARGET", "Pharmacotype"], fontsize=6.0, rotation=0, ha="center")
axD.set_ylabel("Spearman ρ (ZEB1 vs maturity)")
axD.set_title("ZEB1-maturity association by cohort", loc="left", fontsize=7.2)
axD.set_ylim(0, 0.68)
for xi, rr, p in zip(x, rep2["rho"], rep2["p"]):
    axD.text(xi, rr + 0.015, S.stars(p), ha="center", fontsize=7)
S.style_ax(axD); S.plabel(axD, "D")

# ---- E: nearest stage ------------------------------------------------------
axE = axd["E"]
pivot = near.pivot_table(index="immunophenotype", columns="nearest_label", values="N",
                         aggfunc="sum", fill_value=0)
order_ip = [ip for ip in ["Immature", "Pre_ab_Cortical", "TCR_Mature"] if ip in pivot.index]
pivot = pivot.reindex(order_ip)
stage_cols = [c for c in ["CD34+", "ISP", "SP8"] if c in pivot.columns]
palette = {"CD34+": C["ball"], "ISP": C["gold"], "SP8": C["tall"]}
bottom = np.zeros(len(pivot)); x = np.arange(len(pivot))
for sc_ in stage_cols:
    axE.bar(x, pivot[sc_].values, bottom=bottom, label=sc_, color=palette[sc_],
            width=0.6, edgecolor="white", linewidth=0.4, zorder=3)
    bottom += pivot[sc_].values
axE.set_xticks(x); axE.set_xticklabels([ip.replace("_", "\n") for ip in pivot.index], fontsize=5.8)
axE.set_ylabel("Adult T-ALL (n)")
axE.set_title("Nearest normal stage", loc="left", fontsize=7.5)
axE.legend(title="Nearest", ncol=1, fontsize=6, title_fontsize=6, loc="upper right")
S.style_ax(axE); S.plabel(axE, "E")

# ---- F: scRNA stage coupling ----------------------------------------------
axF = axd["F"]
sel = ["CD3E", "CD3D", "TRAC", "CD8B", "LEF1", "IL7R", "ZEB1", "ZEB2", "LMO2", "MEF2C", "CD34", "LYL1"]
srr = scr.set_index("gene").reindex(sel).dropna(subset=["rho_vs_stage"]).sort_values("rho_vs_stage")
# Existing stage means; standardization affects display only.
means = scm[["g_" + g for g in srr.index]].T
means.index = srr.index
sd = means.std(axis=1, ddof=0).replace(0, np.nan)
z = means.sub(means.mean(axis=1), axis=0).div(sd, axis=0)
im = axF.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                interpolation="nearest")
axF.set_yticks(range(len(srr))); axF.set_yticklabels(srr.index, fontsize=5.8)
axF.set_xticks(range(len(order))); axF.set_xticklabels(blab, rotation=55,
                                                       ha="right", fontsize=5.1)
for label in axF.get_yticklabels():
    if label.get_text() == "ZEB1":
        label.set_color(C["tall"]); label.set_fontweight("bold")
for yi, rho in enumerate(srr.rho_vs_stage):
    axF.text(1.04, yi, f"{rho:+.2f}", transform=axF.get_yaxis_transform(),
             va="center", fontsize=5.2, clip_on=False)
axF.text(1.04, 1.035, "ρ", transform=axF.transAxes, fontsize=6)
axF.set_title("Stage-resolved gene expression", loc="left", fontsize=7.2)
axF.tick_params(length=0)
for spine in axF.spines.values():
    spine.set_visible(False)
cbax = axF.inset_axes([0.16, -0.32, 0.65, 0.045])
cb = fig.colorbar(im, cax=cbax, orientation="horizontal", ticks=[-2, 0, 2])
cb.ax.tick_params(labelsize=5, length=2, pad=1)
cb.set_label("Within-gene z score", fontsize=5.2, labelpad=1)
cb.outline.set_visible(False)
S.plabel(axF, "F")
source = means.rename_axis("gene").reset_index().melt(
    id_vars="gene", var_name="facs_stage", value_name="mean_log1p")
zlong = z.rename_axis("gene").reset_index().melt(
    id_vars="gene", var_name="facs_stage", value_name="display_zscore")
source = source.merge(zlong, on=["gene", "facs_stage"], validate="one_to_one")
source = source.merge(srr.reset_index(), on="gene", validate="many_to_one")
source.to_csv(os.path.join(S.DATA_DIR, "SourceData_Fig2F_stage_expression.tsv"),
              sep="\t", index=False)

# ---- G: Park/HTA independent thymus atlas (donor-level) --------------------
axG = axd["G"]
park = pd.read_csv(os.path.join(S.DATA_DIR, "F4_park_donor_celltype.tsv"), sep="\t")
keep = {
    "double negative thymocyte": "DN",
    "double-positive, alpha-beta thymocyte": "DP",
    "CD4-positive, alpha-beta T cell": "CD4 SP",
    "CD8-positive, alpha-beta T cell": "CD8 SP",
}
pg = park[park["celltype"].isin(keep) & (park["n_cells"] >= 20)].copy()
pg["stage"] = pg["celltype"].map(keep)
order_g = ["DN", "DP", "CD4 SP", "CD8 SP"]
S.add_grid(axG, "y")
rng = np.random.default_rng(1)
for i, st in enumerate(order_g):
    sub = pg[pg["stage"] == st]
    xj = i + rng.uniform(-0.12, 0.12, len(sub))
    axG.scatter(xj, sub["ZEB1"].astype(float), s=14, color=C["tall"], alpha=0.55,
                edgecolor="#333", linewidth=0.25, zorder=3)
    q = sub["ZEB1"].astype(float)
    axG.plot([i - 0.18, i + 0.18], [q.median(), q.median()], color="#222", lw=1.6, zorder=4)
axG.set_xticks(range(4)); axG.set_xticklabels(order_g)
axG.set_ylabel("Donor-mean ZEB1 (log1p)")
axG.set_ylim(0.05, 2.85)
axG.set_title("ZEB1 across thymocyte stages | independent Park/HTA atlas",
              loc="left", fontsize=7.5)
axG.text(0.01, 0.92, "One point per donor-stage group; at least 20 cells", transform=axG.transAxes,
         fontsize=5.8, color=C["tall"])
S.style_ax(axG); S.plabel(axG, "G")

S.save(fig, "Figure2_development")
