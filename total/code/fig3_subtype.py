"""
Figure 3 | ZEB1 varies with molecular subtype/maturation within T-ALL:
TLX-high is directionally consistent across cohorts; study-defined immature
classes are displayed separately, and the evaluated lesion screen is bounded.
"""
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import _style as S

S.set_style()
C = S.C

sub = pd.read_csv(S.tbl("module3", "M3_3.2_3.4_subtype_effects.tsv"), sep="\t")
counts = pd.read_csv(S.tbl("module3", "M3_3.1_subtype_counts.tsv"), sep="\t")
meta = pd.read_csv(S.tbl("module3", "M3_r2_3.5_meta.tsv"), sep="\t")
emat = pd.read_csv(S.tbl("module3", "M3_r2_3.5_effect_matrix.tsv"), sep="\t")
etp = pd.read_csv(S.tbl("module3", "M3_pharmacotype_ETP_age_model.tsv"), sep="\t")

fig = S.new_fig(183, 164, constrained=False)
axd = {key: fig.add_axes(rect) for key, rect in {
    "A": [.13, .745, .20, .205], "B": [.47, .745, .20, .205],
    "C": [.79, .745, .18, .205], "D": [.13, .415, .20, .205],
    "E": [.47, .415, .20, .205], "F": [.79, .415, .18, .205],
    "G": [.13, .085, .84, .215]}.items()}

colmap = {"TAL": C["ball"], "TLX": C["tall"], "HOXA": C["aml"],
          "LMO2/LYL1": C["purple"], "NKX2_1": C["teal"], "Unknown": C["grey"]}
display_label = lambda value: "LMO-associated" if value == "LMO2/LYL1" else value

# ---- A: composition --------------------------------------------------------
axA = axd["A"]
co = counts.sort_values("N", ascending=True)
cols = [colmap.get(s, C["grey"]) for s in co["subtype"]]
S.add_grid(axA, "x")
axA.barh(co["subtype"], co["N"], color=cols, height=0.7, edgecolor="white", lw=0.35, zorder=3)
for s, n in zip(co["subtype"], co["N"]):
    axA.text(n + 1.5, s, str(int(n)), va="center", fontsize=6.0)
axA.set_xlabel("Patients")
axA.set_yticks(range(len(co)), [display_label(v) for v in co.subtype])
axA.set_xlim(0, max(co["N"]) * 1.22)
axA.set_title("TARGET groups (n=265)", loc="left", fontsize=7.5)
S.style_ax(axA); S.plabel(axA, "A")

# ---- B: ZEB1 by subtype ----------------------------------------------------
axB = axd["B"]
b = sub[sub.gene == "ZEB1"].copy()
b["st"] = b["contrast"].str.replace(" vs rest", "", regex=False)
b = b.sort_values("hedges_g")
y = np.arange(len(b))
cols = [C["up"] if (g > 0 and p < 0.05) else (C["down"] if (g < 0 and p < 0.05) else C["ns"])
        for g, p in zip(b["hedges_g"], b["fdr"])]
S.add_grid(axB, "x")
axB.scatter(b["hedges_g"], y, c=cols, s=26, edgecolor="white", linewidth=.4, zorder=5)
axB.errorbar(b["hedges_g"], y, xerr=[b["hedges_g"] - b["ci95_lo"], b["ci95_hi"] - b["hedges_g"]],
             fmt="none", ecolor="#333", lw=0.8, capsize=1.8, zorder=4)
axB.axvline(0, color="#333", lw=0.6, zorder=2)
axB.set_yticks(y); axB.set_yticklabels([display_label(v) for v in b["st"]], fontsize=6.2)
axB.set_xlabel("ZEB1 Hedges g (vs rest)")
axB.set_title("ZEB1 by subtype (TARGET)", loc="left", fontsize=7.5)
for yi, g, p, hi in zip(y, b["hedges_g"], b["fdr"], b["ci95_hi"]):
    axB.text(float(hi) + 0.05, yi, S.stars(p), va="center", ha="left", fontsize=6)
axB.set_xlim(-1.55, 1.55)
S.style_ax(axB); S.plabel(axB, "B")

# ---- C: ZEB1 vs LMO2 axis --------------------------------------------------
axC = axd["C"]
z = sub[sub.gene == "ZEB1"].set_index("contrast")["hedges_g"]
l = sub[sub.gene == "LMO2"].set_index("contrast")["hedges_g"]
axC.axhline(0, color="#B9BEC4", lw=0.6, ls=":"); axC.axvline(0, color="#B9BEC4", lw=0.6, ls=":")
nudge = {"NKX2_1": (10, 11), "TAL": (-16, -12), "HOXA": (10, -11),
         "TLX": (0, 8), "LMO2/LYL1": (-5, -12)}
for contrast in z.index:
    st = contrast.replace(" vs rest", "")
    x, yv = l[contrast], z[contrast]
    axC.scatter(x, yv, s=54, color=colmap.get(st, C["grey"]), edgecolor="#333", linewidth=0.5, zorder=4)
    dx, dy = nudge.get(st, (0, 6))
    axC.annotate(display_label(st), (x, yv), fontsize=5.8, xytext=(dx, dy), textcoords="offset points",
                 ha="center",
                 arrowprops=dict(arrowstyle="-", color="#B8BCC2", lw=0.4) if st in ("NKX2_1", "TAL") else None)
axC.set_xlabel("LMO2 Hedges g (vs rest)")
axC.set_ylabel("ZEB1 Hedges g (vs rest)")
axC.set_title("Joint subtype effects", loc="left", fontsize=7.5)
axC.set_xlim(-1.15, 2.05); axC.set_ylim(-1.15, 1.05)
S.style_ax(axC); S.plabel(axC, "C")

# ---- D: cross-cohort heatmap ----------------------------------------------
axD = axd["D"]
order_sub = ["TLX", "TAL", "HOXA", "Immature", "NKX"]
cohorts = ["TARGET-ALL-P2", "GSE62156", "GSE26713", "GSE110636"]
M = emat.set_index("subtype").reindex(order_sub)[cohorts]
im = axD.imshow(M.values.astype(float), cmap="RdBu_r", vmin=-1.3, vmax=1.3, aspect="auto")
axD.set_xticks(range(len(cohorts)))
axD.set_xticklabels(["TARGET", "GSE62156", "GSE26713", "GSE110636"], rotation=35, ha="right", fontsize=5.8)
axD.set_yticks(range(len(order_sub))); axD.set_yticklabels(["TLX", "TAL", "HOXA", "Study-defined\nimmature*", "NKX"], fontsize=6.2)
for i in range(len(order_sub)):
    for j in range(len(cohorts)):
        v = M.values[i, j]
        if not np.isnan(v):
            axD.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.4,
                     color="white" if abs(v) > 0.75 else "#222")
cb = fig.colorbar(im, ax=axD, fraction=0.046, pad=0.03)
cb.set_label("ZEB1 g (vs rest)", fontsize=6); cb.ax.tick_params(labelsize=5.5)
axD.set_title("Cross-cohort effects*", loc="left", fontsize=7.1)
S.plabel(axD, "D")

# ---- E: meta forest --------------------------------------------------------
axE = axd["E"]
mrows = {"TLX vs rest": ("TLX", C["tall"])}
y0 = 0; yt, ytl = [], []
S.add_grid(axE, "x")
for contrast, (lab, col) in mrows.items():
    mrow = meta[meta.contrast == contrast].iloc[0]
    key = "TLX"
    pts = emat.set_index("subtype").loc[key, cohorts].astype(float).dropna()
    axE.scatter(pts.values, [y0] * len(pts), s=18, color=col, alpha=0.45, zorder=3)
    axE.plot([mrow["ci95_lo"], mrow["ci95_hi"]], [y0, y0], color=col, lw=1.7, zorder=4)
    axE.plot(mrow["hedges_g"], y0, "D", color=col, ms=7.5, zorder=5)
    axE.text(mrow["hedges_g"], y0 + 0.22, f"g={mrow['hedges_g']:.2f}, I²={mrow['I2']:.0f}%",
             ha="center", fontsize=5.6)
    yt.append(y0); ytl.append(f"{lab}\n(k=4)")
    y0 -= 1
axE.axvline(0, color="#333", lw=0.6, ls="--", zorder=2)
axE.set_yticks(yt); axE.set_yticklabels(ytl, fontsize=6.2)
sens = pd.read_csv(S.DATA_DIR + "/validation/tlx_meta_sensitivity.tsv", sep="\t")
hk = sens[sens.method.eq("REML_modified_HK")].iloc[0]
axE.plot([hk.ci95_lo, hk.ci95_hi], [-.65, -.65], color=C["ball"], lw=1.7)
axE.scatter([hk.estimate], [-.65], marker="D", color=C["ball"], s=40, zorder=5)
axE.set_yticks([0, -.65], ["DL normal\n(k=4)", "REML\nmodified HK"], fontsize=5.7)
axE.set_ylim(-1.1, 0.75)
axE.set_xlabel("Random-effects Hedges g")
axE.set_title("TLX meta-analysis (k=4)", loc="left", fontsize=7.5)
S.style_ax(axE); S.plabel(axE, "E")

# ---- F: ETP ----------------------------------------------------------------
axF = axd["F"]
r = etp[etp.term == "etpnotETP"].iloc[0]
S.add_grid(axF, "y")
axF.scatter([0], [r["estimate"]], color=C["up"], s=42, edgecolor="white", linewidth=.4, zorder=5)
axF.errorbar([0], [r["estimate"]], yerr=[[r["estimate"] - r["ci95_lo"]], [r["ci95_hi"] - r["estimate"]]],
             fmt="none", ecolor="#333", lw=0.9, capsize=3, zorder=4)
axF.axhline(0, color="#333", lw=0.6, zorder=2)
axF.set_xticks([0]); axF.set_xticklabels(["non-ETP vs ETP"], fontsize=6.5)
axF.set_ylabel("ZEB1 Δ log2 (age-adj.)")
axF.set_title("ETP comparison", loc="left", fontsize=7.5)
axF.text(0, r["estimate"] + 0.06, S.stars(r["p"]), ha="center", fontsize=8)
axF.set_xlim(-0.9, 0.9)
S.style_ax(axF); S.plabel(axF, "F")

# G: original fine labels, reusing the user's completed summaries and patients.
base = Path(S.DATA_DIR) / "validation/lmo2_subgroups_existing"
patient = pd.read_csv(base / "TALL_clinical_annotation_full.csv")
summary = pd.read_csv(base / "TALL_subtype_ZEB1_stats.csv").set_index("Subtype")
summary.reset_index().to_csv(Path(S.DATA_DIR) / "SourceData_Fig3G_subtype_summary.tsv", sep="\t", index=False)
levels = ["LMO1/2", "LMO2_LYL1", "TAL1", "TAL2", "TLX1", "TLX3", "HOXA", "NKX2_1", "Unknown"]
colors = [C["aml"], C["purple"], C["ball"], "#56B4E9", C["tall"], "#CC79A7", C["green"], C["teal"], C["grey"]]
axG = axd["G"]
rng = np.random.default_rng(20260924)
S.add_grid(axG, "y")
for i, (label, color) in enumerate(zip(levels, colors)):
    values = patient.loc[patient.MolecularSubtype.eq(label), "ZEB1"]
    row = summary.loc[label]
    assert len(values) == int(row.n)
    axG.scatter(i + rng.uniform(-.2, .2, len(values)), values, s=7,
                color=color, alpha=.6, linewidths=0, zorder=3)
    axG.plot([i, i], [row.Q1, row.Q3], color="#333", lw=2.4, zorder=4)
    axG.plot([i-.15, i+.15], [row.Median_ZEB1]*2, color="#111", lw=1.3, zorder=5)
axG.set_xticks(range(len(levels)), [f"{label.replace('_', '/')}\nn={int(summary.loc[label,'n'])}" for label in levels], fontsize=5.8)
axG.set_ylabel("ZEB1 log2(TPM + 1)")
axG.set_title("Original TARGET subtypes | patient values, median and IQR", loc="left", fontsize=7.5)
S.style_ax(axG); S.plabel(axG, "G")
patient[["Sample_ID", "MolecularSubtype", "ZEB1"]].to_csv(
    Path(S.DATA_DIR) / "SourceData_Fig3G_original_subtypes.tsv", sep="\t", index=False)
S.save(fig, "Figure3_subtype")
