"""Supplementary Figure S15: fixed three-gene patient pseudobulk check."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import _style as S

S.set_style()
C = S.C
V = os.path.join(S.DATA_DIR, "validation", "three_gene_patient")
patients = pd.read_csv(os.path.join(V, "patient_three_gene_score.tsv"), sep="\t")
assoc = pd.read_csv(os.path.join(V, "three_gene_associations.tsv"), sep="\t")
assert patients.groupby("cohort").size().to_dict() == {"GSE227122": 10, "GSE248287": 15}

fig = S.new_fig(183, 177, constrained=False)
gs = fig.add_gridspec(3, 2, left=.12, right=.96, top=.94, bottom=.08,
                      height_ratios=[1.15, .72, 1], hspace=.72, wspace=.38)
for j, cohort in enumerate(["GSE227122", "GSE248287"]):
    ax = fig.add_subplot(gs[0, j])
    d = patients[patients.cohort.eq(cohort)]
    r = assoc[(assoc.cohort.eq(cohort)) & (assoc.feature.eq("three_gene_score"))].iloc[0]
    color = C["ball"] if cohort == "GSE227122" else C["tall"]
    ax.scatter(d.ZEB1, d.three_gene_score, s=34, color=color,
               edgecolor="#333", linewidth=.4, alpha=.82)
    ax.axhline(0, color="#B8BCC2", lw=.55, ls=":")
    ax.set_xlabel("Patient ZEB1 log2(1+CPM)")
    ax.set_ylabel("Three-gene mean z score")
    ax.set_title(f"{cohort} (n={len(d)} patients)", loc="left", fontsize=7.5)
    ax.text(.03 if j == 0 else .97, .97 if j == 0 else .04,
            f"ρ={r.rho:+.2f}; permutation P={r.permutation_p:.2f}\n"
            f"leave-one-out ρ: {r.loo_rho_min:+.2f} to {r.loo_rho_max:+.2f}",
            transform=ax.transAxes, va="top" if j == 0 else "bottom",
            ha="left" if j == 0 else "right", fontsize=5.6)
    S.style_ax(ax); S.panel_label(ax, "A" if j == 0 else "B", x=-.16)

ax = fig.add_subplot(gs[1, :])
order = ["three_gene_score", "GATA3", "TCF7", "BCL11B"]
cohorts = ["GSE227122", "GSE248287"]
matrix = assoc.pivot(index="feature", columns="cohort", values="rho").loc[order, cohorts]
im = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-.75, vmax=.75, aspect="auto")
ax.set_xticks(range(2), ["Candidate malignant (GSE227122)",
                         "Author malignant (GSE248287)"], fontsize=6)
ax.set_yticks(range(4), ["3-gene score", "GATA3", "TCF7", "BCL11B"], fontsize=6)
ax.tick_params(length=0)
for yi in range(4):
    for xi in range(2):
        ax.text(xi, yi, f"{matrix.iloc[yi, xi]:+.2f}", ha="center", va="center", fontsize=6,
                color="white" if abs(matrix.iloc[yi, xi]) > .4 else C["ink"])
cb = fig.colorbar(im, ax=ax, fraction=.024, pad=.02)
cb.set_label("Spearman ρ", fontsize=5.8); cb.ax.tick_params(labelsize=5)
ax.set_title("Direction differs by cohort; no test passes within-cohort BH FDR<0.05",
             loc="left", fontsize=7)
S.panel_label(ax, "C", x=-.07)
qc_dir = os.path.join(S.DATA_DIR, "validation", "scrna_patient_qc")
qc = pd.read_csv(os.path.join(qc_dir, "patient_qc_scores.tsv"), sep="\t")
thresholds = pd.read_csv(os.path.join(qc_dir, "cell_threshold_sensitivity.tsv"), sep="\t")
ax = fig.add_subplot(gs[2, 0])
for cohort, color in zip(cohorts, [C["ball"], C["tall"]]):
    d = qc[qc.cohort.eq(cohort)]
    ax.scatter(d.median_umi, d.zeb1_detection_fraction * 100, s=22, color=color,
               alpha=.8, edgecolor="white", lw=.3, label=cohort)
ax.set_xscale("log")
ax.set_xlabel("Median UMI per malignant cell (log axis)")
ax.set_ylabel("Cells detecting ZEB1 (%)")
ax.set_title("Patient-level depth and detection", loc="left", fontsize=7)
ax.legend(frameon=False, fontsize=5.5, loc="best")
S.style_ax(ax); S.panel_label(ax, "D", x=-.16)
ax = fig.add_subplot(gs[2, 1])
for cohort, color in zip(cohorts, [C["ball"], C["tall"]]):
    d = thresholds[thresholds.cohort.eq(cohort)].sort_values("min_cells")
    ax.plot(range(len(d)), d.rho, "o-", color=color, markersize=3, lw=.9, label=cohort)
    for i, (_, row) in enumerate(d.iterrows()):
        ax.annotate(f"n={row.n_patients}", (i, row.rho), xytext=(0, 7),
                    textcoords="offset points", ha="center", fontsize=5)
ax.axhline(0, color="#AAA", lw=.6, ls=":")
ax.set_xticks(range(5), [30, 100, 250, 500, 1000], fontsize=6)
ax.set_ylim(-.5, .7)
ax.set_xlabel("Minimum malignant cells per patient")
ax.set_ylabel("Score versus ZEB1: Spearman ρ")
ax.set_title("Cell-count threshold sensitivity", loc="left", fontsize=7)
S.style_ax(ax); S.panel_label(ax, "E", x=-.16)
S.save_supp(fig, "FigureS15_three_gene_patient")
