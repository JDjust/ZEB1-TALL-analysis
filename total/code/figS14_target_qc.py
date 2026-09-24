"""Supplementary Figure 14: TARGET STAR count/QC sensitivity of Figure 4."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import _style as S

S.set_style()
C = S.C
V = os.path.join(S.DATA_DIR, "validation")
sample = pd.read_csv(os.path.join(V, "target_global_axis_sample_qc.tsv"), sep="\t")
partial = pd.read_csv(os.path.join(V, "target_global_axis_partial.tsv"), sep="\t")
summary = pd.read_csv(os.path.join(V, "qc_sensitivity_summary.tsv"), sep="\t")
key = pd.read_csv(os.path.join(V, "qc_sensitivity_keygenes.tsv"), sep="\t")
assert len(sample) == 265 and set(summary.model) == {
    "original", "source", "assigned", "nofeature", "mito"}

fig = S.new_fig(183, 204, constrained=False)
gs = fig.add_gridspec(3, 2, left=.14, right=.97, top=.95, bottom=.08,
                      hspace=.58, wspace=.55)

ax = fig.add_subplot(gs[0, 0])
metrics = ["expression_PC1", "median_log2tpm", "N_noFeature",
           "assigned_fraction_all_reported", "mitochondrial_count_fraction",
           "detected_genes_count_gt0"]
labels = ["Expression PC1", "Median log2TPM", "STAR no-feature counts",
          "STAR assigned fraction", "Mitochondrial fraction", "Detected genes"]
p = partial.set_index("metric").loc[metrics]
y = np.arange(len(p))
cols = [C["tall"] if v > 0 else C["ball"] for v in p.partial_rho_rank_residual]
ax.hlines(y, 0, p.partial_rho_rank_residual, color=cols, lw=.8)
ax.scatter(p.partial_rho_rank_residual, y, c=cols, s=24, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.7)
ax.invert_yaxis(); ax.axvline(0, color="#333", lw=.6)
ax.set_xlim(-.7, .9)
ax.set_xlabel("Partial rank ρ with ZEB1")
ax.set_title("Adjusted expression/QC associations", loc="left", fontsize=7)
S.style_ax(ax); S.panel_label(ax, "A", x=-.18)

ax = fig.add_subplot(gs[0, 1])
palette = {"TAL": C["ball"], "TLX": C["tall"], "HOXA": C["aml"],
           "LMO2/LYL1": C["purple"], "NKX2_1": C["teal"], "Unknown": C["grey"]}
for subtype, group in sample.groupby("subtype"):
    ax.scatter(np.log10(group.N_noFeature), group.ZEB1_log, s=12,
               color=palette.get(subtype, C["grey"]), alpha=.55, label=subtype,
               edgecolors="none")
ax.set_xlabel("log10 STAR no-feature count")
ax.set_ylabel("ZEB1 log2(TPM+1)")
ax.set_title("N_noFeature association is descriptive", loc="left", fontsize=7)
ax.legend(loc="lower right", fontsize=4.7, ncol=2)
S.style_ax(ax); S.panel_label(ax, "B", x=-.17)

ax = fig.add_subplot(gs[1, 0])
k = key[key.model.isin(["original", "nofeature"])].copy()
order = ["GATA3", "TCF7", "BCL11B", "IL7R", "ZEB2", "LMO2", "LYL1", "CD34"]
for i, gene in enumerate(order):
    for model, offset, col in [("original", -.12, C["grey"]),
                               ("nofeature", .12, C["tall"])]:
        row = k[(k.gene == gene) & (k.model == model)].iloc[0]
        # Descriptive Wald-scale interval from limma coefficient and moderated t.
        se = abs(row.logFC / row.t) if row.t != 0 else np.nan
        ax.plot([row.logFC - 1.96*se, row.logFC + 1.96*se],
                [i+offset, i+offset], color=col, lw=1)
        ax.scatter(row.logFC, i+offset, s=18, color=col, zorder=3)
ax.axvline(0, color="#333", lw=.6)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=6)
ax.invert_yaxis()
ax.set_xlabel("ZEB1 coefficient (log2TPM per SD)")
ax.set_title("Key gene sensitivity to no-feature adjustment", loc="left", fontsize=7)
ax.scatter([], [], color=C["grey"], label="original")
ax.scatter([], [], color=C["tall"], label="+ no-feature")
ax.legend(loc="lower right", fontsize=5.3)
S.style_ax(ax); S.panel_label(ax, "C", x=-.18)

ax = fig.add_subplot(gs[1, 1])
order_models = ["original", "source", "assigned", "nofeature", "mito"]
s = summary.set_index("model").loc[order_models]
labels = ["Subtype + age", "+ blood/marrow", "+ assigned fraction",
          "+ no-feature", "+ mitochondrial fraction"]
y = np.arange(len(s))
ax.scatter(s.fdr05, y, c=[C["tall"] if m == "nofeature" else C["grey"] for m in order_models],
           s=28, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.6)
ax.invert_yaxis(); ax.set_xlim(0, 13090)
for yi, n in zip(y, s.fdr05):
    ax.text(n+120, yi, f"{n:,}", va="center", fontsize=5.3)
ax.set_xlabel("Genes with BH FDR<0.05 (of 13,090)")
ax.set_title("Broad association attenuates, but remains broad", loc="left", fontsize=7)
S.style_ax(ax); S.panel_label(ax, "D", x=-.18)

matched = pd.read_csv(os.path.join(V, "matched_gene_sets", "random_sets.tsv.gz"), sep="\t")
competitive = pd.read_csv(os.path.join(V, "matched_gene_sets", "competitive_summary.tsv"), sep="\t")
for j, cohort in enumerate(["TARGET", "Pharmacotype"]):
    ax = fig.add_subplot(gs[2, j])
    for scheme, color, label in [("mean_detection", C["grey"], "Mean + detection"),
                                 ("mean_detection_sd", C["ball"], "+ expression SD")]:
        values = np.sort(matched.loc[matched.cohort.eq(cohort) & matched.matching.eq(scheme), "rho"])
        ax.plot(values, np.arange(1, len(values)+1)/len(values), color=color, lw=1.1, label=label)
    stats = competitive[competitive.cohort.eq(cohort)]
    observed = stats.observed_rho.iloc[0]
    ax.axvline(observed, color=C["tall"], ls="--", lw=1.1, label="Observed 3-gene score")
    ax.text(.03, .97, f"Observed ρ={observed:.2f}\nCompetitive P={stats.competitive_upper_tail_p.min():.2f}–"
            f"{stats.competitive_upper_tail_p.max():.2f}", transform=ax.transAxes, va="top", fontsize=5.5)
    ax.set_xlim(-1, 1); ax.set_ylim(0, 1.03)
    ax.set_xlabel("Gene-set score versus ZEB1: Spearman ρ")
    ax.set_ylabel("Cumulative fraction of matched sets")
    ax.set_title(f"{cohort}: 10,000 sets per matching scheme", loc="left", fontsize=6.5)
    ax.legend(loc="lower right", frameon=False, fontsize=4.8)
    S.style_ax(ax); S.panel_label(ax, "E" if j == 0 else "F", x=-.18)

fig.text(.14, .025, "Exploratory specificity and QC checks; these comparisons do not identify a technical or causal mechanism.",
         fontsize=6, color="#555")
S.save_supp(fig, "FigureS14_TARGET_global_axis_QC")
