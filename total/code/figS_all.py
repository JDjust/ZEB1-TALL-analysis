# -*- coding: utf-8 -*-
"""
Supplementary Figures S1–S12 for the ZEB1 / T-ALL reappraisal (S13 is generated separately).

Output: total/figures/Supplementary/*.pdf  +  png/*.png
Style matches the seven main figures (_style.py).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import _style as S

S.set_style()
C = S.C
MM = S.MM
D = S.DATA_DIR


def _p(path_mod, name):
    return S.tbl(path_mod, name)


def _sd(name):
    return __import__("os").path.join(D, name)


def _short(s, n=28):
    s = str(s).replace("HALLMARK_", "").replace("REACTOME_", "").replace("GOBP_", "")
    s = s.replace("_", " ")
    return s if len(s) <= n else s[: n - 1] + "…"


def _cohort_label(x):
    return {"TARGET-ALL-P2": "TARGET", "StJude_Pharmacotype": "Pharmacotype",
            "GSE13159": "MILE"}.get(x, x)


# ===========================================================================
# S1  QC and sensitivity (M1)
# ===========================================================================
def fig_s1():
    counts = pd.read_csv(_p("module1", "M1_1.1_sample_counts.tsv"), sep="\t")
    pair = pd.read_csv(_p("module1", "M1_1.2_1.5_pairwise_vs_TALL.tsv"), sep="\t")
    sens = pd.read_csv(_p("module1", "M1_sensitivity_ZEB1_probe_212764_at.tsv"), sep="\t")
    loo = pd.read_csv(_p("module1", "M1_1.16_leave_one_cohort_out.tsv"), sep="\t")
    meta = pd.read_csv(_p("module1", "M1_1.14_random_effects_meta.tsv"), sep="\t")
    plat = pd.read_csv(_p("module1", "M1_1.18_platform_effects.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 135 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.42)

    ax = fig.add_subplot(gs[0, 0])
    d = counts.groupby("disease", as_index=False)["N"].sum().sort_values("N")
    cols = [C["tall"] if x == "T-ALL" else C["grey"] for x in d["disease"]]
    ax.barh(d["disease"], d["N"], color=cols, height=0.7, edgecolor="white", lw=0.3)
    ax.set_xlabel("Samples (MILE, GSE13159)")
    ax.set_title("Cohort composition", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.22)

    ax = fig.add_subplot(gs[0, 1])
    # mean-probe vs canonical 212764_at
    meanp = pair[pair["contrast"].str.contains("B-ALL|AML|Normal")].copy()
    meanp = meanp.drop_duplicates("contrast")
    sp = sens.copy()
    sp["cmp"] = sp["contrast"].str.replace("T-ALL vs ", "")
    meanp["cmp"] = meanp["contrast"].str.replace("T-ALL vs ", "")
    m = meanp.merge(sp, on="cmp", suffixes=("_mean", "_probe"))
    y = np.arange(len(m))
    ax.hlines(y, m["hedges_g_mean"], m["hedges_g_probe"], color="#B8BCC2", lw=1.0)
    ax.scatter(m["hedges_g_mean"], y, s=22, color=C["grey"], zorder=3, label="4-probe mean")
    ax.scatter(m["hedges_g_probe"], y, s=28, color=C["tall"], zorder=3, label="212764_at only")
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(m["cmp"], fontsize=6.2)
    ax.set_xlabel("Hedges g (T-ALL vs comparator)")
    ax.set_title("Probe sensitivity", loc="left")
    ax.legend(loc="lower right", fontsize=5.8)
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.32)

    ax = fig.add_subplot(gs[1, 0])
    g0 = float(meta["hedges_g"].iloc[0])
    lo0, hi0 = float(meta["ci95_lo"].iloc[0]), float(meta["ci95_hi"].iloc[0])
    rows = [("All three (RE)", g0, lo0, hi0)]
    for _, r in loo.iterrows():
        rows.append((f"leave out {_cohort_label(r['left_out'])}", r["hedges_g"],
                     r["ci95_lo"], r["ci95_hi"]))
    y = np.arange(len(rows))[::-1]
    for yi, (_, g, lo, hi) in zip(y, rows):
        ax.plot([lo, hi], [yi, yi], color=C["ink"], lw=1.1)
        ax.scatter([g], [yi], s=28, color=C["tall"] if yi == y[0] else C["ball"], zorder=3)
    ax.axvline(0, color="#999", lw=0.6, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.2)
    ax.set_xlabel("Random-effects g (T-ALL vs B-ALL)")
    ax.set_title("Leave-one-cohort-out", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "C", x=-0.42)

    ax = fig.add_subplot(gs[1, 1])
    p = plat.drop_duplicates("cohort").copy()
    y = np.arange(len(p))
    ax.hlines(y, p["ci95_lo"], p["ci95_hi"], color="#B8BCC2", lw=1.1)
    cols = [C["tall"] if "microarray" in str(t).lower() else C["ball"] for t in p["platform"]]
    ax.scatter(p["hedges_g"], y, s=28, color=cols, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{_cohort_label(c)}\n{t}" for c, t in zip(p["cohort"], p["platform"])],
                       fontsize=5.8)
    ax.set_xlabel("Hedges g (T-ALL vs B-ALL)")
    ax.set_title("Platform replication", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "D", x=-0.38)

    S.save_supp(fig, "FigureS1_QC_sensitivity")


# ===========================================================================
# S2  Park/HTA thymus + GSE206710  (§3.2)
# ===========================================================================
def fig_s2():
    park = pd.read_csv(_sd("SourceData_ParkHTA_donor_celltype.tsv"), sep="\t")
    st = pd.read_csv(_sd("SourceData_FigS2_stage_means.tsv"), sep="\t")
    cd3 = pd.read_csv(_sd("SourceData_FigS2_CD3_paired.tsv"), sep="\t")
    sn = pd.read_csv(_p("module2", "M2_r2_GSE206710_stage_n.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 195 * MM))
    gs = GridSpec(3, 2, figure=fig, hspace=0.65, wspace=0.42)

    cmap = {
        "double negative thymocyte": "DN",
        "double-positive, alpha-beta thymocyte": "DP",
        "CD4-positive, alpha-beta T cell": "CD4SP",
        "CD8-positive, alpha-beta T cell": "CD8SP",
    }
    order = ["DN", "DP", "CD4SP", "CD8SP"]
    p = park[park["celltype"].isin(cmap)].copy()
    p = p[p["n_cells"] >= 20]
    p["stage"] = p["celltype"].map(cmap)
    p["stage"] = pd.Categorical(p["stage"], categories=order, ordered=True)
    p = p.dropna(subset=["stage", "ZEB1"])

    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(1)
    for i, stg in enumerate(order):
        sub = p[p.stage == stg]
        x = i + rng.uniform(-0.12, 0.12, len(sub))
        ax.scatter(x, sub["ZEB1"], s=14, color=C["tall"], edgecolor="#333", lw=0.25, zorder=3)
        ax.hlines(sub["ZEB1"].median(), i - 0.22, i + 0.22, color=C["ink"], lw=1.1, zorder=4)
    ax.set_xticks(range(4))
    n_don = [p[p.stage == s]["donor"].nunique() for s in order]
    ax.set_xticklabels([f"{s}\nn={n} donors" for s, n in zip(order, n_don)], fontsize=5.8)
    ax.set_ylabel("Donor-mean ZEB1")
    ax.set_title("Park/HTA T cells (donor n, not cells)", loc="left")
    ax.text(0.03, 0.05, "DN/DP high, SP lower",
            transform=ax.transAxes, fontsize=5.6, color="#555")
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    dn = p[p.stage == "DN"].set_index("donor")["ZEB1"]
    cd4 = p[p.stage == "CD4SP"].set_index("donor")["ZEB1"]
    both = dn.index.intersection(cd4.index)
    x = np.array([0, 1])
    for d in both:
        ax.plot(x, [dn[d], cd4[d]], "-o", color=C["ink"], lw=0.8, ms=3.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["DN", "CD4SP"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("ZEB1 (donor mean)")
    ax.set_title(f"Paired DN → CD4SP (n = {len(both)} donors)", loc="left")
    if len(both) >= 3:
        w = wilcoxon(dn.loc[both], cd4.loc[both])
        ax.text(0.97, 0.92, f"Wilcoxon P = {w.pvalue:.1e}\nDN higher",
                transform=ax.transAxes, ha="right", va="top", fontsize=5.8)
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.18)

    ax = fig.add_subplot(gs[1, 0])
    order206 = ["DN_early", "DP_immature", "DP_mature"]
    st = st.set_index("inferred_stage").reindex(order206)
    nmap = sn.set_index("inferred_stage")["n_cells"]
    xx = np.arange(len(order206))
    ax.plot(xx, st["g_ZEB1"], "-o", color=C["tall"], lw=1.4, ms=5, label="ZEB1")
    ax.plot(xx, st["g_LMO2"], "-s", color=C["purple"], lw=1.0, ms=4, label="LMO2")
    ax.plot(xx, st["g_ZEB2"], "-^", color=C["ball"], lw=1.0, ms=4, label="ZEB2")
    ax.set_xticks(xx)
    ax.set_xticklabels([f"{s.replace('_', ' ')}\nn={int(nmap.get(s, 0))}" for s in order206],
                       fontsize=5.6)
    ax.set_ylabel("Mean log-normalised expression")
    ax.set_title("GSE206710 (3 donors; second atlas)", loc="left")
    ax.legend(fontsize=5.6, loc="upper right")
    ax.text(0.97, 0.28, "ZEB1 vs stage ρ = −0.50, P = 0.67",
            transform=ax.transAxes, ha="right", fontsize=5.5, color="#555")
    S.style_ax(ax)
    S.panel_label(ax, "C", x=-0.16)

    ax = fig.add_subplot(gs[1, 1])
    x = np.array([0, 1])
    for _, r in cd3.iterrows():
        ax.plot(x, [r["ZEB1_CD3neg"], r["ZEB1_CD3pos"]], "-o", color=C["ink"],
                lw=0.9, ms=4)
        ax.text(-0.08, r["ZEB1_CD3neg"], r["donor"], ha="right", va="center", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(["CD3−", "CD3+"])
    ax.set_xlim(-0.45, 1.35)
    ax.set_ylabel("ZEB1 (donor mean)")
    ax.set_title("GSE206710 paired CD3− / CD3+", loc="left")
    ax.text(0.97, 0.92, "Wilcoxon P = 0.25\nΔ = +0.079 (CD3− higher)",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.8)
    S.style_ax(ax)
    S.panel_label(ax, "D", x=-0.18)

    ax = fig.add_subplot(gs[2, 0])
    confusion = pd.read_csv(_sd("validation/thymus_reference/heldout_confusion.tsv"), sep="\t", index_col=0)
    fractions = confusion.div(confusion.sum(axis=1), axis=0)
    ax.imshow(fractions, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(4):
        for j in range(4):
            ax.text(j,i,str(confusion.iloc[i,j]),ha="center",va="center",fontsize=6,
                    color="white" if fractions.iloc[i,j]>.5 else "#333")
    ax.set_xticks(range(4), confusion.columns, fontsize=5.5)
    ax.set_yticks(range(4), confusion.index, fontsize=5.5)
    ax.set_xlabel("Held-out prediction");ax.set_ylabel("Author reference stage")
    ax.set_title("Leave-one-donor-out: 70 profiles / 20 donors", loc="left", fontsize=6.5)
    S.panel_label(ax,"E",x=-.18)
    ax = fig.add_subplot(gs[2, 1])
    affinity = pd.read_csv(_sd("validation/thymus_reference/query_patient_reference_affinity.tsv"), sep="\t")
    affinity=affinity.set_index("patient").loc[sorted(affinity.patient, key=lambda p:int(p[1:]))]
    mat=affinity[["affinity_"+stage for stage in order]].to_numpy()
    im=ax.imshow(mat,cmap="YlGnBu",vmin=0,vmax=1,aspect="auto")
    ax.set_xticks(range(4),order,fontsize=5.5)
    ax.set_yticks(range(len(affinity)),[p+(" *" if flag else "") for p,flag in
        zip(affinity.index,affinity.outside_heldout_95pct_distance)],fontsize=5.5)
    cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
    cb.set_label("Relative affinity (not probability)",fontsize=5)
    ax.set_title("T-ALL patient profiles: reference affinity",loc="left",fontsize=6.5)
    ax.set_xlabel("* Beyond held-out 95th-percentile distance",fontsize=5.3)
    S.panel_label(ax,"F",x=-.18)
    S.save_supp(fig, "FigureS2_thymus_atlas")


# ===========================================================================
# S3  Subtype extras, CIMP, cis-genetic (M3)
# ===========================================================================
def fig_s3():
    mut = pd.read_csv(_p("module3", "M3_3.11_3.13_mutation_effects.tsv"), sep="\t")
    fus = pd.read_csv(_p("module3", "M3_3.24_ZEB1_contrasts.tsv"), sep="\t")
    cimp_n = pd.read_csv(_p("module3", "M3_3.19_cimp_counts.tsv"), sep="\t")
    cnv = pd.read_csv(_p("module3", "M3_3.14_cnv_burden.tsv"), sep="\t")
    loc = pd.read_csv(_p("module3", "M3_3.16_ZEB1_locus_CNV.tsv"), sep="\t")
    zf = pd.read_csv(_p("module3", "M3_3.17_3.18_ZEB1_mutation_fusion.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 135 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.45)

    ax = fig.add_subplot(gs[0, 0])
    row = fus[fus["contrast"].str.contains("CIMP")].iloc[0]
    ax.bar([0, 1], [row["mean_b"], row["mean_a"]], color=[C["ball"], C["tall"]],
           width=0.55, edgecolor="white", lw=0.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"CIMP-low\nn={int(row['n_b'])}", f"CIMP-high\nn={int(row['n_a'])}"],
                       fontsize=6.2)
    ax.set_ylabel("ZEB1 (VST)")
    ax.set_ylim(0, max(row["mean_a"], row["mean_b"]) * 1.24)
    ax.set_title("ZEB1 by CIMP group (GSE272023)", loc="left")
    ax.text(0.5, 0.92, f"g = {row['hedges_g']:.2f} [{row['ci95_lo']:.2f}, {row['ci95_hi']:.2f}]\nP = {row['wilcox_p']:.2f}",
            transform=ax.transAxes, ha="center", fontsize=6)
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    m = mut.copy()
    m["lab"] = m["contrast"].str.replace(" mut vs wt", "")
    m = m.sort_values("hedges_g")
    y = np.arange(len(m))
    ax.hlines(y, m["ci95_lo"], m["ci95_hi"], color="#B8BCC2", lw=1.0)
    cols = [C["tall"] if p < 0.05 else C["grey"] for p in m["wilcox_p"]]
    ax.scatter(m["hedges_g"], y, s=26, color=cols, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(m["lab"], fontsize=6.2)
    ax.set_xlabel("Hedges g (mut vs wt)")
    ax.set_title("Driver mutations vs ZEB1 (TARGET)", loc="left")
    ax.text(0.97, 0.06, "none pass FDR < 0.05", transform=ax.transAxes,
            ha="right", fontsize=5.8, color="#555")
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.28)

    ax = fig.add_subplot(gs[1, 0])
    f = fus[fus["contrast"].str.contains("vs no_fusion")].drop_duplicates("contrast").copy()
    f["lab"] = f["group_a"]
    f = f.sort_values("hedges_g")
    y = np.arange(len(f))
    ax.hlines(y, f["ci95_lo"], f["ci95_hi"], color="#B8BCC2", lw=1.0)
    cols = [C["down"] if (g < 0 and p < 0.05) else (C["up"] if (g > 0 and p < 0.05) else C["grey"])
            for g, p in zip(f["hedges_g"], f["wilcox_p"])]
    ax.scatter(f["hedges_g"], y, s=22, color=cols, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} (n={int(n)})" for a, n in zip(f["lab"], f["n_a"])], fontsize=5.8)
    ax.set_xlabel("Hedges g vs no fusion")
    ax.set_title("Fusion class vs ZEB1", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "C", x=-0.42)

    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("Cis-genetic summary", loc="left")
    lines = [
        ("ZEB1 coding mutations", f"{int(zf.loc[zf.gene=='ZEB1','n_mut'].iloc[0])} / {int(zf.loc[zf.gene=='ZEB1','n'].iloc[0])}"),
        ("ZEB1 fusions", f"{int(zf.loc[zf.gene=='ZEB1','n_fusion'].iloc[0])} / {int(zf.loc[zf.gene=='ZEB1','n'].iloc[0])}"),
        ("ZEB1 locus CN vs RNA", f"ρ = {float(loc['rho'].iloc[0]):.2f}, P = {float(loc['p'].iloc[0]):.2f}"),
        ("CNV burden vs ZEB1", f"ρ = {float(cnv['rho'].iloc[0]):.2f}, P = {float(cnv['p'].iloc[0]):.2f}"),
        ("CIMP-high / CIMP-low", f"{int(cimp_n.loc[cimp_n.cimp_g=='CIMP-high','N'].iloc[0])} / "
         f"{int(cimp_n.loc[cimp_n.cimp_g=='CIMP-low','N'].iloc[0])}"),
        ("CIMP effect (95% CI)", f"{row['hedges_g']:.2f} [{row['ci95_lo']:.2f}, {row['ci95_hi']:.2f}]"),
    ]
    yy = 0.88
    for lab, val in lines:
        ax.text(0.04, yy, lab, fontsize=6.4, transform=ax.transAxes, color="#444")
        ax.text(0.96, yy, val, fontsize=6.6, transform=ax.transAxes, ha="right",
                fontweight="bold", color=C["ink"])
        yy -= 0.13
    ax.text(0.04, 0.04, "Denominator: RNA-linked records;\nvariant callability not verified per sample.",
            fontsize=6.2, style="italic", color=C["tall"], transform=ax.transAxes)
    S.panel_label(ax, "D", x=-0.05)

    S.save_supp(fig, "FigureS3_subtype_CIMP_cis")


# ===========================================================================
# S4  LMO2–ZEB1 subtype-conditional cohort heterogeneity
# ===========================================================================
def fig_s4():
    edge = pd.read_csv(_sd("SourceData_FigS4_LMO2_ZEB1_false_edge.tsv"), sep="\t")
    adj = pd.read_csv(_p("module4", "M4_r2_keygenes_adjusted.tsv"), sep="\t")

    fig = plt.figure(figsize=(183 * MM, 183 * MM), constrained_layout=True)
    gs = GridSpec(3, 2, figure=fig, hspace=0.16, wspace=0.16,
                  height_ratios=[1.05, 0.85, 1.0])

    # A observed effect matrix; rows are independent cohorts, not causal paths
    ax = fig.add_subplot(gs[0, 0])
    matrix = edge.set_index("cohort")[["rho_raw", "rho_adj"]]
    order = ["TARGET-ALL-P2", "GSE62156", "GSE26713", "GSE110636"]
    matrix = matrix.reindex(order)
    im = ax.imshow(matrix.values, cmap="RdBu_r", vmin=-.65, vmax=.65, aspect="auto")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Unadjusted", "Subtype residual"], fontsize=6)
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels([_cohort_label(c) for c in matrix.index], fontsize=6)
    for i in range(len(matrix)):
        for j in range(2):
            ax.text(j, i, f"{matrix.iloc[i, j]:+.2f}", ha="center", va="center",
                    fontsize=6.5, color="white" if abs(matrix.iloc[i, j]) > .38 else "#222")
    cb = fig.colorbar(im, ax=ax, fraction=.055, pad=.03)
    cb.set_label("Spearman ρ", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("Observed associations differ across cohorts", loc="left", pad=8)
    S.panel_label(ax, "A", x=-0.02, y=1.08)

    # B forest raw vs adj
    ax = fig.add_subplot(gs[0, 1])
    e = edge.set_index("cohort").reindex(order).reset_index()
    e["lab"] = e["cohort"].map(_cohort_label) + e["n"].map(lambda n: f"\nn={int(n)}")
    y = np.arange(len(e))[::-1]
    ax.hlines(y + 0.14, e["ci_raw_lo"], e["ci_raw_hi"], color=C["grey"], lw=1.0)
    ax.hlines(y - 0.14, e["ci_adj_lo"], e["ci_adj_hi"], color=C["tall"], lw=1.0)
    ax.scatter(e["rho_raw"], y + 0.14, s=22, color=C["grey"], zorder=3, label="unadjusted")
    ax.scatter(e["rho_adj"], y - 0.14, s=26, color=C["tall"], zorder=3, label="subtype-residual")
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(e["lab"], fontsize=6)
    ax.set_xlabel("Spearman ρ (LMO2 vs ZEB1)")
    ax.set_title("Cohort estimates and intervals", loc="left", pad=25)
    ax.legend(loc="lower center", fontsize=5.6, ncol=2, bbox_to_anchor=(0.5, 1.01))
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.32)

    # C TARGET key-gene t unadj vs adj (genome-supported)
    ax = fig.add_subplot(gs[1, :])
    order = ["GATA3", "TCF7", "BCL11B", "IL7R", "ZEB2", "LYL1", "MEF2C", "LMO2"]
    a = adj.set_index("gene").reindex(order).dropna(subset=["t_unadj"])
    y = np.arange(len(a))
    ax.hlines(y, a["t_unadj"], a["t_adj"], color="#B8BCC2", lw=1.0)
    ax.scatter(a["t_unadj"], y, s=20, color=C["grey"], zorder=2, label="unadjusted")
    cols = [C["down"] if f >= 0.05 else C["ink"] for f in a["fdr_adj"]]
    ax.scatter(a["t_adj"], y, s=26, color=cols, zorder=3, label="subtype+age adj.")
    ax.axvline(0, color="#333", lw=0.6, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(a.index, fontsize=6.2)
    ax.set_xlabel("moderated t (TARGET)")
    ax.set_title("Selected-gene statistics before and after adjustment", loc="left")
    ax.legend(loc="upper right", fontsize=5.5)
    ax.annotate("LMO2 FDR = 0.18", xy=(a.loc["LMO2", "t_adj"],
                list(a.index).index("LMO2")),
                xytext=(6, -8), textcoords="offset points", fontsize=5.6, color=C["down"])
    S.style_ax(ax)
    S.panel_label(ax, "C", x=-0.02)

    # Reuse the user's completed TARGET subgroup export (no regression refit).
    subdir = "validation/lmo2_subgroups_existing/"
    patient = pd.read_csv(_sd(subdir + "TALL_clinical_annotation_full.csv"))
    results = pd.read_csv(_sd(subdir + "TALL_LMO2_ZEB1_correlation_by_group.csv")).set_index("Cohort")
    is_lmo = patient.MolecularSubtype.isin(["LMO1/2", "LMO2_LYL1"])
    palette = {"LMO1/2": C["aml"], "LMO2_LYL1": C["purple"],
               "TAL1": C["ball"], "TAL2": "#56B4E9", "TLX1": C["tall"],
               "TLX3": "#CC79A7", "HOXA": C["green"], "NKX2_1": C["teal"],
               "Unknown": C["grey"]}
    for j, (group, mask, letter) in enumerate([
            ("LMO2-associated", is_lmo, "D"),
            ("non-LMO2-associated", ~is_lmo, "E")]):
        ax = fig.add_subplot(gs[2, j])
        d = patient.loc[mask]
        r = results.loc[group]
        assert len(d) == int(r["n"])
        for subtype, rows in d.groupby("MolecularSubtype", sort=True):
            ax.scatter(rows.LMO2, rows.ZEB1, s=13, alpha=.8,
                       color=palette[subtype], edgecolor="white", linewidth=.25,
                       label=subtype.replace("_", "/"), zorder=3)
        ax.set_xlabel("LMO2 log2(TPM + 1)")
        ax.set_ylabel("ZEB1 log2(TPM + 1)")
        ax.set_title(f"TARGET: {group} (n={len(d)})", loc="left", fontsize=6.8)
        ax.text(.03, .97, f"ρ = {r.Spearman_R:+.2f}; P = {r.P_value:.3f}",
                transform=ax.transAxes, va="top", fontsize=5.8)
        ax.set_xlim(patient.LMO2.min()-.4, patient.LMO2.max()+.4)
        ax.set_ylim(patient.ZEB1.min()-.4, patient.ZEB1.max()+1.0)
        ax.legend(loc="lower center", bbox_to_anchor=(.5, -.48),
                  ncol=2 if j == 0 else 4, fontsize=4.8, columnspacing=.6,
                  handletextpad=.2)
        S.style_ax(ax); S.panel_label(ax, letter, x=-.16)
    S.save_supp(fig, "FigureS4_LMO2_ZEB1_false_edge")


# ===========================================================================
# S5  Complete GSEA (M4)
# ===========================================================================
def fig_s5():
    hall = pd.read_csv(_p("module4", "M4_4.6_hallmark_fgsea.tsv"), sep="\t")
    rea = pd.read_csv(_p("module4", "M4_4.7_reactome_fgsea.tsv"), sep="\t")
    go = pd.read_csv(_p("module4", "M4_4.8_gobp_fgsea.tsv"), sep="\t")
    mod = pd.read_csv(_p("module4", "M4_4.28_module_vs_ZEB1.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 145 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.52, wspace=0.72)

    ax = fig.add_subplot(gs[0, 0])
    h = hall.copy()
    h = h.reindex(h["NES"].abs().sort_values(ascending=False).index).head(16).sort_values("NES")
    y = np.arange(len(h))
    cols = [C["tall"] if p < 0.05 else C["ns"] for p in h["padj"]]
    ax.hlines(y, 0, h["NES"], color="#C7CBD1", lw=0.9)
    ax.scatter(h["NES"], y, s=16, color=cols, edgecolor="#333", lw=0.3, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([_short(p, 26) for p in h["pathway"]], fontsize=5.0)
    ax.set_xlabel("NES")
    ax.set_title("Hallmark (0 / 50 FDR < 0.05)", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.55)

    ax = fig.add_subplot(gs[0, 1])
    r = rea.copy()
    r = r.reindex(r["NES"].abs().sort_values(ascending=False).index).head(12).sort_values("NES")
    y = np.arange(len(r))
    dots = ax.scatter(r["NES"], y, c=-np.log10(r["padj"].clip(lower=1e-300)),
                      cmap="viridis", s=24, edgecolor="#333", lw=.25, zorder=3)
    cb = fig.colorbar(dots, ax=ax, fraction=.04, pad=.02)
    cb.set_label("-log10 BH FDR", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([_short(p, 32) for p in r["pathway"]], fontsize=4.8)
    ax.set_xlabel("NES")
    ax.set_title("Reactome enrichment", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.72)

    ax = fig.add_subplot(gs[1, 0])
    g = go.copy()
    g = g.reindex(g["NES"].abs().sort_values(ascending=False).index).head(12).sort_values("NES")
    y = np.arange(len(g))
    dots = ax.scatter(g["NES"], y, c=-np.log10(g["padj"].clip(lower=1e-300)),
                      cmap="viridis", s=24, edgecolor="#333", lw=.25, zorder=3)
    cb = fig.colorbar(dots, ax=ax, fraction=.04, pad=.02)
    cb.set_label("-log10 BH FDR", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([_short(p, 32) for p in g["pathway"]], fontsize=4.8)
    ax.set_xlabel("NES")
    ax.set_title("GO biological processes", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "C", x=-0.72)

    ax = fig.add_subplot(gs[1, 1])
    y = np.arange(len(mod))
    ax.scatter(mod["rho"], y, s=36, color=C["purple"], edgecolor="white", lw=.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{m} (n={int(n)})" for m, n in zip(mod["module"], mod["n_genes"])],
                       fontsize=6.2)
    ax.set_xlabel("Spearman ρ with ZEB1")
    ax.set_xlim(0.70, 0.95)
    ax.set_title("Co-expression module associations", loc="left")
    ax.text(0.97, 0.12, "all ρ 0.84–0.87", transform=ax.transAxes,
            ha="right", fontsize=5.8, color="#555")
    S.style_ax(ax)
    S.panel_label(ax, "D", x=-0.32)

    S.save_supp(fig, "FigureS5_complete_GSEA")


# ===========================================================================
# S6  Single-cell QC and corrected malignant-cell analysis (§3.5)
# ===========================================================================
def fig_s6():
    comp = pd.read_csv(_p("module5", "M5_r2_composition_by_ZEB1bin.tsv"), sep="\t")
    disc = pd.read_csv(_p("module5", "M5_5.41_patient_ZEB1_associations.tsv"), sep="\t")
    qc = pd.read_csv(_sd("gse248287_malignant_job1699/patient_cell_qc.tsv"), sep="\t")
    comp_sens = pd.read_csv(_sd("SourceData_Fig5G_composition_sensitivity.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 140 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.45)

    ax = fig.add_subplot(gs[0, 0])
    c = comp.set_index("zeb_bin")
    stages = [x for x in c.columns]
    pal = [C["purple"], C["aml"], C["ball"]]
    bottom = np.zeros(len(c))
    x = np.arange(len(c))
    for i, stg in enumerate(stages):
        v = c[stg].to_numpy()
        ax.bar(x, v, bottom=bottom, color=pal[i % 3], width=0.55, label=stg,
               edgecolor="white", lw=0.3)
        bottom = bottom + v
    ax.set_xticks(x)
    ax.set_xticklabels(c.index.str.replace("_", " "), fontsize=6.2)
    ax.set_ylabel("Cells")
    ax.set_title("Lineage composition by ZEB1 bin", loc="left")
    ax.legend(fontsize=5.6, loc="upper right")
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.18)

    ax = fig.add_subplot(gs[0, 1])
    stability = pd.read_csv(_sd("validation/harmony_label_sensitivity/patient_scores_and_overlap.tsv"), sep="\t")
    patient_order = sorted(stability.patient.unique(), key=lambda p: int(p[1:]))
    matrix = stability.pivot(index="patient", columns=["run_id", "resolution"], values="jaccard").loc[patient_order]
    matrix = matrix.sort_index(axis=1)
    im = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(4), ["Run 1\nr=0.6", "Run 1\nr=1.0", "Run 2\nr=0.6", "Run 2\nr=1.0"], fontsize=5.5)
    ax.set_yticks(range(len(matrix)), matrix.index, fontsize=5.6)
    for i in range(len(matrix)):
        for j in range(4):
            v = matrix.iloc[i,j]
            ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=5.5,
                    color="white" if v>.65 else "#222")
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.03)
    cb.set_label("Malignant-label Jaccard",fontsize=5.5)
    ax.set_title("Extended Harmony: patient label overlap", loc="left")
    S.panel_label(ax, "B", x=-0.18)

    ax = fig.add_subplot(gs[1, 0])
    q = qc.sort_values("malignant_fraction")
    y = np.arange(len(q))
    ax.barh(y, q["malignant_fraction"], color=C["tall"], height=0.7,
            edgecolor="white", lw=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(q["patient"], fontsize=5.6)
    ax.set_xlabel("Malignant fraction in sequenced mixture")
    ax.set_title("GSE248287 diagnosis cell composition", loc="left")
    ax.set_xlim(0, 1)
    S.style_ax(ax)
    S.panel_label(ax, "C", x=-0.42)

    ax = fig.add_subplot(gs[1, 1])
    m = comp_sens
    ax.scatter(m["rho_mixed"], m["rho_malignant"], s=32, c=C["tall"],
               edgecolor="#333", lw=0.35, zorder=3)
    for _, r in m.iterrows():
        offset = (4, 10) if r["feature"] == "TCF7" else ((4, -12) if r["feature"] == "GATA3prog" else (3, 3))
        ax.annotate(r["feature"], (r["rho_mixed"], r["rho_malignant"]),
                    fontsize=4.7, xytext=offset, textcoords="offset points")
    ax.axhline(0, color="#999", lw=0.6, ls=":")
    ax.axvline(0, color="#999", lw=0.6, ls=":")
    ax.set_xlabel("ρ, mixed cell types (previous)")
    ax.set_ylabel("ρ, malignant-only pseudobulk")
    ax.set_xlim(-1, 1.08)
    ax.set_title("Mixed versus malignant-cell associations", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "D", x=-0.18)
    S.save_supp(fig, "FigureS6_scrna_composition_sensitivity")


# ===========================================================================
# S7  CIMP / methylation / occupancy extras (M6)
# ===========================================================================
def fig_s7():
    probes = pd.read_csv(_p("module6", "M6_6.17_ZEB1_probe_CIMP_stats.tsv"), sep="\t")
    occ = pd.read_csv(_p("module6", "M6_r2_GSE154675_occupancy.tsv"), sep="\t")
    audit_occ = pd.read_csv(_sd("validation/chromatin_windows/window_summaries.tsv"),sep="\t")
    occ = audit_occ[audit_occ.window.eq("legacy_12kb")].rename(columns={"mean_covered":"hg38_prom"})
    k27 = pd.read_csv(_p("module6", "M6_r2_GSE70734_H3K27ac_ZEB1.tsv"), sep="\t")
    summ = pd.read_csv(_p("module6", "M6_6.17_ZEB1_promoter_beta_summary.tsv"), sep="\t")

    fig = plt.figure(figsize=(183 * MM, 158 * MM), layout="constrained")
    gs = GridSpec(2, 3, figure=fig, wspace=0.12, hspace=.10, height_ratios=[1,1.3])

    ax = fig.add_subplot(gs[0, 0])
    pr = probes[probes["region"] == "promoter"].dropna(subset=["mean_CIMP_pos"]).copy()
    pr = pr.sort_values("pos")
    x = np.arange(len(pr))
    ax.scatter(x, pr["mean_CIMP_neg"], s=12, color=C["ball"], label="CIMP−", zorder=3)
    ax.scatter(x, pr["mean_CIMP_pos"], s=12, color=C["tall"], label="CIMP+", zorder=3)
    ax.plot(x, pr["mean_CIMP_neg"], color=C["ball"], lw=0.6, alpha=0.6)
    ax.plot(x, pr["mean_CIMP_pos"], color=C["tall"], lw=0.6, alpha=0.6)
    ax.axhline(0.5, color="#999", lw=0.6, ls=":")
    ax.set_ylim(0, 0.7)
    ax.set_xlabel("Promoter probes (TSS / 5′UTR / 1st exon)")
    ax.set_ylabel("Mean β")
    ax.set_title("Promoter methylation", loc="left")
    ax.legend(fontsize=5.6)
    ax.set_xticks([])
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.16)
    _ = summ, k27

    ax = fig.add_subplot(gs[0, 1])
    gb = probes[probes["region"] == "gene_body"].dropna(subset=["delta_pos_minus_neg"]).copy()
    gb = gb.sort_values("delta_pos_minus_neg")
    y = np.arange(len(gb))
    cols = [C["tall"] if (q < 0.05 if pd.notna(q) else False) else C["grey"] for q in gb["q_bh"]]
    ax.hlines(y, 0, gb["delta_pos_minus_neg"], color=cols, lw=.8)
    ax.scatter(gb["delta_pos_minus_neg"], y, color=cols, s=18, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(gb["probe"], fontsize=5.0)
    ax.set_xlabel("Δβ (CIMP+ − CIMP−)")
    ax.set_title("Gene-body methylation", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.32)

    ax = fig.add_subplot(gs[0, 2])
    o = occ.copy()
    lines = sorted(o["line"].unique())
    abs_ = ["LMO2", "TAL1", "LDB1", "GATA2"]
    M = np.array([[float(o.loc[(o.line == ln) & (o.antibody == ab), "hg38_prom"].iloc[0])
                   if ((o.line == ln) & (o.antibody == ab)).any() else np.nan
                   for ab in abs_] for ln in lines])
    im = ax.imshow(M, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(abs_)))
    ax.set_xticklabels(abs_, fontsize=6.2)
    ax.set_yticks(range(len(lines)))
    ax.set_yticklabels(lines, fontsize=6.2)
    ax.set_title("Shared-region ChIP signal", loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("signal", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.text(0.03, -0.22, "Occupancy ≠ repression", transform=ax.transAxes,
            fontsize=5.6, style="italic", color="#555")
    S.panel_label(ax, "C", x=-0.22)

    ax = fig.add_subplot(gs[1, :])
    bins = pd.read_csv(_sd("validation/chromatin_windows/locus_200bp_bins.tsv"), sep="\t")
    bins["track"] = bins.line + " / " + bins.antibody
    order = [line + " / " + ab for line in ["ARR","DU528","HSB2","CCRFCEM"] for ab in abs_]
    hm = bins.pivot(index="track",columns="relative_mid_bp",values="mean_missing_as_zero").loc[order]
    transformed = np.log1p(hm.to_numpy())
    im = ax.imshow(transformed,aspect="auto",cmap="magma",extent=[-10,10,len(order)-.5,-.5])
    ax.set_yticks(range(len(order)),order,fontsize=5.4)
    ax.set_xlabel("Distance from canonical ZEB1 TSS (kb; GRCh38 chr10:31,319,216)")
    ax.set_title("Local ChIP profiles across the ZEB1 / ZEB1-AS1 overlap",loc="left",fontsize=7.5,pad=32)
    for y in [3.5,7.5,11.5]: ax.axhline(y,color="white",lw=.9)
    for x in [-1,1]: ax.axvline(x,color="#CCCCCC",lw=.7,ls="--")
    ax.axvline(0,color="white",lw=.6)
    top = ax.inset_axes([0,1.015,1,.085]); top.set_xlim(-10,10);top.set_ylim(0,2);top.set_axis_off()
    top.axvspan(-2.768,9.232,color="#BFC5CC",alpha=.4)
    top.annotate("",xy=(9.5,1.5),xytext=(-.721,1.5),arrowprops=dict(arrowstyle="->",color=C["tall"],lw=1.2))
    top.text(4.5,1.6,"ZEB1 (+)",ha="center",fontsize=5.5,color=C["tall"])
    top.annotate("",xy=(-9.5,.4),xytext=(1.269,.4),arrowprops=dict(arrowstyle="->",color=C["ball"],lw=1.2))
    top.text(-5,.5,"ZEB1-AS1 (−)",ha="center",fontsize=5.5,color=C["ball"])
    cb=fig.colorbar(im,ax=ax,fraction=.024,pad=.015)
    cb.set_label("log(1 + mean track signal)",fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    S.plabel(ax,"D",dy=32)

    S.save_supp(fig, "FigureS7_CIMP_occupancy")


# ===========================================================================
# S8  Eighteen drugs profiled, seventeen analyzable (M7)
# ===========================================================================
def fig_s8():
    drug = pd.read_csv(_p("module7", "M7_7.7_7.8_drug_spearman.tsv"), sep="\t")
    miss = pd.read_csv(_p("module7", "M7_r2_missingness.tsv"), sep="\t")
    klass = pd.read_csv(_p("module7", "M7_r2_ZEB1_rho_with_class.tsv"), sep="\t")
    lm = pd.read_csv(_p("module7", "M7_r2_core_drug_lm.tsv"), sep="\t")

    fig = plt.figure(figsize=(183 * MM, 168 * MM),layout="constrained")
    gs = GridSpec(2, 2, figure=fig, hspace=0.16, wspace=0.16,height_ratios=[1,1.15])

    ax = fig.add_subplot(gs[0, 0])
    z = drug[drug.gene == "ZEB1_log"].sort_values("rho")
    assert len(z) == 17 and z["n"].min() == 14 and z["n"].max() == 94
    y = np.arange(len(z))
    cols = [C["tall"] if f < 0.05 else C["grey"] for f in z["fdr"]]
    ax.hlines(y, 0, z["rho"], color="#C7CBD1", lw=0.9)
    ax.scatter(z["rho"], y, s=[max(8, n / 2.2) for n in z["n"]], color=cols,
               edgecolor="#333", lw=0.35, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d}  n={int(n)}" for d, n in zip(z["drug_short"], z["n"])],
                       fontsize=5.4)
    ax.set_xlabel("Spearman ρ (ZEB1 vs LC50)")
    ax.set_title("ZEB1–drug correlations", loc="left")
    S.style_ax(ax)
    S.plabel(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    ms = miss.sort_values("frac_missing")
    y = np.arange(len(ms))
    ax.barh(y, ms["frac_missing"], color=C["grey"], height=0.7, edgecolor="white", lw=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(ms["drug_short"], fontsize=5.4)
    ax.set_xlabel("Fraction missing LC50")
    ax.set_title("LC50 availability | 18 drugs", loc="left")
    ax.axvline(0.5, color=C["tall"], lw=0.7, ls="--")
    S.style_ax(ax)
    S.plabel(ax, "B")

    ax = fig.add_subplot(gs[1, 0])
    k = klass.sort_values(["class","drug_short"])
    gene_order = ["ZEB1_log","ZEB2_log","LMO2_log"]
    existing = drug[drug.gene.isin(gene_order)].copy()
    existing.to_csv(_sd("SourceData_FigS8C_gene_drug_context.tsv"),sep="\t",index=False)
    matrix = existing.pivot(index="drug_short",columns="gene",values="rho").loc[k.drug_short,gene_order]
    assert matrix.shape==(17,3) and np.isfinite(matrix.to_numpy()).all()
    im = ax.imshow(matrix,cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto")
    ax.set_yticks(range(len(k)),[f"{name} · {category}" for name,category in zip(k.drug_short,k['class'])],fontsize=5.1)
    ax.set_xticks([0,1,2],["ZEB1","ZEB2","LMO2"],fontsize=6)
    for yi,row in enumerate(matrix.to_numpy()):
        for xi,val in enumerate(row):
            ax.text(xi,yi,f"{val:+.2f}",ha="center",va="center",fontsize=5.1,color="white" if abs(val)>.7 else "#222")
    classes=k['class'].tolist()
    for i in range(1,len(classes)):
        if classes[i]!=classes[i-1]: ax.axhline(i-.5,color="white",lw=1)
    ax.tick_params(length=0)
    ax.set_title("Three-gene drug-response matrix",loc="left",fontsize=7)
    cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.03)
    cb.set_label("Spearman ρ",fontsize=6);cb.ax.tick_params(labelsize=5)
    S.plabel(ax,"C")

    ax = fig.add_subplot(gs[1, 1])
    zlm = lm[lm["model"] == "ZEB1"].copy()
    zlm = zlm.sort_values("p", ascending=False)
    y = np.arange(len(zlm))
    cols = [C["tall"] if f < 0.05 else C["grey"] for f in zlm["fdr"]]
    ax.scatter(np.clip(zlm["p"], 1e-4, 1), y, s=18, color=cols, zorder=3)
    ax.axvline(0.05, color=C["tall"], lw=0.7, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(zlm["drug_short"], fontsize=5.2)
    ax.set_xlabel("OLS P (ZEB1 term)")
    ax.set_title("ZEB1 coefficient (univariate OLS)", loc="left")
    S.style_ax(ax)
    S.plabel(ax, "D")

    S.save_supp(fig, "FigureS8_complete_18drug")


# ===========================================================================
# S9  Complete DepMap (M8)
# ===========================================================================
def fig_s9():
    dep = pd.read_csv(_sd("SourceData_Fig7C_ZEB1_dependency.tsv"), sep="\t")
    tgt = pd.read_csv(_sd("SourceData_Fig7DE_drug_target_GE.tsv"), sep="\t")
    assert dep["ZEB1_GE"].notna().sum() == 8 and tgt["n_TALL"].eq(8).all()

    fig = plt.figure(figsize=(183 * MM, 135 * MM),layout="constrained")
    gs = GridSpec(2, 2, figure=fig, hspace=0.12, wspace=0.16,height_ratios=[1,.75])

    left = gs[0,0].subgridspec(1,2,width_ratios=[4,1],wspace=.05)
    ax = fig.add_subplot(left[0,0])
    d = dep.dropna(subset=["ZEB1_GE"]).sort_values("ZEB1_GE")
    y = np.arange(len(d))
    cols = [C["tall"] if v < -0.5 else C["grey"] for v in d["ZEB1_GE"]]
    ax.hlines(y, d["ZEB1_GE"], 0, color="#CBD0D6", lw=1)
    ax.scatter(d["ZEB1_GE"], y, color=cols, s=28, edgecolor="white", lw=.3, zorder=3)
    ax.axvline(-0.5, color=C["down"], lw=0.8, ls="--")
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(d["CellLineName"], fontsize=5.6)
    ax.set_xlabel("ZEB1 CRISPR gene effect")
    ax.set_ylim(-.5,len(d)-.5)
    ax.set_title("ZEB1 dependency | 8 T-ALL models", loc="left",fontsize=7)
    S.style_ax(ax)
    S.plabel(ax,"A")
    prob = fig.add_subplot(left[0,1])
    vals=d['ZEB1_dependency'].to_numpy()
    assert np.isfinite(vals).all() and ((vals>=0)&(vals<=1)).all()
    prob.imshow(vals[:,None],origin='lower',cmap='YlOrRd',vmin=0,vmax=1,aspect='auto')
    for yi,v in enumerate(vals):
        prob.text(0,yi,f"{v:.2f}",ha='center',va='center',fontsize=5.5,color='white' if v>.65 else '#222')
    prob.set_xticks([]); prob.set_yticks([])
    prob.set_title("P(dep.)",fontsize=6)
    prob.set_xlabel("0–1",fontsize=5.5)
    d.to_csv(_sd('SourceData_FigS9A_effect_and_probability.tsv'),sep='\t',index=False)

    ax = fig.add_subplot(gs[0, 1])
    t = dep[dep["ZEB1_RNA"].notna() & dep["ZEB1_GE"].notna()]
    ax.scatter(t["ZEB1_RNA"], t["ZEB1_GE"], s=28, color=C["tall"], edgecolor="#333", lw=0.35)
    for _, r in t.iterrows():
        offsets={'CCRF-HSB-2-DM':(-3,-13),'KE-37':(-3,7),'PF-382':(-8,9),
                 'PEER':(5,-4),'DND-41':(5,2),'JURKAT':(-30,-12),
                 'SUP-T1':(4,5),'TALL-1':(4,-9)}
        ax.annotate(str(r["CellLineName"]), (r["ZEB1_RNA"], r["ZEB1_GE"]),
                    fontsize=5.2, xytext=offsets.get(str(r['CellLineName']),(3,2)), textcoords="offset points")
    ax.margins(x=.22,y=.16)
    ax.axhline(-0.5, color=C["down"], lw=0.7, ls="--")
    ax.set_xlabel("ZEB1 RNA log2(TPM+1)")
    ax.set_ylabel("ZEB1 gene effect")
    ax.set_title("Expression vs self-dependency", loc="left")
    S.style_ax(ax)
    S.plabel(ax,"B")

    ax = fig.add_subplot(gs[1, :])
    tt = tgt.sort_values("mean_TALL").set_index("gene")
    matrix = tt[["mean_TALL", "mean_BALL", "mean_myeloid"]].T
    limit = max(abs(float(matrix.min().min())), abs(float(matrix.max().max())))
    im = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-limit, vmax=limit,
                   aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(tt)), tt.index, rotation=90, fontsize=5.4)
    ax.set_yticks(range(3), ["T-ALL (8)", "B-ALL (13)", "Myeloid (45)"], fontsize=6)
    ax.tick_params(length=0)
    ax.set_title("Selected-gene dependency across lineages", loc="left",fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=.025, pad=.02)
    cb.set_label("Mean CRISPR gene effect", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    S.plabel(ax,"C")
    S.save_supp(fig, "FigureS9_complete_DepMap")


# ===========================================================================
# S10  Perturbation extras + CollecTRI TF activity  (§3.4)
# ===========================================================================
def fig_s10():
    kd = pd.read_csv(_p("module9", "M9_keygene_effects.tsv"), sep="\t")
    lfc = pd.read_csv(_p("module9", "M9_9.24_evidence_LFC.tsv"), sep="\t")
    dc = pd.read_csv(_sd("F6_decoupler_ulm_focus.tsv"), sep="\t")
    mo = pd.read_csv(_p("module9", "M9_gse188225_keygenes.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 140 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.45)

    ax = fig.add_subplot(gs[0, 0])
    k = kd[kd.experiment == "gse110635_combined"].copy().sort_values("log2FoldChange")
    y = np.arange(len(k))
    cols = [C["down"] if (v < 0 and p < 0.05) else (C["up"] if (v > 0 and p < 0.05) else C["ns"])
            for v, p in zip(k["log2FoldChange"], k["padj"])]
    ax.hlines(y, 0, k["log2FoldChange"], color=cols, lw=.7)
    ax.scatter(k["log2FoldChange"], y, c=cols, s=24, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(k["gene"], fontsize=6.2)
    ax.set_xlabel("log2FC (TLX1 KD vs NTC)")
    ax.set_title("GSE110635: TLX1 KD raises ZEB1", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.22)

    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(lfc["gse110635_s1"], lfc["gse110635_s2"],
               s=22, color=C["ball"], edgecolor="#333", lw=0.3, zorder=3)
    keep = {"TLX1", "ZEB1", "GATA3", "IL7R", "DNTT", "RAG1", "RUNX1"}
    for _, r in lfc.iterrows():
        if r["gene"] in keep and pd.notna(r["gse110635_s1"]) and pd.notna(r["gse110635_s2"]):
            ax.annotate(r["gene"], (r["gse110635_s1"], r["gse110635_s2"]),
                        fontsize=5.2, xytext=(3, 2), textcoords="offset points")
    ax.axhline(0, color="#999", lw=0.6, ls=":")
    ax.axvline(0, color="#999", lw=0.6, ls=":")
    ax.set_xlabel("siTLX1-1 log2FC")
    ax.set_ylabel("siTLX1-2 log2FC")
    ax.set_title("Both siRNA sequences show ZEB1 up", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.18)

    ax = fig.add_subplot(gs[1, 0])
    audit = pd.read_csv(_sd("validation/regulon_coverage/activity_contrasts.tsv"), sep="\t")
    cov = pd.read_csv(_sd("validation/regulon_coverage/coverage_all.tsv"), sep="\t").set_index("source")
    order = dc.sort_values("delta_kd_minus_ctrl").tf.tolist()
    h = audit.pivot(index="tf", columns="contrast", values="delta_kd_minus_ctrl").loc[order, ["siTLX1_1", "siTLX1_2", "pooled"]]
    limit = float(np.abs(h.to_numpy()).max())
    im = ax.imshow(h.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax.set_yticks(np.arange(len(order)), [f"{tf}  {int(cov.loc[tf, 'detected_any'])}/{int(cov.loc[tf, 'network_targets'])}" for tf in order], fontsize=5.3)
    ax.set_xticks([0, 1, 2], ["siRNA 1", "siRNA 2", "Pooled"], fontsize=6)
    ax.set_title("Activity change / target coverage", loc="left")
    ax.set_xlabel("Labels: targets detected in ≥1 sample / network", fontsize=5.5)
    cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.025)
    cb.set_label("Δ ULM (KD − control)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    S.panel_label(ax, "C", x=-0.22)

    ax = fig.add_subplot(gs[1, 1])
    effects = pd.read_csv(_sd("SourceData_Fig6E_perturbation_forest.tsv"), sep="\t")
    mm = effects[effects.dataset.eq("GSE188225")].copy()
    positions = np.arange(len(mm))
    ax.axhline(0, color="#999", lw=.7, ls=":")
    ax.errorbar(positions, mm.log2FC,
                yerr=np.vstack([mm.log2FC-mm.ci_lo, mm.ci_hi-mm.log2FC]),
                fmt="o", color=C["tall"], capsize=3, markersize=4, lw=1)
    ax.set_xticks(positions, ["Young\nDN2", "Young\nDN3a", "Old\nDN3a"], fontsize=6)
    ax.set_xlim(-.5, 2.5)
    ax.set_ylabel("Zeb1 log2FC (Lmo2-TG versus WT)")
    ax.set_title("GSE188225: age/stage-specific contrasts", loc="left")
    ax.text(.03, .98, "One study; contrasts are not independent cohorts",
            transform=ax.transAxes, va="top", fontsize=5.2)
    S.style_ax(ax)
    S.panel_label(ax, "D", x=-0.18)

    S.save_supp(fig, "FigureS10_perturbation_TF")


# ===========================================================================
# S11  GDSC MTX + pharmacotype proxies + DepMap targets  (§3.3)
# ===========================================================================
def fig_s11():
    pts = pd.read_csv(_sd("F7_mtx_points.tsv"), sep="\t")
    sens = pd.read_csv(_sd("validation/mtx_screen_sensitivity.tsv"), sep="\t")
    paired = pd.read_csv(_sd("validation/mtx_model_overlap.tsv"), sep="\t")
    loo = pd.read_csv(_sd("validation/mtx_leave_one_model_out.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 140 * MM))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.55)

    for j, screen in enumerate(["GDSC1", "GDSC2"]):
        ax = fig.add_subplot(gs[0, j])
        t = pts[(pts.source == screen) & (pts.lineage == "T-ALL")].dropna(subset=["ZEB1", "value"])
        r = sens.set_index("screen").loc[screen]
        ax.scatter(t.ZEB1, t.value, s=30, color=C["tall"], edgecolor="#333", lw=0.4)
        for _, point in t.iterrows():
            if point.cell_line in {"LOUCY", "TALL-1", "SUP-T1", "MOLT-16"}:
                ax.annotate(point.cell_line, (point.ZEB1, point.value), xytext=(3, 2),
                            textcoords="offset points", fontsize=4.5)
        ax.set_xlabel("ZEB1 RNA (log2 TPM+1)")
        ax.set_ylabel("MTX LN_IC50")
        ax.set_title(f"{screen}: n={int(r.n_models)}, ρ={r.rho:.2f}", loc="left")
        ax.text(0.03, 0.97, f"permutation P={r.permutation_p:.3g}\n"
                f"bootstrap 95% CI [{r.bootstrap_ci_lo:.2f}, {r.bootstrap_ci_hi:.2f}]",
                transform=ax.transAxes, fontsize=5.1, va="top")
        S.style_ax(ax); S.panel_label(ax, "A" if j == 0 else "B", x=-0.22)

    ax = fig.add_subplot(gs[1, 0])
    shared = paired.dropna(subset=["GDSC1", "GDSC2"])
    ax.scatter(shared.GDSC1, shared.GDSC2, s=30, color=C["purple"],
               edgecolor="#333", lw=0.4)
    for _, point in shared.iterrows():
        if point.cell_line in {"LOUCY", "SUP-T1", "MOLT-16"}:
            ax.annotate(point.cell_line, (point.GDSC1, point.GDSC2), xytext=(3, 2),
                        textcoords="offset points", fontsize=4.5)
    ax.set_xlabel("GDSC1 MTX LN_IC50")
    ax.set_ylabel("GDSC2 MTX LN_IC50")
    ax.set_title("Shared cell lines: n=8, ρ=0.90", loc="left")
    S.style_ax(ax); S.panel_label(ax, "C", x=-0.22)

    ax = fig.add_subplot(gs[1, 1])
    for i, screen in enumerate(["GDSC1", "GDSC2"]):
        s = loo[loo.screen.eq(screen)]
        ax.scatter(s.rho_omit, np.full(len(s), i), color=C["ball"] if i == 0 else C["tall"],
                   s=19, alpha=.7, edgecolor="#333", lw=.2)
        ax.plot([s.rho_omit.min(), s.rho_omit.max()], [i, i], color="#999", lw=.8)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["GDSC1", "GDSC2"])
    ax.set_xlim(.4, 1.05)
    ax.set_xlabel("ρ after omitting one cell line")
    ax.set_title("Leave-one-model sensitivity", loc="left")
    ax.text(.03, .5, "Screens overlap in 8 models; not independent validation",
            transform=ax.transAxes, fontsize=5.1)
    S.style_ax(ax); S.panel_label(ax, "D", x=-0.22)

    S.save_supp(fig, "FigureS11_MTX_PRISM")


# ===========================================================================
# S12  PDX L-IC + mouse + survival note  (M10, §3.6)
# ===========================================================================
def fig_s12():
    lic = pd.read_csv(_p("module10", "M10_10.2_LIC_vs_DP.tsv"), sep="\t")
    ko = pd.read_csv(_p("module10", "M10_r2_GSE287751_DN1_KO_vs_ctrl.tsv"), sep="\t")
    hto = pd.read_csv(_sd("validation/gse287751/author_condition_contrasts.tsv"), sep="\t")
    lib = pd.read_csv(_p("module10", "M10_r2_GSE260948_library_means.tsv"), sep="\t")
    mrd = pd.read_csv(_sd("SourceData_FigS12_MRD.tsv"), sep="\t")

    fig = plt.figure(figsize=(180 * MM, 210 * MM))
    gs = GridSpec(3, 2, figure=fig, hspace=0.65, wspace=0.55)

    ax = fig.add_subplot(gs[0, 0])
    keys = ["ZEB1", "LMO2", "CD1A"]
    values = pd.read_csv(_p("module10", "M10_10.1_keygene_values.tsv"), sep="\t")
    sub = values[values.kind.eq("fpkm") & values.gene.isin(keys) & values.group.isin(["LIC", "DP"])].copy()
    sub["log2_1pFPKM"] = np.log2(1 + sub.value)
    sub.to_csv(_sd("SourceData_FigS12A_library_expression.tsv"), sep="\t", index=False)
    pats = ["REL13", "REL14"]
    x = np.arange(len(keys))
    for i, pat in enumerate(pats):
        color = C["tall"] if i == 0 else C["ball"]
        for j, group in enumerate(["LIC", "DP"]):
            for k, gene in enumerate(keys):
                v = sub[(sub.patient == pat) & (sub.group == group) & (sub.gene == gene)].log2_1pFPKM.to_numpy()
                center = k + (i-.5)*.35 + (j-.5)*.12
                ax.scatter(center + np.linspace(-.025,.025,len(v)), v, s=14,
                           marker="o" if group == "LIC" else "^", color=color, alpha=.8,
                           label=f"{pat} {group}" if k==0 else None)
    ax.set_xticks(x)
    ax.set_xticklabels(keys, fontsize=6.2)
    ax.set_ylabel("log2(1+FPKM)")
    ax.set_title("PDX compartments: individual libraries", loc="left")
    ax.legend(fontsize=4.7, ncol=2, loc="upper right", frameon=False)
    ax.margins(y=.25)
    S.style_ax(ax)
    S.panel_label(ax, "A", x=-0.16)

    ax = fig.add_subplot(gs[0, 1])
    pick = ["Lmo2", "Zeb1", "Gata3", "Tcf7", "Il7r", "Spi1", "Mef2c"]
    k = ko[ko["gene"].isin(pick)].copy()
    k = k.set_index("gene").reindex(pick)
    y = np.arange(len(k))
    cols = [C["purple"] if g == "Lmo2" else (C["tall"] if g == "Zeb1" else C["grey"])
            for g in k.index]
    ax.hlines(y, 0, k["log2FC_KO_vs_ctrl"], color=cols, lw=1)
    ax.scatter(k["log2FC_KO_vs_ctrl"], y, color=cols, s=24, zorder=3)
    ax.axvline(0, color="#333", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(k.index, fontsize=6.2)
    ax.set_xlabel("log2FC (Lmo2 KO vs ctrl, DN1)")
    ax.set_title("Lmo2 KO: DN1 gene responses", loc="left")
    S.style_ax(ax)
    S.panel_label(ax, "B", x=-0.18)

    ax = fig.add_subplot(gs[1, 0])
    genes = ["Lmo2", "Zeb1", "Tcf7", "Gata3", "Bcl11b", "Il7r"]
    contexts = [("preNotch", "unprimed"), ("preNotch", "flt3l"),
                ("d3", "unprimed"), ("d3", "flt3l")]
    matrix = hto.pivot(index="gene", columns=["timepoint", "priming"],
                       values="delta_KO_minus_EV").loc[genes, contexts]
    limit = float(np.abs(matrix.to_numpy()).max())
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_yticks(range(len(genes)), genes, fontsize=6)
    ax.set_xticks(range(4), ["preN\nU", "preN\nF", "d3\nU", "d3\nF"], fontsize=5.5)
    ax.tick_params(length=0)
    for yi in range(len(genes)):
        for xi in range(4):
            value = matrix.iloc[yi, xi]
            ax.text(xi, yi, f"{value:+.2f}", ha="center", va="center", fontsize=5,
                    color="white" if abs(value) > limit * .55 else C["ink"])
    cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.03)
    cb.set_label("KO - EV: log2(1+CPM)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("Author-labelled scRNA perturbation", loc="left")
    S.panel_label(ax, "C", x=-0.16)

    ax = fig.add_subplot(gs[1, 1])
    cox = pd.read_csv(_sd("F7_target_cox.tsv"), sep="\t")
    z = cox[(cox.term == "ZEB1_z") & (cox.model.isin(["ZEB1_unadj", "ZEB1_age"]))].copy()
    z["lab"] = z["endpoint"] + np.where(z.model.eq("ZEB1_unadj"), " unadj.", " + age")
    # Subtype-adjusted models have zero-event strata and separation symptoms.
    order = ["OS unadj.", "OS + age", "EFS unadj.", "EFS + age"]
    z = z.set_index("lab").reindex(order).reset_index()
    y = np.arange(len(z))[::-1]
    ax.axvline(1, color="#333", lw=0.6)
    ax.hlines(y, z["ci_lo"], z["ci_hi"], color="#B8BCC2", lw=1.1)
    ax.scatter(z["hr"], y, s=26, color=C["grey"], edgecolor="#333", lw=0.35, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(z["lab"], fontsize=6.0)
    ax.set_xlabel("HR per SD ZEB1 (TARGET n=262)")
    ax.set_title("OS / EFS Cox, continuous ZEB1", loc="left")
    ax.set_xlim(0.3, 2.6)
    ax.text(0.97, 0.06, "OS 14 events; EFS 27\nsubtype model unstable",
            transform=ax.transAxes, ha="right", fontsize=5.5, color="#555")
    S.style_ax(ax)
    S.panel_label(ax, "D", x=-0.42)
    _ = mrd, lib

    paired = pd.read_csv(_sd("validation/gse144035/paired_model_changes.tsv"), sep="\t")
    for j, gene in enumerate(["Zeb1", "Lmo2"]):
        ax = fig.add_subplot(gs[2, j])
        for row in paired[paired.gene.eq(gene)].itertuples(index=False):
            col = C["tall"] if row.status == "Evolving" else C["ball"]
            ax.plot([0, 1], [row.minus_log2_1pCPM, row.plus_log2_1pCPM], 'o-', color=col, lw=.9, ms=3, alpha=.8)
        ax.set_xticks([0, 1], ["−Dox", "+Dox"])
        ax.set_xlim(-.15, 1.15)
        ax.set_ylabel(f"{gene}: log2(1+CPM)")
        ax.set_title(f"GSE144035: model-matched {gene}", loc="left")
        ax.text(.02, .02, "3 evolving / 3 independent models\nOne library per model × condition", transform=ax.transAxes, fontsize=5.2)
        S.style_ax(ax)
        S.panel_label(ax, ["E", "F"][j], x=-.22)

    S.save_supp(fig, "FigureS12_PDX_survival")


if __name__ == "__main__":
    for fn in (fig_s1, fig_s2, fig_s3, fig_s4, fig_s5, fig_s6,
               fig_s7, fig_s8, fig_s9, fig_s10, fig_s11, fig_s12):
        print("→", fn.__name__)
        fn()
    print("all supplementary figures written to figures/Supplementary/")
