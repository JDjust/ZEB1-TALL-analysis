# -*- coding: utf-8 -*-
"""
Compute the report.md §3.1–3.6 analyses from already-local public tables
and write Source Data used by the supplementary figures.

3.1 LMO2–ZEB1 false edge across TARGET + independent GEO cohorts
3.2 GSE206710 three-donor thymus atlas (Park/HTA fallback)
3.3 MTX / antifolate / HDAC / BCL2: pharmacotype proxies + DepMap target
    essentiality (power-limited; no local PRISM matrix)
3.4 TF-activity layer from Module 4 C3-TFT scores (decoupler-py fallback)
3.5 Two-cohort single-cell mini-meta (patient-level)
3.6 Survival: TARGET follow-up is not in the local bundle; MRD is the
    available clinical endpoint (already in Figure 7B)
"""
from __future__ import annotations

import os
import csv
import numpy as np
import pandas as pd
from scipy import stats
import _style as S


def residualize(y, groups):
    """OLS residuals of y on dummy-coded groups (drop-first)."""
    y = np.asarray(y, dtype=float)
    g = pd.Series(groups).astype(str).fillna("NA")
    dummies = pd.get_dummies(g, drop_first=True).astype(float)
    X = np.column_stack([np.ones(len(y)), dummies.to_numpy()])
    mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    resid = np.full(len(y), np.nan)
    if mask.sum() < (X.shape[1] + 3):
        return resid
    beta, *_ = np.linalg.lstsq(X[mask], y[mask], rcond=None)
    resid[mask] = y[mask] - X[mask] @ beta
    return resid


def spearman_ci(x, y, n_boot=2000, seed=1):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 6:
        return dict(n=n, rho=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        r, _ = stats.spearmanr(x[i], y[i])
        if np.isfinite(r):
            boots.append(r)
    lo, hi = np.percentile(boots, [2.5, 97.5]) if boots else (np.nan, np.nan)
    return dict(n=n, rho=float(rho), p=float(p), ci_lo=float(lo), ci_hi=float(hi))


def fisher_meta(rows):
    """Inverse-variance meta-analysis of Spearman ρ via Fisher's z."""
    zs, ws = [], []
    for r in rows:
        n, rho = r["n"], r["rho"]
        if n < 6 or not np.isfinite(rho):
            continue
        rho = float(np.clip(rho, -0.999, 0.999))
        z = 0.5 * np.log((1 + rho) / (1 - rho))
        se = 1.0 / np.sqrt(n - 3)
        zs.append(z)
        ws.append(1.0 / se ** 2)
    if not ws:
        return dict(n_cohorts=0, n_total=0, rho=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    w = np.array(ws)
    z = np.array(zs)
    zbar = np.sum(w * z) / np.sum(w)
    se = 1.0 / np.sqrt(np.sum(w))
    zstat = zbar / se
    p = 2 * stats.norm.sf(abs(zstat))
    rho = np.tanh(zbar)
    return dict(
        n_cohorts=len(ws),
        n_total=int(sum(r["n"] for r in rows if r["n"] >= 6 and np.isfinite(r["rho"]))),
        rho=float(rho),
        p=float(p),
        ci_lo=float(np.tanh(zbar - 1.96 * se)),
        ci_hi=float(np.tanh(zbar + 1.96 * se)),
    )


def write(df, name):
    path = os.path.join(S.DATA_DIR, name)
    df.to_csv(path, sep="\t", index=False)
    print("wrote", os.path.relpath(path, S.TOTAL_DIR), "nrow=", len(df))
    return path


# ---------------------------------------------------------------------------
# 3.1  LMO2–ZEB1 false edge
# ---------------------------------------------------------------------------
kg = pd.read_csv(S.tbl("module3", "independent_keygenes.tsv"), sep="\t")
kg = kg[kg["disease"].fillna("T-ALL").str.contains("T-ALL", case=False, na=True)].copy()
kg["level1"] = kg["level1"].replace({"LMO2/LYL1": "Immature", "NKX2_1": "NKX"})
keep_sub = {"TAL", "TLX", "HOXA", "Immature", "NKX"}
kg = kg[kg["level1"].isin(keep_sub)].copy()
kg = kg.dropna(subset=["ZEB1", "LMO2"])

edge_rows = []
net_rows = []
genes_net = ["ZEB1", "LMO2", "GATA3", "TCF7", "BCL11B", "LYL1"]

for cohort, d in kg.groupby("cohort", sort=False):
    raw = spearman_ci(d["ZEB1"], d["LMO2"], seed=11)
    zres = residualize(d["ZEB1"].to_numpy(), d["level1"])
    lres = residualize(d["LMO2"].to_numpy(), d["level1"])
    adj = spearman_ci(zres, lres, seed=12)
    n_sub = int(d["level1"].nunique())
    edge_rows.append({
        "cohort": cohort, "n": raw["n"], "n_subtype": n_sub,
        "rho_raw": raw["rho"], "p_raw": raw["p"],
        "ci_raw_lo": raw["ci_lo"], "ci_raw_hi": raw["ci_hi"],
        "rho_adj": adj["rho"], "p_adj": adj["p"],
        "ci_adj_lo": adj["ci_lo"], "ci_adj_hi": adj["ci_hi"],
        "delta_rho": adj["rho"] - raw["rho"] if np.isfinite(adj["rho"]) else np.nan,
    })
    # partial-correlation style network on residuals vs raw (TARGET highlighted)
    present = [g for g in genes_net if g in d.columns]
    if len(present) >= 3:
        mat = d[present].apply(pd.to_numeric, errors="coerce")
        res_mat = pd.DataFrame({g: residualize(mat[g].to_numpy(), d["level1"]) for g in present})
        for a, b in [("LMO2", "ZEB1"), ("GATA3", "ZEB1"), ("TCF7", "ZEB1"),
                     ("BCL11B", "ZEB1"), ("LYL1", "ZEB1"), ("LMO2", "LYL1")]:
            if a not in present or b not in present:
                continue
            r0 = spearman_ci(mat[a], mat[b], seed=13)
            r1 = spearman_ci(res_mat[a], res_mat[b], seed=14)
            net_rows.append({
                "cohort": cohort, "edge": f"{a}–{b}", "a": a, "b": b,
                "rho_raw": r0["rho"], "p_raw": r0["p"],
                "rho_adj": r1["rho"], "p_adj": r1["p"], "n": r0["n"],
            })

edge = pd.DataFrame(edge_rows)
net = pd.DataFrame(net_rows)
write(edge, "SourceData_FigS4_LMO2_ZEB1_false_edge.tsv")
write(net, "SourceData_FigS4_partial_network.tsv")
print("\n=== 3.1 LMO2–ZEB1 false edge ===")
print(edge.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


# ---------------------------------------------------------------------------
# 3.2  GSE206710 second normal thymus reference
# ---------------------------------------------------------------------------
st = pd.read_csv(S.tbl("module2", "M2_r2_GSE206710_stage_gene_means.tsv"), sep="\t")
sm = pd.read_csv(S.tbl("module2", "M2_r2_GSE206710_sample_means.tsv"), sep="\t")
don = pd.read_csv(S.tbl("module2", "M2_r2_GSE206710_donor_fraction.tsv"), sep="\t")
order = {"DN_early": 1, "DP_immature": 2, "DP_mature": 3}
st = st.copy()
st["stage_rank"] = st["inferred_stage"].map(order)
rho_stage, p_stage = stats.spearmanr(st["stage_rank"], st["g_ZEB1"])
# paired CD3− vs CD3+ by donor
wide = don.pivot(index="donor", columns="fraction", values="ZEB1")
cd3_rows = []
for donor, row in wide.iterrows():
    cd3_rows.append({
        "donor": donor,
        "ZEB1_CD3neg": float(row.get("CD3neg", np.nan)),
        "ZEB1_CD3pos": float(row.get("CD3pos", np.nan)),
        "delta_neg_minus_pos": float(row.get("CD3neg", np.nan) - row.get("CD3pos", np.nan)),
    })
cd3 = pd.DataFrame(cd3_rows)
if len(cd3) >= 3:
    wstat = stats.wilcoxon(cd3["ZEB1_CD3neg"], cd3["ZEB1_CD3pos"], alternative="two-sided")
    cd3_p = float(wstat.pvalue)
else:
    cd3_p = np.nan
# sample-level ZEB1 vs T-diff / stage
sm = sm.copy()
sm["stage_rank"] = sm["stage"].map(order)
rho_tdiff, p_tdiff = stats.spearmanr(sm["g_ZEB1"], sm["sig_Tdiff"], nan_policy="omit")
sum2 = pd.DataFrame([{
    "atlas": "GSE206710",
    "n_donors": int(don["donor"].nunique()),
    "n_libraries": int(sm["sample"].nunique()) if "sample" in sm else int(len(sm)),
    "n_cells_DN": int(st.loc[st.inferred_stage == "DN_early", "g_ZEB1"].notna().sum()),
    "ZEB1_vs_stage_rho": float(rho_stage),
    "ZEB1_vs_stage_p": float(p_stage),
    "ZEB1_vs_Tdiff_rho": float(rho_tdiff),
    "ZEB1_vs_Tdiff_p": float(p_tdiff),
    "CD3_paired_wilcoxon_p": cd3_p,
    "mean_delta_CD3neg_minus_pos": float(cd3["delta_neg_minus_pos"].mean()),
    "note": "Second normal reference (3 donors). ZEB1 remains detectable in CD3+; no monotonic rise with stage.",
}])
write(st, "SourceData_FigS2_stage_means.tsv")
write(sm, "SourceData_FigS2_library_means.tsv")
write(cd3, "SourceData_FigS2_CD3_paired.tsv")
write(sum2, "SourceData_FigS2_summary.tsv")
print("\n=== 3.2 GSE206710 ===")
print(sum2.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


# ---------------------------------------------------------------------------
# 3.3  MTX / antifolate / HDAC / BCL2  (pharmacotype + DepMap targets)
# ---------------------------------------------------------------------------
drug = pd.read_csv(S.tbl("module7", "M7_r2_ZEB1_rho_with_class.tsv"), sep="\t")
class_focus = {
    "antimetabolite": "antifolate/antimetabolite proxy (MTX absent)",
    "HDAC": "HDAC inhibitor",
    "BCL2": "BCL2 inhibitor",
}
d3 = drug.copy()
d3["focus"] = d3["class"].map(class_focus).fillna("other")
write(d3, "SourceData_FigS11_pharmacotype_classes.tsv")

# DepMap: ZEB1 RNA vs CRISPR gene-effect of MTX/HDAC/BCL2 targets in T-ALL
dep_note = "PRISM / GDSC MTX matrices were not in the local bundle."
dep_rows = []
try:
    model = pd.read_csv(os.path.join(S.ROOT, "data", "Model.csv"), usecols=[
        "ModelID", "CellLineName", "OncotreePrimaryDisease", "OncotreeSubtype",
        "OncotreeLineage", "DepmapModelType",
    ])
    is_tall = (
        model["OncotreePrimaryDisease"].fillna("").str.contains("T-Lymphoblastic", case=False)
        | model["OncotreeSubtype"].fillna("").str.contains("T-ALL|T-Lymphoblastic", case=False)
        | model["DepmapModelType"].fillna("").str.contains("TALL|T-ALL", case=False)
    )
    tall_ids = set(model.loc[is_tall, "ModelID"])
    # gene-effect header
    ge_path = os.path.join(S.ROOT, "data", "CRISPRGeneEffect.csv")
    expr_path = os.path.join(S.ROOT, "data", "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    want = ["ZEB1", "DHFR", "TYMS", "GART", "ATIC", "FPGS",
            "HDAC1", "HDAC2", "HDAC3", "BCL2", "BCL2L1", "MCL1"]

    def gene_col_index(path, genes):
        """Map column index -> gene; first column is ModelID even if unnamed."""
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            header = next(csv.reader(fh))
        idx = {0: "ModelID"}
        for i, h in enumerate(header):
            if i == 0:
                continue
            g = h.split(" ")[0].split("(")[0]
            if g in genes:
                idx[i] = g
        return idx

    if os.path.exists(ge_path) and os.path.exists(expr_path) and tall_ids:
        ge_idx = gene_col_index(ge_path, want)
        ex_idx = gene_col_index(expr_path, ["ZEB1"])
        ge = pd.read_csv(ge_path, usecols=sorted(ge_idx))
        ge.columns = [ge_idx[i] for i in sorted(ge_idx)]
        ge = ge[ge["ModelID"].isin(tall_ids)].copy()
        ge = ge.rename(columns={g: g + "_GE" for g in want if g in ge.columns})
        ex = pd.read_csv(expr_path, usecols=sorted(ex_idx))
        ex.columns = [ex_idx[i] for i in sorted(ex_idx)]
        ex = ex[ex["ModelID"].isin(tall_ids)].copy()
        ex = ex.rename(columns={"ZEB1": "ZEB1_RNA"})
        m = ge.merge(ex, on="ModelID", how="inner")
        name_map = model.set_index("ModelID")["CellLineName"]
        m["CellLineName"] = m["ModelID"].map(name_map)
        for gene in want:
            ge_col = gene + "_GE"
            if ge_col not in m.columns:
                continue
            sub = m[["ZEB1_RNA", ge_col, "CellLineName"]].dropna()
            sp = spearman_ci(sub["ZEB1_RNA"], sub[ge_col], seed=21)
            role = {
                "DHFR": "MTX primary target",
                "TYMS": "antifolate / 5-FU axis",
                "GART": "purine/antifolate",
                "ATIC": "purine/antifolate",
                "FPGS": "folate polyglutamation",
                "HDAC1": "HDAC inhibitor target",
                "HDAC2": "HDAC inhibitor target",
                "HDAC3": "HDAC inhibitor target",
                "BCL2": "venetoclax target",
                "BCL2L1": "BCL-XL",
                "MCL1": "MCL1",
                "ZEB1": "self-dependency",
            }.get(gene, "")
            dep_rows.append({
                "gene": gene, "role": role, "n": sp["n"],
                "rho_ZEB1RNA_vs_GE": sp["rho"], "p": sp["p"],
                "ci_lo": sp["ci_lo"], "ci_hi": sp["ci_hi"],
                "mean_GE": float(sub[ge_col].mean()) if len(sub) else np.nan,
            })
        write(m, "SourceData_FigS11_DepMap_TALL_targets.tsv")
        n_overlap = int(m["ZEB1_RNA"].notna().sum())
        dep_note = (
            f"Local DepMap CRISPR∩RNA T-ALL n={n_overlap}. "
            "PRISM / GDSC MTX dose-response matrices were not present locally; "
            "MTX is tested via DHFR essentiality (power-limited)."
        )
except Exception as exc:
    dep_note = f"DepMap target analysis skipped: {exc}"
    print("WARN 3.3 DepMap:", exc)

dep_sum = pd.DataFrame(dep_rows) if dep_rows else pd.DataFrame(
    columns=["gene", "role", "n", "rho_ZEB1RNA_vs_GE", "p", "ci_lo", "ci_hi", "mean_GE"]
)
if len(dep_sum):
    write(dep_sum, "SourceData_FigS11_DepMap_ZEB1_vs_targets.tsv")
note11 = pd.DataFrame([{"note": dep_note, "MTX_in_pharmacotype": "absent",
                        "TALL_n_CRISPR_RNA": int(dep_sum["n"].max()) if len(dep_sum) else 0}])
write(note11, "SourceData_FigS11_methods_note.tsv")
print("\n=== 3.3 MTX / DepMap targets ===")
print(dep_note)
if len(dep_sum):
    print(dep_sum.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


# ---------------------------------------------------------------------------
# 3.4  TF activity (C3-TFT / Module 4 substitute for decoupleR)
# ---------------------------------------------------------------------------
tf = pd.read_csv(S.tbl("module4", "M4_4.20_tf_activity_vs_ZEB1.tsv"), sep="\t")
kd = pd.read_csv(S.tbl("module9", "M9_keygene_effects.tsv"), sep="\t")
write(tf, "SourceData_FigS10_TF_activity.tsv")
# signed concordance: TLX1 KD raises ZEB1; patient TLX is ZEB1-high (opposite)
kd_on = kd[kd["experiment"] == "gse110635_combined"].copy()
write(kd_on, "SourceData_FigS10_TLX1KD_keygenes.tsv")
print("\n=== 3.4 TF activity / TLX1 KD ===")
print(tf.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


# ---------------------------------------------------------------------------
# 3.5  Two-cohort patient-level mini-meta
# ---------------------------------------------------------------------------
disc = pd.read_csv(S.tbl("module5", "M5_5.41_patient_ZEB1_associations.tsv"), sep="\t")
vali = pd.read_csv(S.tbl("module5", "M5_5.42_ZEB1_patient_assoc.tsv"), sep="\t")
alias = {
    "g_ZEB2": "ZEB2", "g_LMO2": "LMO2", "g_GATA3": "GATA3", "g_TCF7": "TCF7",
    "g_IL7R": "IL7R", "g_CD34": "CD34", "g_BCL11B": "BCL11B", "g_LYL1": "LYL1",
    "sig_ETP": "ETP signature", "sig_Stemness": "Stemness",
    "sig_Tcell_diff": "T-cell differentiation", "sig_Tdiff": "T-cell differentiation",
    "sig_GATA3prog": "GATA3 program", "sig_LMO2prog": "LMO2 program",
    "sig_TGFb": "TGF-β", "sig_Prolif": "Proliferation",
}
disc = disc.copy()
vali = vali.copy()
disc["feature_std"] = disc["feature"].map(alias).fillna(disc["feature"])
vali["feature_std"] = vali["feature"].map(alias).fillna(vali["feature"])
shared = sorted(set(disc["feature_std"]) & set(vali["feature_std"]))
meta_rows = []
for feat in shared:
    a = disc.loc[disc.feature_std == feat].iloc[0]
    b = vali.loc[vali.feature_std == feat].iloc[0]
    rA = dict(n=int(a["n"]) if pd.notna(a["n"]) else 0,
              rho=float(a["rho"]) if pd.notna(a["rho"]) else np.nan)
    rB = dict(n=int(b["n"]) if pd.notna(b["n"]) else 0,
              rho=float(b["rho"]) if pd.notna(b["rho"]) else np.nan)
    mm = fisher_meta([rA, rB])
    same = np.sign(rA["rho"]) == np.sign(rB["rho"]) if np.isfinite(rA["rho"]) and np.isfinite(rB["rho"]) else False
    meta_rows.append({
        "feature": feat,
        "rho_GSE227122": rA["rho"], "p_GSE227122": float(a["p"]) if pd.notna(a["p"]) else np.nan, "n_GSE227122": rA["n"],
        "rho_GSE248287": rB["rho"], "p_GSE248287": float(b["p"]) if pd.notna(b["p"]) else np.nan, "n_GSE248287": rB["n"],
        "same_sign": bool(same),
        "meta_rho": mm["rho"], "meta_p": mm["p"],
        "meta_ci_lo": mm["ci_lo"], "meta_ci_hi": mm["ci_hi"],
        "n_total": mm["n_total"],
    })
mini = pd.DataFrame(meta_rows).sort_values("meta_rho")
write(mini, "SourceData_FigS6_mini_meta.tsv")
print("\n=== 3.5 scRNA mini-meta ===")
print(mini[["feature", "rho_GSE227122", "rho_GSE248287", "same_sign", "meta_rho", "meta_p"]].to_string(
    index=False, float_format=lambda v: f"{v:.3g}"))


# ---------------------------------------------------------------------------
# 3.6  Survival note (TARGET follow-up not local)
# ---------------------------------------------------------------------------
mrd = pd.read_csv(S.tbl("module7", "M7_7.18_ZEB1_MRD.tsv"), sep="\t")
surv = pd.DataFrame([{
    "analysis": "TARGET T-ALL Cox OS/EFS",
    "status": "not computed",
    "reason": (
        "TARGET-ALL-P2 clinical follow-up Excel files live on the analysis "
        "server (Clinical_Supplement) and are not in the local total/ bundle. "
        "No OS/EFS table is available to residualize ZEB1 on subtype/age/MRD."
    ),
    "available_clinical_endpoint": "Day-15 and Day-42 MRD (Figure 7B / Module 7)",
    "MRD15_rho": float(mrd.loc[mrd.index[0], "rho"]) if "rho" in mrd.columns else np.nan,
    "interpretation": (
        "ZEB1 is not claimed as an independent prognostic factor. "
        "The available clinical endpoint (MRD) is weak and unstable."
    ),
}])
write(surv, "SourceData_FigS12_survival_note.tsv")
write(mrd, "SourceData_FigS12_MRD.tsv")
print("\n=== 3.6 survival ===")
print(surv["status"].iloc[0], "|", surv["available_clinical_endpoint"].iloc[0])
print("done.")
