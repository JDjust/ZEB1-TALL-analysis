#!/usr/bin/env python3
"""Module 8 redo on DepMap Public 26Q1 (expression + CRISPR + CNV + dependency).

ZEB1 stays continuous. n is cell lines. No MTX story.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1")
EXPR_FILE = "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CNV_FILE = "OmicsCNGeneWGS.csv"
DEP_FILE = "CRISPRGeneDependency.csv"
M7 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module7/tables")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8")
PROC, TAB = ROOT / "processed", ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

DRUG_TARGETS = {
    "Venetoclax": ["BCL2", "BCL2L1", "MCL1"],
    "Ruxolitinib": ["JAK1", "JAK2", "JAK3", "TYK2"],
    "Dasatinib": ["ABL1", "SRC", "LCK", "FYN", "YES1"],
    "Panobinostat": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
    "Vorinostat": ["HDAC1", "HDAC2", "HDAC3"],
    "Trametinib": ["MAP2K1", "MAP2K2"],
    "Bortezomib": ["PSMB5", "PSMB1"],
    "Ibrutinib": ["BTK", "BLK"],
    "Nelarabine": ["DCK", "NT5C2"],
}


def log(m):
    print(m, flush=True)


def hedges_g(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 4 or n2 < 4:
        return np.nan
    sp = math.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return 0.0
    d = (a.mean() - b.mean()) / sp
    j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    return j * d


def fdr_bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    mask = np.isfinite(p)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = np.empty(n, dtype=float)
    prev = 1.0
    for i, k in enumerate(order[::-1]):
        rank = n - i
        val = min(prev, pv[k] * n / rank)
        ranked[k] = val
        prev = val
    out[mask] = ranked
    return out


def gene_token(col):
    return re.split(r"[\s(]", str(col), maxsplit=1)[0]


def pick_gene_col(cols, gene):
    for c in cols:
        if gene_token(c) == gene:
            return c
    return None


def load_models():
    m = pd.read_csv(DATA / "Model.csv")
    m["lineage"] = np.where(
        m["OncotreePrimaryDisease"].eq("T-Lymphoblastic Leukemia/Lymphoma"), "T-ALL",
        np.where(m["OncotreePrimaryDisease"].eq("B-Lymphoblastic Leukemia/Lymphoma"), "B-ALL",
                 np.where(m["OncotreeLineage"].eq("Lymphoid"), "other_lymphoid",
                          np.where(m["OncotreeLineage"].eq("Myeloid"), "myeloid", "other"))),
    )
    m["is_tall"] = m["lineage"].eq("T-ALL")
    log("models " + str(m.lineage.value_counts().to_dict()))
    m.to_csv(PROC / "depmap_models.tsv", sep="\t", index=False)
    return m


def load_subset_csv(path, keep_ids):
    keep = set(keep_ids)
    chunks = []
    for ch in pd.read_csv(path, chunksize=250):
        if "ModelID" in ch.columns:
            key = "ModelID"
        else:
            key = ch.columns[0]
        sub = ch[ch[key].astype(str).isin(keep)].copy()
        if "IsDefaultEntryForModel" in sub.columns:
            flag = sub["IsDefaultEntryForModel"]
            if flag.isin([True, "True", "true", 1, "1"]).any():
                sub = sub[flag.isin([True, "True", "true", 1, "1"])]
        chunks.append(sub)
    if not chunks:
        return pd.DataFrame(columns=["ModelID"])
    out = pd.concat(chunks, ignore_index=True)
    if "ModelID" not in out.columns:
        out = out.rename(columns={str(out.columns[0]): "ModelID"})
    drop = [c for c in ["SequencingID", "ModelConditionID", "IsDefaultEntryForMC",
                        "IsDefaultEntryForModel"] if c in out.columns]
    out = out.drop(columns=drop, errors="ignore")
    if out.columns.duplicated().any():
        out = out.loc[:, ~out.columns.duplicated()].copy()
    out = out.drop_duplicates("ModelID")
    return out


def main():
    models = load_models()
    keep = models.loc[models["lineage"].isin(["T-ALL", "B-ALL", "myeloid", "other_lymphoid"]), "ModelID"]
    log(f"load expression for {len(keep)} hematopoietic models")
    expr_path = DATA / EXPR_FILE
    if not expr_path.exists():
        alt = DATA / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
        expr_path = alt if alt.exists() else expr_path
    log(f"expression file {expr_path.name} exists={expr_path.exists()}")
    expr = load_subset_csv(expr_path, keep)
    log(f"expr {expr.shape}")
    extra = [
        "ZEB1", "ZEB2", "LMO2", "JAK1", "JAK2", "BCL2", "MCL1", "CDK6", "NOTCH1",
        "GATA3", "BCL11B", "IL7R", "TCF7", "MYC", "LCK", "ABL1", "HDAC1",
    ]
    ren = {}
    for g in extra:
        c = pick_gene_col(expr.columns, g)
        if c and c not in ren:
            ren[c] = g
    if "ZEB1" not in ren.values():
        raise SystemExit("ZEB1 not in expression matrix")
    expr = expr.rename(columns=ren)
    keep_expr = ["ModelID"] + [g for g in extra if g in expr.columns]
    lab = models.merge(expr[keep_expr], on="ModelID", how="left")
    lab.to_csv(TAB / "M8_8.1_8.3_models_expression.tsv", sep="\t", index=False)
    log(f"models-with-or-without-expr {lab.lineage.value_counts().to_dict()}")
    log(f"T-ALL with ZEB1 expr n={lab.loc[lab.lineage.eq('T-ALL'), 'ZEB1'].notna().sum()}")

    log("load CRISPR GeneEffect")
    ge = load_subset_csv(DATA / "CRISPRGeneEffect.csv", keep)
    log(f"crispr {ge.shape}")
    gene_cols = [c for c in ge.columns if c != "ModelID"]
    ge_z = pick_gene_col(gene_cols, "ZEB1")
    ge_z2 = pick_gene_col(gene_cols, "ZEB2")
    ge_l = pick_gene_col(gene_cols, "LMO2")
    slim = ge[["ModelID"]].copy()
    if ge_z:
        slim["ZEB1_GE"] = ge[ge_z]
    if ge_z2:
        slim["ZEB2_GE"] = ge[ge_z2]
    if ge_l:
        slim["LMO2_GE"] = ge[ge_l]
    lab = lab.merge(slim, on="ModelID", how="left")
    lab.to_csv(TAB / "M8_8.1_8.5_line_table.tsv", sep="\t", index=False)

    tall_ids = set(lab.loc[lab.lineage.eq("T-ALL"), "ModelID"])
    ball_ids = set(lab.loc[lab.lineage.eq("B-ALL"), "ModelID"])
    ge_t = ge[ge.ModelID.isin(tall_ids)].set_index("ModelID")
    ge_b = ge[ge.ModelID.isin(ball_ids)].set_index("ModelID")
    zeb1_map = lab.set_index("ModelID")["ZEB1"]

    # 8.4-8.5 self dependency
    self = lab.loc[lab.lineage.eq("T-ALL"), ["ModelID", "CellLineName", "ZEB1", "ZEB1_GE"]].copy()
    self["dependent_0.5"] = self["ZEB1_GE"] < -0.5
    self.to_csv(TAB / "M8_8.4_8.5_ZEB1_self.tsv", sep="\t", index=False)
    n_dep = int(self["dependent_0.5"].sum())
    log(f"T-ALL ZEB1 GeneEffect mean={self.ZEB1_GE.mean():.3f} n_dep={n_dep}/{len(self)}")

    # 8.6/8.8 genome-wide Spearman GeneEffect vs ZEB1 expression in T-ALL
    zmap = pd.to_numeric(zeb1_map, errors="coerce")
    zmap = zmap.groupby(level=0).first()
    common = [i for i in ge_t.index if i in zmap.index and np.isfinite(zmap.loc[i])]
    log(f"T-ALL CRISPR n={len(ge_t)} with finite ZEB1 expr n={len(common)}")
    rows = []
    min_n = 6
    if len(common) >= min_n:
        x = zmap.loc[common].astype(float)
        mat = ge_t.loc[common, gene_cols]
        xv = x.to_numpy().reshape(-1)
        for c in gene_cols:
            y = pd.to_numeric(mat[c], errors="coerce").to_numpy().reshape(-1)
            if y.size != xv.size:
                continue
            ok = np.isfinite(xv) & np.isfinite(y)
            if ok.sum() < min_n:
                continue
            rho, p = spearmanr(xv[ok], y[ok])
            if not np.isfinite(rho):
                continue
            rows.append({"gene": gene_token(c), "col": c, "rho": float(rho), "p": float(p), "n": int(ok.sum())})
    gw = pd.DataFrame(rows, columns=["gene", "col", "rho", "p", "n"])
    if len(gw):
        gw["fdr"] = fdr_bh(gw["p"].to_numpy())
        gw = gw.sort_values("p")
    else:
        gw["fdr"] = pd.Series(dtype=float)
        log("genome-wide empty: too few T-ALL lines with CRISPR+ZEB1")
    gw.to_csv(TAB / "M8_8.6_8.10_GE_vs_ZEB1.tsv", sep="\t", index=False)
    log(f"genome-wide n_genes={len(gw)} FDR<0.05={int((gw.fdr<0.05).sum()) if len(gw) else 0}")

    # 8.9/8.11/8.15 mean essentiality T-ALL vs B-ALL
    tmean = ge_t[gene_cols].mean(axis=0)
    bmean = ge_b[gene_cols].mean(axis=0)
    ess = pd.DataFrame({
        "gene": [gene_token(c) for c in gene_cols],
        "col": gene_cols,
        "mean_TALL": tmean.to_numpy(),
        "mean_BALL": bmean.reindex(gene_cols).to_numpy(),
        "n_TALL": len(ge_t),
        "n_BALL": len(ge_b),
    })
    mye_ids = set(lab.loc[lab.lineage.eq("myeloid"), "ModelID"])
    ge_m = ge[ge.ModelID.isin(mye_ids)].set_index("ModelID")
    ess["mean_myeloid"] = ge_m[gene_cols].mean(axis=0).reindex(gene_cols).to_numpy() if len(ge_m) else np.nan
    ess["delta_T_minus_B"] = ess["mean_TALL"] - ess["mean_BALL"]
    ess["delta_T_minus_M"] = ess["mean_TALL"] - ess["mean_myeloid"]
    # Wilcoxon on a subset of genes that are strong in either lineage to keep runtime OK
    strong = ess[(ess.mean_TALL < -0.35) | (ess.mean_BALL < -0.35) | (ess.delta_T_minus_B.abs() > 0.25)]
    pvals = []
    gs = []
    for _, r in strong.iterrows():
        a = pd.to_numeric(ge_t[r["col"]], errors="coerce").to_numpy()
        b = pd.to_numeric(ge_b[r["col"]], errors="coerce").to_numpy()
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) < 6 or len(b) < 6:
            pvals.append(np.nan); gs.append(np.nan); continue
        try:
            pvals.append(mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:
            pvals.append(np.nan)
        gs.append(hedges_g(a, b))
    strong = strong.copy()
    strong["p_T_vs_B"] = pvals
    strong["g_T_vs_B"] = gs
    strong["fdr_T_vs_B"] = fdr_bh(np.array(pvals, dtype=float))
    strong.sort_values("mean_TALL").to_csv(TAB / "M8_8.9_8.15_essentiality.tsv", sep="\t", index=False)

    # 8.7 ZEB1 Q1 vs Q4 GeneEffect (n is small; still report)
    tall = lab.loc[lab.lineage.eq("T-ALL")].dropna(subset=["ZEB1"])
    q1, q3 = tall["ZEB1"].quantile(0.25), tall["ZEB1"].quantile(0.75)
    low_ids = set(tall.loc[tall.ZEB1 <= q1, "ModelID"])
    high_ids = set(tall.loc[tall.ZEB1 >= q3, "ModelID"])
    ge_low = ge[ge.ModelID.isin(low_ids)].set_index("ModelID")
    ge_high = ge[ge.ModelID.isin(high_ids)].set_index("ModelID")
    log(f"8.7 Q1 n={len(ge_low)} Q4 n={len(ge_high)}")
    rows7 = []
    test_cols = list(dict.fromkeys(list(strong["col"]) + [c for c in gene_cols if gene_token(c) in extra]))
    for c in test_cols:
        a = pd.to_numeric(ge_low[c], errors="coerce").to_numpy()
        b = pd.to_numeric(ge_high[c], errors="coerce").to_numpy()
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) < 3 or len(b) < 3:
            continue
        try:
            p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        except Exception:
            continue
        rows7.append({
            "gene": gene_token(c), "mean_Q1": float(a.mean()), "mean_Q4": float(b.mean()),
            "delta_Q4_minus_Q1": float(b.mean() - a.mean()), "g": hedges_g(b, a),
            "p": p, "n_Q1": len(a), "n_Q4": len(b),
        })
    qtab = pd.DataFrame(rows7)
    if len(qtab):
        qtab["fdr"] = fdr_bh(qtab["p"].to_numpy())
        qtab.sort_values("p").to_csv(TAB / "M8_8.7_Q1Q4_GE.tsv", sep="\t", index=False)
        log(f"8.7 genes={len(qtab)} FDR<0.05={(qtab.fdr < 0.05).sum()}")

    # 8.17 expression-adjusted GE vs ZEB1 for key genes
    tlab = lab.loc[lab.lineage.eq("T-ALL")].set_index("ModelID")
    rows17 = []
    for g in extra:
        ge_col = pick_gene_col(gene_cols, g)
        if ge_col is None or g not in tlab.columns:
            continue
        y = pd.to_numeric(ge_t.reindex(tlab.index)[ge_col], errors="coerce")
        xg = pd.to_numeric(tlab[g], errors="coerce")
        z = pd.to_numeric(tlab["ZEB1"], errors="coerce")
        ok = np.isfinite(y.to_numpy()) & np.isfinite(xg.to_numpy()) & np.isfinite(z.to_numpy())
        if ok.sum() < 8:
            continue
        yy, xx, zz = y.to_numpy()[ok], xg.to_numpy()[ok], z.to_numpy()[ok]
        X = np.column_stack([np.ones(len(xx)), xx])
        coef, *_ = np.linalg.lstsq(X, yy, rcond=None)
        resid = yy - X @ coef
        rho_adj, p_adj = spearmanr(zz, resid)
        rho_raw, p_raw = spearmanr(zz, yy)
        rows17.append({
            "gene": g, "n": int(ok.sum()), "rho_raw": rho_raw, "p_raw": p_raw,
            "rho_expr_adjusted": rho_adj, "p_expr_adjusted": p_adj,
            "lm_intercept": float(coef[0]), "lm_expr": float(coef[1]),
        })
    pd.DataFrame(rows17).to_csv(TAB / "M8_8.17_expr_adjusted.tsv", sep="\t", index=False)

    # 8.16 pan-cancer: ZEB1 GeneEffect by lineage (all models that have GE)
    ge_all_z = None
    if ge_z:
        tmp = ge[["ModelID", ge_z]].rename(columns={ge_z: "ZEB1_GE"})
        pan = models.merge(tmp, on="ModelID", how="inner")
        pan.groupby("OncotreeLineage", dropna=False)["ZEB1_GE"].agg(["mean", "median", "count"]).reset_index().to_csv(
            TAB / "M8_8.16_ZEB1_GE_by_lineage.tsv", sep="\t", index=False
        )
        ge_all_z = pan

    # 8.17 expression-adjusted: residual GE after target expression, vs ZEB1 (top genes only)
    # skip full matrix; do for ZEB1/JAK1/BCL2/CDK6/NOTCH1 if present
    # 8.19 drug-target join
    targets = sorted({g for vs in DRUG_TARGETS.values() for g in vs})
    hit = ess[ess.gene.isin(targets + ["ZEB1", "ZEB2", "LMO2", "CDK6", "NOTCH1", "IL7R", "GATA3", "BCL11B"])]
    hit.to_csv(TAB / "M8_8.19_drug_target_GE.tsv", sep="\t", index=False)
    ess[ess.gene.isin(extra + targets)].to_csv(TAB / "M8_8.14_key_lineage_GE.tsv", sep="\t", index=False)

    m7 = M7 / "M7_7.7_7.8_drug_spearman.tsv"
    pharm = pd.read_csv(m7, sep="\t") if m7.exists() else pd.DataFrame()
    ev = ess.merge(gw[["gene", "rho", "p", "fdr"]], on="gene", how="left", suffixes=("", "_zeb1"))
    ev["ess_TALL"] = ev["mean_TALL"] < -0.5
    ev["tall_specific"] = (ev["mean_TALL"] < -0.35) & (ev["delta_T_minus_B"] < -0.2)
    ev["zeb1_assoc"] = (ev["fdr"] < 0.1) if "fdr" in ev.columns else False
    ev["score"] = ev["ess_TALL"].astype(int) + ev["tall_specific"].astype(int) + ev["zeb1_assoc"].fillna(False).astype(int)
    ev.sort_values(["score", "mean_TALL"]).tail(40)  # not used
    ev = ev.sort_values(["score", "mean_TALL"])
    ev.to_csv(TAB / "M8_8.22_8.23_evidence.tsv", sep="\t", index=False)
    cand = ev[ev.score >= 2].sort_values("mean_TALL")
    cand.to_csv(TAB / "M8_8.22_candidates.tsv", sep="\t", index=False)
    log(f"candidates score>=2: {len(cand)}")
    if len(pharm):
        pharm.to_csv(TAB / "M8_8.21_pharmacotype_spearman.tsv", sep="\t", index=False)

    # 8.13 co-dependency among top T-ALL essentials
    top = ess.nsmallest(40, "mean_TALL")["col"].tolist()
    sub = ge_t[top].apply(pd.to_numeric, errors="coerce")
    cm = sub.corr(method="spearman")
    cm.index = [gene_token(c) for c in cm.index]
    cm.columns = [gene_token(c) for c in cm.columns]
    cm.to_csv(TAB / "M8_8.13_codependency_corr.tsv", sep="\t")

    # 26Q1 extras: dependency probability and ZEB1 CNV
    dep_path = DATA / DEP_FILE
    if dep_path.exists() and ge_z:
        log("load CRISPRGeneDependency")
        dep = load_subset_csv(dep_path, keep)
        dcol = pick_gene_col([c for c in dep.columns if c != "ModelID"], "ZEB1")
        if dcol:
            dslim = dep[["ModelID", dcol]].rename(columns={dcol: "ZEB1_dependency"})
            lab = lab.merge(dslim, on="ModelID", how="left")
            lab.to_csv(TAB / "M8_8.1_8.5_line_table.tsv", sep="\t", index=False)
            tall_dep = lab.loc[lab.lineage.eq("T-ALL"), ["ModelID", "CellLineName", "ZEB1", "ZEB1_GE", "ZEB1_dependency"]]
            tall_dep.to_csv(TAB / "M8_r2_8.5_ZEB1_dependency.tsv", sep="\t", index=False)
            log(f"T-ALL with dependency n={tall_dep.ZEB1_dependency.notna().sum()} mean={tall_dep.ZEB1_dependency.mean():.3f}")
    cnv_path = DATA / CNV_FILE
    if cnv_path.exists():
        log("load CNV")
        cnv = load_subset_csv(cnv_path, keep)
        ccol = pick_gene_col([c for c in cnv.columns if c != "ModelID"], "ZEB1")
        if ccol:
            cslim = cnv[["ModelID", ccol]].rename(columns={ccol: "ZEB1_CNV"})
            lab = lab.merge(cslim, on="ModelID", how="left")
            lab.to_csv(TAB / "M8_8.1_8.5_line_table.tsv", sep="\t", index=False)
            lab.loc[lab.lineage.eq("T-ALL"), ["ModelID", "CellLineName", "ZEB1", "ZEB1_GE", "ZEB1_CNV"]].to_csv(
                TAB / "M8_r2_8.18_ZEB1_CNV.tsv", sep="\t", index=False
            )
            log(f"T-ALL with CNV n={lab.loc[lab.lineage.eq('T-ALL'), 'ZEB1_CNV'].notna().sum()}")
    n_t = int(((lab.lineage.eq("T-ALL")) & lab.get("ZEB1", pd.Series(dtype=float)).notna() & lab.get("ZEB1_GE", pd.Series(dtype=float)).notna()).sum()) if "ZEB1_GE" in lab.columns else 0
    Path(TAB / "M8_r2_n_note.txt").write_text(
        f"DepMap 26Q1 T-ALL with ZEB1 RNA and CRISPR n={n_t}\n"
        f"T-ALL models={(lab.lineage.eq('T-ALL')).sum()} B-ALL={(lab.lineage.eq('B-ALL')).sum()}\n"
    )
    lab.to_csv(PROC / "line_table_26q1.tsv", sep="\t", index=False)
    log("MODULE8 26Q1 PREPARE DONE")


if __name__ == "__main__":
    main()
