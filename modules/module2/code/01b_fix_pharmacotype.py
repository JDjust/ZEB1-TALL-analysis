#!/usr/bin/env python3
"""Pharmacotype maturity panel for Module 2.18.

FPKM columns are St. Jude Sample IDs (SJALL*), not Patient IDs (ALL*).
Reuse Module 1's already-matched T-ALL labels.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
M1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
M3 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed")
PROC = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/processed")

EARLY = ["CD34", "KIT", "IL7R", "SPI1", "LMO2", "LYL1", "MEF2C", "CD44", "HES1"]
LATE = ["CD1A", "CD1B", "CD3E", "CD3D", "CD4", "CD8A", "RAG1", "BCL11B", "TCF7", "LEF1", "DNTT"]
PANEL = sorted(set(EARLY + LATE + ["ZEB1", "ZEB2"]))


def nid(s) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def log(m):
    print(m, flush=True)


def nmatch(a, b) -> int:
    return len(set(map(nid, a)) & set(map(nid, b)))


def main():
    kg = pd.read_csv(M1 / "pharmacotype_keygenes.tsv", sep="\t")
    kg_t = kg[kg["disease"].eq("T-ALL")].copy()
    lab = pd.read_csv(M3 / "pharmacotype_tall_subtype.tsv", sep="\t")
    log("M1 T-ALL n=" + str(len(kg_t)))
    log("M1 sample head " + str(kg_t["sample"].head(6).tolist()))
    if "k_sample" in kg_t.columns:
        log("M1 k_sample head " + str(kg_t["k_sample"].head(6).tolist()))
    if "Sample ID" in kg_t.columns:
        log("M1 Sample ID head " + str(kg_t["Sample ID"].head(6).tolist()))
    log("M3 n=" + str(len(lab)) + " etp=" + str(lab["etp"].value_counts().to_dict()))

    zpath = DATA / "StJude_ALL_Pharmacotype/pharmacotyping_ped_rnaseq_fpkm.zip"
    with zipfile.ZipFile(zpath) as z:
        info = max(z.infolist(), key=lambda x: x.file_size)
        log("zip entry " + info.filename)
        with z.open(info.filename) as fh:
            peek = fh.read(800).decode("utf-8", "replace")
        sep = "," if peek.splitlines()[0].count(",") > peek.splitlines()[0].count("\t") else "\t"
        log("FPKM sep=" + repr(sep) + " peek_cols=" + peek.splitlines()[0][:180])
        with z.open(info.filename) as fh:
            mat = pd.read_csv(fh, sep=sep)
    log("FPKM shape " + str(mat.shape) + " first cols " + str(list(mat.columns)[:8]))

    gene_col = "GeneName" if "GeneName" in mat.columns else mat.columns[0]
    for c in mat.columns[:6]:
        if mat[c].astype(str).str.upper().isin({"ZEB1", "LMO2"}).any():
            gene_col = c
            break
    log("gene_col " + str(gene_col))
    mat["gene_u"] = mat[gene_col].astype(str).str.upper().str.replace(r"\.\d+$", "", regex=True)
    sub = mat[mat["gene_u"].isin(PANEL)].copy()
    log("panel genes found " + str(sorted(sub["gene_u"].unique().tolist())))
    num = [c for c in sub.columns if c not in {gene_col, "gene_u"} and pd.api.types.is_numeric_dtype(sub[c])]
    expr = sub.set_index("gene_u")[num].T.reset_index().rename(columns={"index": "fpkm_id"})
    expr["k"] = expr["fpkm_id"].map(nid)
    log("expr samples " + str(len(expr)) + " example " + str(expr["fpkm_id"].head(8).tolist()))

    kg_t["k_patient"] = kg_t["Patient ID"].map(nid) if "Patient ID" in kg_t.columns else kg_t["sample"].map(nid)
    kg_t["k_sample"] = kg_t["Sample ID"].map(nid) if "Sample ID" in kg_t.columns else (
        kg_t["k_sample"].map(nid) if "k_sample" in kg_t.columns else kg_t["sample"].map(nid)
    )
    kg_t["k_row"] = kg_t["sample"].map(nid)
    lab["k_patient"] = lab["Patient ID"].map(nid) if "Patient ID" in lab.columns else lab.iloc[:, 0].map(nid)

    cand = {
        "k_sample": nmatch(expr["k"], kg_t["k_sample"]),
        "k_patient": nmatch(expr["k"], kg_t["k_patient"]),
        "k_row": nmatch(expr["k"], kg_t["k_row"]),
    }
    log("overlap " + str(cand))
    key = max(cand, key=cand.get)
    log("using " + key + " n=" + str(cand[key]))
    meta_cols = [c for c in [key, "sample", "Patient ID", "disease", "Immunophenotype"] if c in kg_t.columns]
    meta_cols = list(dict.fromkeys(meta_cols))
    df = expr.merge(kg_t.drop_duplicates(key)[meta_cols], left_on="k", right_on=key, how="inner")
    log("matched expr-kg " + str(len(df)))
    if len(df) < 20:
        raise SystemExit("Pharmacotype FPKM still unmatched after Sample ID / Patient ID tries")

    lab_keep = [c for c in ["k_patient", "subtype", "etp", "age_years"] if c in lab.columns]
    if "Patient ID" in df.columns:
        df["k_patient"] = df["Patient ID"].map(nid)
    elif "k_patient" not in df.columns:
        df["k_patient"] = df["sample"].map(nid)
    df = df.merge(lab.drop_duplicates("k_patient")[lab_keep], on="k_patient", how="left")
    log("etp after merge " + str(df["etp"].value_counts(dropna=False).to_dict() if "etp" in df.columns else None))

    present_early = [g for g in EARLY if g in df.columns]
    present_late = [g for g in LATE if g in df.columns]
    log("early " + str(present_early) + " late " + str(present_late))
    logm = np.log2(df[present_early + present_late].astype(float) + 1)
    z = (logm - logm.mean()) / logm.std(ddof=0)
    df["maturity_literature"] = z[present_late].mean(axis=1) - z[present_early].mean(axis=1)
    for g in ["ZEB1", "ZEB2", "LMO2"]:
        df[g + "_log"] = np.log2(df[g].astype(float) + 1)
    df["cohort"] = "StJude_Pharmacotype"
    df.to_csv(PROC / "pharmacotype_maturity_panel.tsv", sep="\t", index=False)
    r = spearmanr(df["ZEB1_log"], df["maturity_literature"], nan_policy="omit")
    log(f"Pharmacotype ZEB1 vs maturity rho={r.statistic:.3f} p={r.pvalue:.3g} n={len(df)}")
    log("DONE")


if __name__ == "__main__":
    main()
