#!/usr/bin/env python3
"""Module 4 prepare: ZEB1-continuous transcriptome in TARGET and Pharmacotype.

ZEB1 stays continuous. Q1 vs Q4 is display-only DEG, not a hard biological class.
No MTX/LMO2 causal story is imposed.
"""
from __future__ import annotations

import math
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
M1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
M3 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4")
PROC = ROOT / "processed"
PROC.mkdir(parents=True, exist_ok=True)
STAR = DATA / "TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts"


def log(m):
    print(m, flush=True)


def nid(s) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def hedges_g(a, b) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 8 or n2 < 8:
        return np.nan
    sp = math.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return 0.0
    d = (a.mean() - b.mean()) / sp
    j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    return j * d


def usi_from_barcode(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x)
    m = re.search(r"TARGET-\d+-([A-Z0-9]+)", s)
    return m.group(1) if m else s


def fdr(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    mask = np.isfinite(p)
    if mask.any():
        _, q, _, _ = multipletests(p[mask], method="fdr_bh")
        out[mask] = q
    return out


def load_target_labels() -> pd.DataFrame:
    lab = pd.read_csv(M3 / "target_tall_subtype.tsv", sep="\t")
    if "file_name" not in lab.columns:
        fmap = pd.read_csv(M1 / "target_gdc_file_map.tsv", sep="\t")
        fmap["usi"] = fmap["sample_submitter_id"].map(usi_from_barcode)
        fmap = fmap[fmap["sample_type"].astype(str).str.startswith("Primary")].copy()
        fmap["pref"] = np.where(fmap["sample_type"].astype(str).str.contains("Bone Marrow"), 1, 2)
        fmap = fmap.sort_values(["usi", "pref"]).drop_duplicates("usi")
        lab = lab.merge(fmap[["usi", "file_name", "sample_submitter_id", "sample_type"]], on="usi", how="inner")
    if "ZEB1" not in lab.columns:
        raise SystemExit("ZEB1 missing from target_tall_subtype")
    lab["ZEB1_log"] = np.log2(lab["ZEB1"].astype(float) + 1)
    lab["q"] = pd.qcut(lab["ZEB1_log"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    log(f"TARGET labels n={len(lab)} files={lab['file_name'].nunique()} subtype={lab['subtype'].value_counts().to_dict() if 'subtype' in lab.columns else None}")
    return lab


def extract_target_matrix(lab: pd.DataFrame) -> pd.DataFrame:
    """protein_coding TPM, genes x samples (usi)."""
    want_files = {str(x): usi for usi, x in zip(lab["usi"], lab["file_name"])}
    cols = {}
    gene_type = None
    n = 0
    files = sorted(STAR.glob("*.tsv"))
    log(f"STAR files {len(files)}; wanted {len(want_files)}")
    for i, p in enumerate(files, 1):
        usi = want_files.get(p.name)
        if usi is None:
            continue
        names, vals, types = [], [], []
        with p.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#") or line.startswith("gene_id") or line.startswith("N_"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7:
                    continue
                names.append(parts[1])
                types.append(parts[2])
                try:
                    vals.append(float(parts[6]))
                except ValueError:
                    vals.append(np.nan)
        s = pd.Series(vals, index=pd.Index(names, name="gene"))
        if gene_type is None:
            gene_type = pd.Series(types, index=s.index, name="gene_type")
        cols[usi] = s
        n += 1
        if n % 40 == 0:
            log(f"  extracted {n}/{len(want_files)}")
    mat = pd.DataFrame(cols)
    if gene_type is not None:
        keep = gene_type.eq("protein_coding")
        mat = mat.loc[keep[keep].index.intersection(mat.index)]
    mat = mat[~mat.index.duplicated(keep="first")]
    log(f"TARGET protein-coding matrix {mat.shape}")
    return mat


def filter_expr(mat: pd.DataFrame, min_tpm=1.0, min_frac=0.2) -> pd.DataFrame:
    frac = (mat > min_tpm).mean(axis=1)
    out = mat.loc[frac >= min_frac]
    log(f"filter TPM>{min_tpm} in >= {min_frac:.0%} samples: {out.shape}")
    return out


def stats_vs_zeb1(logmat: pd.DataFrame, zeb1: pd.Series, q: pd.Series, cohort: str) -> pd.DataFrame:
    common = logmat.columns.intersection(zeb1.index)
    X = logmat.loc[:, common]
    z = zeb1.loc[common]
    qq = q.loc[common]
    q1 = qq[qq.astype(str).eq("Q1")].index
    q4 = qq[qq.astype(str).eq("Q4")].index
    rows = []
    z_np = z.to_numpy(float)
    for gene, row in X.iterrows():
        y = row.to_numpy(float)
        mask = np.isfinite(y) & np.isfinite(z_np)
        if mask.sum() < 20:
            continue
        rho, p_rho = spearmanr(z_np[mask], y[mask])
        a = row.reindex(q4).to_numpy(float)
        b = row.reindex(q1).to_numpy(float)
        a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
        if len(a) >= 8 and len(b) >= 8:
            try:
                p_w = mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p_w = np.nan
            g = hedges_g(a, b)
            lfc = float(np.nanmean(a) - np.nanmean(b))
        else:
            p_w, g, lfc = np.nan, np.nan, np.nan
        rows.append({
            "gene": gene, "cohort": cohort, "n": int(mask.sum()),
            "rho_zeb1": float(rho), "p_rho": float(p_rho),
            "lfc_Q4_minus_Q1": lfc, "hedges_g_Q4_vs_Q1": g, "p_wilcox_Q4_vs_Q1": p_w,
            "mean_Q1": float(np.nanmean(b)) if len(b) else np.nan,
            "mean_Q4": float(np.nanmean(a)) if len(a) else np.nan,
        })
    out = pd.DataFrame(rows)
    out["fdr_rho"] = fdr(out["p_rho"])
    out["fdr_wilcox"] = fdr(out["p_wilcox_Q4_vs_Q1"])
    out = out.sort_values("p_rho")
    log(f"{cohort} genes {len(out)} FDR_rho<0.05 {(out.fdr_rho<0.05).sum()} FDR_Q1Q4<0.05 {(out.fdr_wilcox<0.05).sum()}")
    return out


def pharmacotype_matrix(lab_ph: pd.DataFrame) -> pd.DataFrame:
    zpath = DATA / "StJude_ALL_Pharmacotype/pharmacotyping_ped_rnaseq_fpkm.zip"
    with zipfile.ZipFile(zpath) as z:
        info = max(z.infolist(), key=lambda x: x.file_size)
        with z.open(info.filename) as fh:
            peek = fh.read(400).decode("utf-8", "replace")
        sep = "," if peek.splitlines()[0].count(",") > peek.splitlines()[0].count("\t") else "\t"
        with z.open(info.filename) as fh:
            mat = pd.read_csv(fh, sep=sep)
    gene_col = "GeneName" if "GeneName" in mat.columns else mat.columns[0]
    mat["gene_u"] = mat[gene_col].astype(str).str.upper().str.replace(r"\.\d+$", "", regex=True)
    mat = mat[~mat["gene_u"].duplicated(keep="first")]
    ids = set(lab_ph["sample"].map(str))
    num = [c for c in mat.columns if str(c) in ids]
    log(f"Pharmacotype FPKM matched columns {len(num)} / labels {len(lab_ph)}")
    expr = mat.set_index("gene_u")[num]
    expr.columns = [str(c) for c in expr.columns]
    return expr.astype(float)


def main():
    lab = load_target_labels()
    lab.to_csv(PROC / "target_samples.tsv", sep="\t", index=False)

    tmat = extract_target_matrix(lab)
    tmat = filter_expr(tmat)
    tlog = np.log2(tmat.astype(float) + 1)
    # keep ZEB1 even if filtered
    if "ZEB1" not in tlog.index:
        log("WARNING ZEB1 not in filtered TARGET matrix")
    tlog.to_csv(PROC / "target_log2tpm_filtered.tsv.gz", sep="\t", compression="gzip")
    zeb1 = lab.set_index("usi")["ZEB1_log"]
    q = lab.set_index("usi")["q"]
    st_t = stats_vs_zeb1(tlog, zeb1, q, "TARGET")
    st_t.to_csv(PROC / "target_zeb1_gene_stats.tsv", sep="\t", index=False)
    st_t.head(500).to_csv(PROC / "target_zeb1_gene_stats_top500.tsv", sep="\t", index=False)

    lab_ph = pd.read_csv(M3 / "pharmacotype_tall_subtype.tsv", sep="\t")
    lab_ph["ZEB1_log"] = np.log2(lab_ph["ZEB1"].astype(float) + 1)
    lab_ph["q"] = pd.qcut(lab_ph["ZEB1_log"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    lab_ph.to_csv(PROC / "pharmacotype_samples.tsv", sep="\t", index=False)
    pmat = pharmacotype_matrix(lab_ph)
    pmat = filter_expr(pmat, min_tpm=1.0, min_frac=0.2)
    plog = np.log2(pmat.astype(float) + 1)
    plog.to_csv(PROC / "pharmacotype_log2fpkm_filtered.tsv.gz", sep="\t", compression="gzip")
    st_p = stats_vs_zeb1(plog, lab_ph.set_index("sample")["ZEB1_log"], lab_ph.set_index("sample")["q"], "Pharmacotype")
    st_p.to_csv(PROC / "pharmacotype_zeb1_gene_stats.tsv", sep="\t", index=False)

    both = st_t[["gene", "rho_zeb1", "fdr_rho", "lfc_Q4_minus_Q1", "fdr_wilcox"]].merge(
        st_p[["gene", "rho_zeb1", "fdr_rho", "lfc_Q4_minus_Q1", "fdr_wilcox"]],
        on="gene", suffixes=("_TARGET", "_Pharmacotype"), how="inner"
    )
    both.to_csv(PROC / "cross_cohort_gene_stats.tsv", sep="\t", index=False)
    n_sig = ((both.fdr_rho_TARGET < 0.05) & (both.fdr_rho_Pharmacotype < 0.05) &
             (np.sign(both.rho_zeb1_TARGET) == np.sign(both.rho_zeb1_Pharmacotype))).sum()
    log(f"cross-cohort same-direction FDR<0.05 genes {n_sig} / {len(both)}")
    log("MODULE4 PREPARE DONE")


if __name__ == "__main__":
    main()
