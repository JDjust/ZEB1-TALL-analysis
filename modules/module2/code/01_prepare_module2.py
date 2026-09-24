#!/usr/bin/env python3
"""Prepare Module 2: ZEB1 along normal T-cell development and adult T-ALL projection.

Primary matrix: GSE142522 Table S1 (DESeq2 log2-transformed counts; Cieslak et al.
PMID 32667968). Sorted thymocyte stages have n=1, so stage DEG/GAM are descriptive.
Independent replication uses TARGET STAR TPM and Pharmacotype FPKM marker panels.
"""
from __future__ import annotations

import math
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
M1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
M3 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2")
PROC = ROOT / "processed"
PROC.mkdir(parents=True, exist_ok=True)

STAGE_ORDER = ["CD34plus", "ISP", "EC", "LC", "SP4", "SP8"]
STAGE_INDEX = {s: i + 1 for i, s in enumerate(STAGE_ORDER)}
STAGE_LABEL = {
    "CD34plus": "CD34+",
    "ISP": "ISP",
    "EC": "Early cortical",
    "LC": "Late cortical",
    "SP4": "SP4",
    "SP8": "SP8",
}
THYMUS = ["Thy_2ng", "Thy_5ng", "Thy_10ng"]
KEY_GENES = [
    "ZEB1", "ZEB2", "LMO2", "IL7R", "TAL1", "LYL1", "HOXA9", "HOXA10",
    "TLX1", "TLX3", "PTCRA", "RAG1", "RAG2", "CD34", "KIT", "SPI1", "MEF2C",
    "CD1A", "CD1B", "CD3E", "CD3D", "CD4", "CD8A", "CD8B", "NOTCH1", "MYC",
    "BCL11B", "GATA3", "TCF7", "LEF1", "DNTT", "CD44", "CD7", "CD5", "CD2",
    "LCK", "SATB1", "RUNX1", "HES1", "PTEN", "TRAC",
]
EARLY_SIG = ["CD34", "KIT", "IL7R", "SPI1", "LMO2", "LYL1", "MEF2C", "CD44", "HES1"]
LATE_SIG = ["CD1A", "CD3E", "CD3D", "CD4", "CD8A", "RAG1", "BCL11B", "TCF7", "LEF1", "DNTT"]
HEAT_GENES = [
    "CD34", "KIT", "SPI1", "MEF2C", "LMO2", "LYL1", "IL7R", "HES1",
    "TAL1", "PTCRA", "CD1A", "RAG1", "DNTT", "CD4", "CD8A", "CD3E",
    "BCL11B", "TCF7", "NOTCH1", "MYC", "HOXA9", "TLX1", "ZEB1", "ZEB2",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def hedges_g(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float]:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return (np.nan, np.nan, np.nan, np.nan)
    sp = math.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return (0.0, 0.0, 0.0, 0.0)
    d = (a.mean() - b.mean()) / sp
    j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    g = j * d
    se = math.sqrt((n1 + n2) / (n1 * n2) + d * d / (2.0 * (n1 + n2))) * j
    return g, se, g - 1.96 * se, g + 1.96 * se


def wilcox_p(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import mannwhitneyu

    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def spearman(x, y) -> tuple[float, float]:
    from scipy.stats import spearmanr

    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return (np.nan, np.nan)
    r = spearmanr(x[m], y[m])
    return float(r.statistic), float(r.pvalue)


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=0).replace(0, np.nan)
    return mat.sub(mu, axis=0).div(sd, axis=0)


def parse_table_s1() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates = [
        DATA / "GSE142522/suppl/GSE142522_Table_S1.tsv.gz",
        DATA / "GSE142522/GSE142522_Table_S1.tsv.gz",
        DATA / "GSE142522/suppl/GSE142522_Table_S1.tsv",
        DATA / "GSE142522/GSE142522_Table_S1.tsv",
        DATA / "GSE142522/suppl/GSE142522_Table_S1.xlsx",
    ]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        raise SystemExit("GSE142522 Table S1 not found")
    log(f"reading {path}")
    if path.suffix.lower() == ".xlsx":
        raw = pd.read_excel(path, header=None, engine="openpyxl")
    else:
        raw = pd.read_csv(path, sep="\t", header=None)
    subtype_row = raw.iloc[4].tolist()
    onco_row = raw.iloc[5].tolist()
    header = [str(x).strip() if pd.notna(x) else "" for x in raw.iloc[6].tolist()]
    data = raw.iloc[7:].copy()
    data.columns = header
    data = data.rename(columns={"ID": "gene"})
    data["gene"] = data["gene"].astype(str).str.strip()
    data["gene_u"] = data["gene"].str.upper().str.replace(r"\.\d+$", "", regex=True)
    data["Gene_type"] = data["Gene_type"].astype(str)
    sample_cols = header[6:]
    meta_rows = []
    for i, s in enumerate(sample_cols, start=6):
        sub = str(subtype_row[i]).strip() if i < len(subtype_row) and pd.notna(subtype_row[i]) else ""
        onco = str(onco_row[i]).strip() if i < len(onco_row) and pd.notna(onco_row[i]) else ""
        if onco in {"NA", "nan", "None", ""}:
            onco = "unknown"
        if s in STAGE_ORDER:
            group, stage = "sorted_thymocyte", s
        elif s in THYMUS:
            group, stage = "total_thymus", "total_thymus"
        else:
            group, stage = "adult_T-ALL", "T-ALL"
        meta_rows.append(
            {
                "sample": s,
                "group": group,
                "stage": stage,
                "stage_label": STAGE_LABEL.get(stage, stage),
                "stage_index": STAGE_INDEX.get(stage, np.nan),
                "immunophenotype": sub if group == "adult_T-ALL" else "Normal",
                "oncogene": onco if group == "adult_T-ALL" else "Normal",
                "is_tall": group == "adult_T-ALL",
            }
        )
    meta = pd.DataFrame(meta_rows)
    expr = data.set_index("gene_u")[sample_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr[~expr.index.duplicated(keep="first")]
    de = data[["gene", "gene_u", "Gene_type", "Localisation",
               "log2FoldChange(T-ALL_Vs_Normal)", "pvalue", "padj"]].copy()
    de = de.rename(columns={
        "log2FoldChange(T-ALL_Vs_Normal)": "log2FC_tall_vs_normal",
        "pvalue": "pval_tall_vs_normal",
        "padj": "fdr_tall_vs_normal",
    })
    for c in ["log2FC_tall_vs_normal", "pval_tall_vs_normal", "fdr_tall_vs_normal"]:
        de[c] = pd.to_numeric(de[c], errors="coerce")
    de = de.drop_duplicates("gene_u")
    log(f"genes {expr.shape[0]} samples {expr.shape[1]}")
    log(meta.groupby(["group", "immunophenotype"]).size().to_string())
    return meta, expr, de


def keygene_long(meta: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    present = [g for g in KEY_GENES if g in expr.index]
    miss = [g for g in KEY_GENES if g not in expr.index]
    if miss:
        log("missing key genes " + str(miss))
    long = expr.loc[present].T
    long.index.name = "sample"
    long = long.reset_index().merge(meta, on="sample", how="left")
    return long


def stage_gene_stats(expr: pd.DataFrame) -> pd.DataFrame:
    stage_mat = expr[STAGE_ORDER]
    idx = np.array([STAGE_INDEX[s] for s in STAGE_ORDER], float)
    rows = []
    for gene, vals in stage_mat.iterrows():
        y = vals.to_numpy(float)
        if not np.isfinite(y).all():
            continue
        rho, p = spearman(idx, y)
        rows.append({
            "gene": gene,
            "rho_stage": rho,
            "p_stage": p,
            "cd34": float(vals["CD34plus"]),
            "isp": float(vals["ISP"]),
            "ec": float(vals["EC"]),
            "lc": float(vals["LC"]),
            "sp4": float(vals["SP4"]),
            "sp8": float(vals["SP8"]),
            "early_mean": float(vals[["CD34plus", "ISP"]].mean()),
            "cortical_mean": float(vals[["EC", "LC"]].mean()),
            "sp_mean": float(vals[["SP4", "SP8"]].mean()),
            "range": float(np.nanmax(y) - np.nanmin(y)),
        })
    out = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests

    mask = out["p_stage"].notna()
    out["fdr_stage"] = np.nan
    if mask.any():
        _, fdr, _, _ = multipletests(out.loc[mask, "p_stage"], method="fdr_bh")
        out.loc[mask, "fdr_stage"] = fdr
    return out.sort_values("rho_stage", ascending=False)


def scores_for_samples(meta: pd.DataFrame, expr: pd.DataFrame, stage_stats: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    early = [g for g in EARLY_SIG if g in expr.index]
    late = [g for g in LATE_SIG if g in expr.index]
    z = zscore_rows(expr.loc[sorted(set(early + late))])
    lit = (z.loc[late].mean(axis=0) - z.loc[early].mean(axis=0)).rename("maturity_literature")

    ranked = stage_stats.dropna(subset=["rho_stage"]).sort_values("rho_stage", ascending=False)
    up = ranked.head(80)["gene"].tolist()
    down = ranked.tail(80)["gene"].tolist()
    emp_genes = [g for g in up + down if g in expr.index]
    z2 = zscore_rows(expr.loc[emp_genes])
    emp = (
        z2.loc[[g for g in up if g in z2.index]].mean(axis=0)
        - z2.loc[[g for g in down if g in z2.index]].mean(axis=0)
    ).rename("maturity_empirical")

    sc = meta.copy()
    sc = sc.merge(lit.reset_index().rename(columns={"index": "sample"}), on="sample", how="left")
    sc = sc.merge(emp.reset_index().rename(columns={"index": "sample"}), on="sample", how="left")
    for g in ["ZEB1", "ZEB2", "LMO2", "IL7R"]:
        if g in expr.index:
            sc[g] = expr.loc[g, sc["sample"]].to_numpy()
    sc["ratio_zeb"] = sc["ZEB1"] - sc["ZEB2"]
    sc["ratio_lmo"] = sc["LMO2"] - sc["ZEB1"]
    return sc, _sig_table(up, down, early, late)


def _sig_table(up, down, early, late) -> pd.DataFrame:
    rows = [{"gene": g, "set": "empirical_late"} for g in up]
    rows += [{"gene": g, "set": "empirical_early"} for g in down]
    rows += [{"gene": g, "set": "literature_early"} for g in early]
    rows += [{"gene": g, "set": "literature_late"} for g in late]
    return pd.DataFrame(rows)


def pca_table(meta: pd.DataFrame, expr: pd.DataFrame, de: pd.DataFrame) -> pd.DataFrame:
    pc = expr.copy()
    if "Gene_type" in de.columns:
        coding = set(de.loc[de["Gene_type"].str.contains("mRNA|protein", case=False, na=False), "gene_u"])
        if len(coding) > 1000:
            pc = pc.loc[pc.index.intersection(coding)]
    pc = pc.loc[pc.notna().mean(axis=1) > 0.95]
    v = pc.var(axis=1, ddof=1)
    top = v.nlargest(min(2000, len(v))).index
    X = pc.loc[top].T.to_numpy(float)
    X = np.nan_to_num(X, nan=np.nanmedian(X))
    X = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))
    u, s, _vt = np.linalg.svd(X, full_matrices=False)
    var = (s ** 2) / (s ** 2).sum()
    out = meta.copy()
    out["PC1"] = u[:, 0] * s[0]
    out["PC2"] = u[:, 1] * s[1]
    out["PC3"] = u[:, 2] * s[2]
    out.attrs["var"] = var[:3]
    log(f"PCA var {var[0]:.3f} {var[1]:.3f} {var[2]:.3f} genes={len(top)}")
    Path(PROC / "pca_variance.txt").write_text(
        f"PC1\t{var[0]}\nPC2\t{var[1]}\nPC3\t{var[2]}\nn_genes\t{len(top)}\n", encoding="utf-8"
    )
    return out


def nearest_stage(meta: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    stage_mat = expr[STAGE_ORDER]
    v = stage_mat.var(axis=1, ddof=1)
    top = v.nlargest(1500).index
    ref = stage_mat.loc[top]
    tall = meta.loc[meta["is_tall"], "sample"].tolist()
    rows = []
    for s in tall:
        x = expr.loc[top, s]
        cors = {st: float(x.corr(ref[st], method="pearson")) for st in STAGE_ORDER}
        best = max(cors, key=cors.get)
        rows.append({
            "sample": s,
            "nearest_stage": best,
            "nearest_label": STAGE_LABEL[best],
            "nearest_index": STAGE_INDEX[best],
            **{f"cor_{k}": v for k, v in cors.items()},
        })
    out = pd.DataFrame(rows).merge(meta, on="sample", how="left")
    return out


def regressions(sc: pd.DataFrame) -> pd.DataFrame:
    import statsmodels.api as sm

    d = sc[sc["group"].isin(["sorted_thymocyte", "adult_T-ALL"])].copy()
    d["is_tall_n"] = d["is_tall"].astype(int)
    rows = []

    def fit(formula_y, cols, subset, name):
        dd = subset.dropna(subset=[formula_y] + cols)
        if len(dd) < 10:
            return
        y = dd[formula_y].to_numpy(float)
        x = sm.add_constant(dd[cols].to_numpy(float))
        m = sm.OLS(y, x).fit()
        for i, c in enumerate(["intercept"] + cols):
            rows.append({
                "model": name, "term": c, "n": int(m.nobs),
                "coef": float(m.params[i]), "se": float(m.bse[i]),
                "ci95_lo": float(m.conf_int()[i, 0]), "ci95_hi": float(m.conf_int()[i, 1]),
                "p": float(m.pvalues[i]), "r2": float(m.rsquared),
            })

    if "maturity_literature" in d.columns:
        fit("ZEB1", ["is_tall_n"], d, "ZEB1 ~ T-ALL")
        fit("ZEB1", ["maturity_literature"], d, "ZEB1 ~ maturity")
        fit("ZEB1", ["is_tall_n", "maturity_literature"], d, "ZEB1 ~ T-ALL + maturity")
        fit("maturity_literature", ["is_tall_n"], d, "maturity ~ T-ALL")
        fit("maturity_literature", ["ZEB1"], d, "maturity ~ ZEB1")
        fit("maturity_literature", ["is_tall_n", "ZEB1"], d, "maturity ~ T-ALL + ZEB1")
        tall = d[d["is_tall"]]
        fit("ZEB1", ["maturity_literature"], tall, "T-ALL only: ZEB1 ~ maturity")
        if "maturity_empirical" in d.columns:
            fit("ZEB1", ["is_tall_n", "maturity_empirical"], d, "ZEB1 ~ T-ALL + empirical_maturity")
            fit("ZEB1", ["maturity_empirical"], tall, "T-ALL only: ZEB1 ~ empirical_maturity")
    return pd.DataFrame(rows)


def extract_star_panel(genes: list[str]) -> pd.DataFrame:
    rna_dir = DATA / "TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts"
    fmap = pd.read_csv(M1 / "target_gdc_file_map.tsv", sep="\t")
    want = set(genes)
    recs = []
    files = sorted(p for p in rna_dir.iterdir() if p.suffix == ".tsv")
    fmap_idx = fmap.drop_duplicates("file_name").set_index("file_name")
    log(f"TARGET extracting {len(want)} genes from {len(files)} files")
    for i, p in enumerate(files, 1):
        out = {}
        with p.open(encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#") or line.startswith("gene_id") or line.startswith("N_"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                name = parts[1]
                if name in want:
                    try:
                        out[name] = float(parts[6])
                    except ValueError:
                        out[name] = np.nan
                    if len(out) == len(want):
                        break
        rec = {"file_name": p.name, **out}
        if p.name in fmap_idx.index:
            rec.update(fmap_idx.loc[p.name].to_dict())
        recs.append(rec)
        if i % 80 == 0:
            log(f"  {i}/{len(files)}")
    return pd.DataFrame(recs)


def usi_from_barcode(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x)
    m = re.search(r"TARGET-\d+-([A-Z0-9]+)", s)
    if m:
        return m.group(1)
    return s


def prepare_target_panel(genes: list[str]) -> pd.DataFrame:
    df = extract_star_panel(genes)
    df["usi"] = df["sample_submitter_id"].map(usi_from_barcode)
    df["usi2"] = df["case_submitter_id"].map(usi_from_barcode) if "case_submitter_id" in df.columns else ""
    df["usi"] = np.where(df["usi"].astype(str).isin(["", "nan", "None"]), df["usi2"], df["usi"])
    df = df[df["sample_type"].astype(str).str.startswith("Primary")].copy()
    df["pref"] = np.where(df["sample_type"].astype(str).str.contains("Bone Marrow"), 1, 2)
    df = df.sort_values(["usi", "pref"]).drop_duplicates("usi")
    lab = pd.read_csv(M3 / "target_tall_subtype.tsv", sep="\t")
    keep = [c for c in ["usi", "subtype", "etp", "ZEB1", "ZEB2", "LMO2"] if c in lab.columns]
    df = df.merge(lab[keep], on="usi", how="inner", suffixes=("", "_lab"))
    present_early = [g for g in EARLY_SIG if g in df.columns]
    present_late = [g for g in LATE_SIG if g in df.columns]
    logm = np.log2(df[present_early + present_late].astype(float) + 1)
    z = (logm - logm.mean()) / logm.std(ddof=0)
    df["maturity_literature"] = z[present_late].mean(axis=1) - z[present_early].mean(axis=1)
    if "ZEB1" in df.columns:
        df["ZEB1_log"] = np.log2(df["ZEB1"] + 1)
        df["ZEB2_log"] = np.log2(df["ZEB2"] + 1)
        df["LMO2_log"] = np.log2(df["LMO2"] + 1)
    df["cohort"] = "TARGET-ALL-P2"
    log(f"TARGET panel n={len(df)} subtype={df['subtype'].value_counts().to_dict() if 'subtype' in df.columns else None}")
    return df


def _norm_id(s) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def prepare_pharmacotype_panel(genes: list[str]) -> pd.DataFrame:
    zpath = DATA / "StJude_ALL_Pharmacotype/pharmacotyping_ped_rnaseq_fpkm.zip"
    log(f"reading {zpath}")
    with zipfile.ZipFile(zpath) as z:
        info = max(z.infolist(), key=lambda x: x.file_size)
        with z.open(info.filename) as fh:
            peek = fh.read(800).decode("utf-8", "replace")
        sep = "," if peek.splitlines()[0].count(",") > peek.splitlines()[0].count("\t") else "\t"
        with z.open(info.filename) as fh:
            mat = pd.read_csv(fh, sep=sep)
    gene_col = "GeneName" if "GeneName" in mat.columns else mat.columns[0]
    for c in mat.columns[:4]:
        if mat[c].astype(str).str.upper().isin({"ZEB1", "LMO2"}).any():
            gene_col = c
            break
    mat["gene_u"] = mat[gene_col].astype(str).str.upper().str.replace(r"\.\d+$", "", regex=True)
    want = {g.upper() for g in genes}
    sub = mat[mat["gene_u"].isin(want)].copy()
    num = [c for c in sub.columns if c not in {gene_col, "gene_u"} and pd.api.types.is_numeric_dtype(sub[c])]
    expr = sub.set_index("gene_u")[num].T
    expr.index.name = "sample"
    expr = expr.reset_index()
    expr["k"] = expr["sample"].map(_norm_id)
    lab = pd.read_csv(M3 / "pharmacotype_tall_subtype.tsv", sep="\t")
    idcol = "sample" if "sample" in lab.columns else lab.columns[0]
    lab["k"] = lab[idcol].map(_norm_id)
    keep = [c for c in ["k", idcol, "subtype", "etp", "ZEB1", "ZEB2", "LMO2"] if c in lab.columns or c == "k"]
    keep = list(dict.fromkeys(keep))
    df = expr.merge(lab[keep].drop_duplicates("k"), on="k", how="inner")
    if len(df) < 20:
        m1 = pd.read_csv(M1 / "pharmacotype_keygenes.tsv", sep="\t")
        m1 = m1[m1["disease"].eq("T-ALL")].copy()
        sid = "sample" if "sample" in m1.columns else m1.columns[0]
        m1["k"] = m1[sid].map(_norm_id)
        df = expr.merge(m1[["k", sid, "ZEB1", "ZEB2", "LMO2"]].drop_duplicates("k"), on="k", how="inner")
        df = df.merge(lab[keep].drop_duplicates("k"), on="k", how="left", suffixes=("", "_lab"))
    log(f"Pharmacotype FPKM matched n={len(df)}")
    present_early = [g for g in EARLY_SIG if g in df.columns]
    present_late = [g for g in LATE_SIG if g in df.columns]
    logm = np.log2(df[present_early + present_late].astype(float) + 1)
    z = (logm - logm.mean()) / logm.std(ddof=0)
    df["maturity_literature"] = z[present_late].mean(axis=1) - z[present_early].mean(axis=1)
    if "ZEB1" in df.columns:
        if df["ZEB1"].max() > 50:
            df["ZEB1_log"] = np.log2(df["ZEB1"] + 1)
            df["ZEB2_log"] = np.log2(df["ZEB2"] + 1)
            df["LMO2_log"] = np.log2(df["LMO2"] + 1)
        else:
            df["ZEB1_log"] = df["ZEB1"]
            df["ZEB2_log"] = df["ZEB2"]
            df["LMO2_log"] = df["LMO2"]
    df["cohort"] = "StJude_Pharmacotype"
    log(f"Pharmacotype panel n={len(df)}")
    return df


def contrast_table(sc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tall = sc.loc[sc["is_tall"], "ZEB1"].to_numpy()
    stages = sc.loc[sc["group"] == "sorted_thymocyte", "ZEB1"].to_numpy()
    thymus = sc.loc[sc["group"] == "total_thymus", "ZEB1"].to_numpy()
    normal = sc.loc[~sc["is_tall"], "ZEB1"].to_numpy()

    def add(name, a, b, ga, gb):
        g, se, lo, hi = hedges_g(a, b)
        rows.append({
            "contrast": name, "group_a": ga, "n_a": int(np.isfinite(a).sum()),
            "mean_a": float(np.nanmean(a)), "group_b": gb, "n_b": int(np.isfinite(b).sum()),
            "mean_b": float(np.nanmean(b)), "hedges_g": g, "se": se,
            "ci95_lo": lo, "ci95_hi": hi, "wilcox_p": wilcox_p(a, b),
        })

    add("T-ALL vs all normal", tall, normal, "T-ALL", "normal")
    add("T-ALL vs sorted stages", tall, stages, "T-ALL", "sorted_thymocyte")
    add("T-ALL vs total thymus", tall, thymus, "T-ALL", "total_thymus")
    for ip in ["Immature", "Pre_ab_Cortical", "TCR_Mature"]:
        a = sc.loc[sc["immunophenotype"] == ip, "ZEB1"].to_numpy()
        b = sc.loc[sc["is_tall"] & (sc["immunophenotype"] != ip), "ZEB1"].to_numpy()
        add(f"{ip} vs other T-ALL", a, b, ip, "other_T-ALL")
    for gene in ["ZEB1", "ZEB2", "LMO2"]:
        if gene not in sc.columns:
            continue
        a = sc.loc[sc["is_tall"], gene].to_numpy()
        b = sc.loc[sc["group"] == "sorted_thymocyte", gene].to_numpy()
        g, se, lo, hi = hedges_g(a, b)
        rows.append({
            "contrast": f"T-ALL vs stages ({gene})", "group_a": "T-ALL",
            "n_a": int(np.isfinite(a).sum()), "mean_a": float(np.nanmean(a)),
            "group_b": "sorted_thymocyte", "n_b": int(np.isfinite(b).sum()),
            "mean_b": float(np.nanmean(b)), "hedges_g": g, "se": se,
            "ci95_lo": lo, "ci95_hi": hi, "wilcox_p": wilcox_p(a, b), "gene": gene,
        })
    out = pd.DataFrame(rows)
    if "wilcox_p" in out.columns:
        from statsmodels.stats.multitest import multipletests
        mask = out["wilcox_p"].notna()
        if mask.any():
            _, fdr, _, _ = multipletests(out.loc[mask, "wilcox_p"], method="fdr_bh")
            out.loc[mask, "fdr"] = fdr
    return out


def main() -> None:
    meta, expr, de = parse_table_s1()
    meta.to_csv(PROC / "sample_metadata.tsv", sep="\t", index=False)
    expr.to_csv(PROC / "expression_log2.tsv", sep="\t")
    de.to_csv(PROC / "deseq2_tall_vs_normal.tsv", sep="\t", index=False)

    kg = keygene_long(meta, expr)
    kg.to_csv(PROC / "keygenes_samples.tsv", sep="\t", index=False)

    log("stage Spearman across 6 stages (n=1 each, descriptive)")
    st = stage_gene_stats(expr)
    st.to_csv(PROC / "stage_gene_spearman.tsv", sep="\t", index=False)
    log(st.head(8)[["gene", "rho_stage", "p_stage"]].to_string(index=False))
    log("ZEB1/ZEB2/LMO2 stage rho:")
    log(st[st["gene"].isin(["ZEB1", "ZEB2", "LMO2"])][
        ["gene", "rho_stage", "p_stage", "cd34", "isp", "ec", "lc", "sp4", "sp8", "range"]
    ].to_string(index=False))

    sc, sig = scores_for_samples(meta, expr, st)
    # fix scores_for_samples return if I left a bug
    sc.to_csv(PROC / "sample_scores.tsv", sep="\t", index=False)
    sig.to_csv(PROC / "signature_genes.tsv", sep="\t", index=False)

    pca = pca_table(meta, expr, de)
    pca.to_csv(PROC / "pca_samples.tsv", sep="\t", index=False)
    near = nearest_stage(meta, expr)
    near.to_csv(PROC / "nearest_stage.tsv", sep="\t", index=False)
    log("nearest stage counts\n" + near["nearest_label"].value_counts().to_string())

    reg = regressions(sc)
    reg.to_csv(PROC / "regressions.tsv", sep="\t", index=False)
    log(reg.to_string(index=False))

    ct = contrast_table(sc)
    ct.to_csv(PROC / "zeb1_contrasts.tsv", sep="\t", index=False)

    heat = expr.loc[[g for g in HEAT_GENES if g in expr.index], STAGE_ORDER + THYMUS]
    heat.to_csv(PROC / "heatmap_stage_genes.tsv", sep="\t")

    panel = sorted(set(KEY_GENES + EARLY_SIG + LATE_SIG))
    tgt = prepare_target_panel(panel)
    tgt.to_csv(PROC / "target_maturity_panel.tsv", sep="\t", index=False)
    rho, p = spearman(tgt["ZEB1_log"], tgt["maturity_literature"])
    log(f"TARGET ZEB1 vs maturity rho={rho:.3f} p={p:.3g} n={len(tgt)}")

    ph = prepare_pharmacotype_panel(panel)
    ph.to_csv(PROC / "pharmacotype_maturity_panel.tsv", sep="\t", index=False)
    if "ZEB1_log" in ph.columns:
        rho2, p2 = spearman(ph["ZEB1_log"], ph["maturity_literature"])
        log(f"Pharmacotype ZEB1 vs maturity rho={rho2:.3f} p={p2:.3g} n={len(ph)}")

    note = (
        "GSE142521/GSE142522 are ABI SOLiD color-space RNA-seq. Analysis uses the authors'\n"
        "DESeq2 log2-transformed counts (GEO GSE142522 Table S1), not re-quantified FASTQ.\n"
        "Sorted thymocyte stages have n=1 biological sample each. Stage trends, GAM and\n"
        "stage-wise Spearman are therefore descriptive. Adult T-ALL n=41 is the statistical unit\n"
        "for T-ALL vs maturity / subtype tests. Independent replication uses TARGET and\n"
        "Pharmacotype literature maturation scores, not the same SOLiD matrix.\n"
    )
    (PROC / "M2_data_note.txt").write_text(note, encoding="utf-8")
    log("MODULE2 PREPARE DONE")


if __name__ == "__main__":
    main()
