#!/usr/bin/env python3
"""Gene-level ZEB1 GSEA on Pölönen T-ALL RNA-seq (teacher-paper claim check)."""
from __future__ import annotations

import json
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"d:\_bioinformation\ZEB1")
OUT = ROOT / "02_analyses" / "02_teacher_revision_polonen_gsea"
TPM = ROOT / "00_data" / "00_raw" / "TALL_X01_tpm.tsv"
GENE_MAP = ROOT / ".project-artifacts" / "work" / "20261008T051121Z-384f93b5" / "M1" / "Lim_discovery_GENE_ID_MAP.tsv"
MANIFEST = ROOT / ".project-artifacts" / "work" / "20261008T051121Z-384f93b5" / "M2" / "AUTHOR_PARTITION_MANIFEST.tsv"
GMT = ROOT / ".project-artifacts" / "work" / "20261006T010250Z-75920166" / "analysis" / "Q2_FUNCTIONAL" / "inputs" / "h.all.v2025.1.Hs.symbols.gmt"
ZEB1 = "ENSG00000148516"
SEED = 20261009


def stable_id(ens: str) -> str:
    return str(ens).split(".")[0]


def load_matrix() -> tuple[pd.DataFrame, pd.Series]:
    expr = pd.read_csv(TPM, sep="\t", index_col=0)
    expr = expr.astype(np.float32)
    expr.index = expr.index.map(stable_id)
    if expr.index.has_duplicates:
        expr = expr.groupby(level=0, sort=False).max()
    samples = expr.columns.astype(str)
    if ZEB1 not in expr.index:
        raise SystemExit(f"ZEB1 {ZEB1} not in TPM matrix")
    zeb1_tpm = expr.loc[ZEB1].astype(np.float64)
    logx = np.log2(expr + 1.0)
    keep = logx.mean(axis=1) > 0.5
    logx = logx.loc[keep]
    return logx, zeb1_tpm.loc[samples]


def collapse_to_symbols(logx: pd.DataFrame, gene_map: pd.Series) -> pd.DataFrame:
    symbols = gene_map.reindex(logx.index)
    symbols = symbols.dropna()
    symbols = symbols[symbols.astype(str).str.match(r"^[A-Za-z0-9].*") & ~symbols.astype(str).str.startswith("ENSG")]
    logx = logx.loc[symbols.index]
    var = logx.var(axis=1)
    order = var.sort_values(ascending=False).index
    chosen = symbols.loc[order].drop_duplicates(keep="first")
    out = logx.loc[chosen.index].copy()
    out.index = chosen.values
    out = out[~out.index.duplicated(keep="first")]
    out = out.drop(index=["ZEB1"], errors="ignore")
    return out


def spearman_rank(logx: pd.DataFrame, zeb1_log: pd.Series) -> pd.Series:
    z = zeb1_log.loc[logx.columns].to_numpy(dtype=np.float64)
    z = (z - z.mean()) / z.std(ddof=0)
    X = logx.to_numpy(dtype=np.float64)
    X = (X - X.mean(axis=1, keepdims=True)) / np.clip(X.std(axis=1, keepdims=True, ddof=0), 1e-12, None)
    rho = (X @ z) / X.shape[1]
    # Positive score = higher when ZEB1 is low (teacher GSEA direction).
    return pd.Series(-rho, index=logx.index, name="neg_spearman_with_ZEB1")


def median_split_t(logx: pd.DataFrame, zeb1_log: pd.Series) -> pd.Series:
    z = zeb1_log.loc[logx.columns]
    med = float(z.median())
    low = z <= med
    high = ~low
    xl = logx.loc[:, low.index[low]].to_numpy(dtype=np.float64)
    xh = logx.loc[:, high.index[high]].to_numpy(dtype=np.float64)
    tstat, _ = stats.ttest_ind(xl, xh, axis=1, equal_var=False)
    tstat = np.nan_to_num(tstat, nan=0.0)
    return pd.Series(tstat, index=logx.index, name="t_low_minus_high"), int(low.sum()), int(high.sum()), med


def run_prerank(rnk: pd.Series, tag: str) -> pd.DataFrame:
    rnk = rnk.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    rnk.to_csv(OUT / f"rank_{tag}.tsv", sep="\t", header=False)
    res = gp.prerank(
        rnk=rnk,
        gene_sets=str(GMT),
        outdir=str(OUT / f"gsea_{tag}"),
        permutation_num=1000,
        min_size=15,
        max_size=500,
        seed=SEED,
        threads=1,
        verbose=False,
    )
    tab = res.res2d.copy()
    tab.to_csv(OUT / f"hallmark_{tag}.tsv", sep="\t", index=False)
    return tab


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gene_map = pd.read_csv(GENE_MAP, sep="\t")
    gene_map["author_gene_id"] = gene_map["author_gene_id"].map(stable_id)
    gene_map = gene_map.drop_duplicates("author_gene_id").set_index("author_gene_id")["analysis_symbol"]

    man = pd.read_csv(MANIFEST, sep="\t")
    man["sample_id"] = man["sample_id"].astype(str)

    print("loading TPM", flush=True)
    logx_ens, zeb1_tpm = load_matrix()
    print(f"loaded genes={logx_ens.shape[0]} samples={logx_ens.shape[1]}", flush=True)
    common = [s for s in logx_ens.columns if s in set(man["sample_id"])]
    logx_ens = logx_ens[common]
    zeb1_tpm = zeb1_tpm.loc[common]
    zeb1_log = np.log2(zeb1_tpm + 1.0)
    man = man[man["sample_id"].isin(common)].copy()
    man["ZEB1_log2tpm1"] = man["sample_id"].map(zeb1_log)

    logx = collapse_to_symbols(logx_ens, gene_map)
    rnk_s = spearman_rank(logx, zeb1_log)
    rnk_t, n_low, n_high, med = median_split_t(logx, zeb1_log)

    print("GSEA spearman", flush=True)
    tab_s = run_prerank(rnk_s, "spearman")
    print("GSEA median-split", flush=True)
    tab_t = run_prerank(rnk_t, "median_split")

    def col(tab: pd.DataFrame, *names: str) -> str:
        for n in names:
            if n in tab.columns:
                return n
        raise KeyError(str(list(tab.columns)))

    def pick(tab: pd.DataFrame, name: str) -> dict:
        term_col = col(tab, "Term", "Name")
        row = tab.loc[tab[term_col] == name]
        if row.empty:
            return {"term": name, "found": False}
        r = row.iloc[0]
        fdr_col = col(tab, "FDR q-val", "FDR qval", "fdr")
        p_col = col(tab, "NOM p-val", "NOM pval", "pval")
        size = r.get("Size", r.get("size", np.nan))
        try:
            size = int(size)
        except Exception:
            size = None
        return {
            "term": name,
            "found": True,
            "nes": float(r["NES"]),
            "fdr": float(r[fdr_col]),
            "pval": float(r[p_col]),
            "size": size,
        }

    key_terms = [
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ]

    def summarize(tab: pd.DataFrame) -> dict:
        tab = tab.copy()
        nes_col = "NES"
        fdr_col = "FDR q-val" if "FDR q-val" in tab.columns else "FDR q-val"
        if fdr_col not in tab.columns:
            fdr_col = [c for c in tab.columns if "FDR" in c.upper() or "fdr" in c][0]
        term_col = col(tab, "Term", "Name")
        fdr_col = col(tab, "FDR q-val", "FDR qval", "fdr")
        pos = tab.sort_values(nes_col, ascending=False).head(8)
        neg = tab.sort_values(nes_col, ascending=True).head(8)
        keep = [term_col, nes_col, fdr_col]
        return {
            "key": {k: pick(tab, k) for k in key_terms},
            "top_positive": pos[keep].to_dict("records"),
            "top_negative": neg[keep].to_dict("records"),
            "columns": list(tab.columns),
        }

    subtype = (
        man.groupby("subtype")["ZEB1_log2tpm1"]
        .agg(n="count", median="median", mean="mean")
        .sort_values("median")
        .reset_index()
    )
    subtype.to_csv(OUT / "zeb1_by_subtype.tsv", sep="\t", index=False)

    e2f_s = pick(tab_s, "HALLMARK_E2F_TARGETS")
    e2f_t = pick(tab_t, "HALLMARK_E2F_TARGETS")
    claim_holds = bool(e2f_s.get("found") and e2f_s.get("nes", 0) > 0 and e2f_s.get("fdr", 1) < 0.05)

    receipt = {
        "question": "In Pölönen RNA-seq, is ZEB1-low associated with E2F/proliferation Hallmark programs?",
        "n_patients": int(len(common)),
        "n_genes_ranked_spearman": int(rnk_s.dropna().shape[0]),
        "zeb1_ensembl": ZEB1,
        "zeb1_log2tpm1_median": float(zeb1_log.median()),
        "median_split": {"n_low": n_low, "n_high": n_high, "cutoff_log2tpm1": med},
        "claim_holds_spearman_E2F_FDR0.05": claim_holds,
        "E2F_spearman": e2f_s,
        "E2F_median_split": e2f_t,
        "spearman": summarize(tab_s),
        "median_split_gsea": summarize(tab_t),
        "inputs": {
            "tpm": str(TPM),
            "gene_map": str(GENE_MAP),
            "manifest": str(MANIFEST),
            "gmt": str(GMT),
        },
    }
    (OUT / "RECEIPT.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ["n_patients", "claim_holds_spearman_E2F_FDR0.05", "E2F_spearman", "E2F_median_split"]}, indent=2))


if __name__ == "__main__":
    main()
