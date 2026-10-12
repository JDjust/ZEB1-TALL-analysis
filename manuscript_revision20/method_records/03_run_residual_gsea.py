#!/usr/bin/env python3
"""Subtype-residual Hallmark GSEA: positive rank = higher with residual ZEB1-high."""
from __future__ import annotations

import json
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd

ROOT = Path(r"d:\_bioinformation\ZEB1")
OUT = ROOT / "02_analyses" / "03_q1_bioinfo_paper"
TPM = ROOT / "00_data" / "00_raw" / "TALL_X01_tpm.tsv"
GENE_MAP = ROOT / ".project-artifacts" / "work" / "20261008T051121Z-384f93b5" / "M1" / "Lim_discovery_GENE_ID_MAP.tsv"
CLIN = ROOT / "00_data" / "00_raw" / "polonen_syn54032669" / "audit_clinical_1309.tsv"
GMT = ROOT / ".project-artifacts" / "work" / "20261006T010250Z-75920166" / "analysis" / "Q2_FUNCTIONAL" / "inputs" / "h.all.v2025.1.Hs.symbols.gmt"
ZEB1 = "ENSG00000148516"
SEED = 20261009


def stable_id(ens: str) -> str:
    return str(ens).split(".")[0]


def load_log_symbols() -> tuple[pd.DataFrame, pd.Series]:
    expr = pd.read_csv(TPM, sep="\t", index_col=0).astype(np.float32)
    expr.index = expr.index.map(stable_id)
    if expr.index.has_duplicates:
        expr = expr.groupby(level=0, sort=False).max()
    zeb1 = np.log2(expr.loc[ZEB1].astype(np.float64) + 1.0)
    logx = np.log2(expr + 1.0)
    logx = logx.loc[logx.mean(axis=1) > 0.5]
    gene_map = pd.read_csv(GENE_MAP, sep="\t")
    gene_map["author_gene_id"] = gene_map["author_gene_id"].map(stable_id)
    symbols = gene_map.drop_duplicates("author_gene_id").set_index("author_gene_id")["analysis_symbol"]
    symbols = symbols.reindex(logx.index).dropna()
    symbols = symbols[symbols.astype(str).str.match(r"^[A-Za-z0-9].*") & ~symbols.astype(str).str.startswith("ENSG")]
    logx = logx.loc[symbols.index]
    var = logx.var(axis=1)
    chosen = symbols.loc[var.sort_values(ascending=False).index].drop_duplicates(keep="first")
    out = logx.loc[chosen.index].copy()
    out.index = chosen.values
    out = out[~out.index.duplicated(keep="first")].drop(index=["ZEB1"], errors="ignore")
    return out, zeb1


def residualize(mat: pd.DataFrame, subtype: pd.Series) -> pd.DataFrame:
    dummies = pd.get_dummies(subtype, drop_first=True).astype(np.float64)
    X = np.column_stack([np.ones(len(subtype)), dummies.to_numpy()])
    xtx = X.T @ X
    xtx_inv = np.linalg.pinv(xtx)
    # Residual = Y - X (X'X)^-1 X' Y
    Y = mat.to_numpy(dtype=np.float64).T  # samples x genes
    beta = xtx_inv @ (X.T @ Y)
    resid = Y - X @ beta
    return pd.DataFrame(resid.T, index=mat.index, columns=mat.columns)


def spearman_pos(logx: pd.DataFrame, z: pd.Series) -> pd.Series:
    zv = z.loc[logx.columns].to_numpy(dtype=np.float64)
    zv = (zv - zv.mean()) / zv.std(ddof=0)
    X = logx.to_numpy(dtype=np.float64)
    X = (X - X.mean(axis=1, keepdims=True)) / np.clip(X.std(axis=1, keepdims=True, ddof=0), 1e-12, None)
    rho = (X @ zv) / X.shape[1]
    return pd.Series(rho, index=logx.index, name="spearman_with_residual_ZEB1")


def pick(tab: pd.DataFrame, name: str) -> dict:
    term_col = "Term" if "Term" in tab.columns else "Name"
    row = tab.loc[tab[term_col] == name]
    if row.empty:
        return {"term": name, "found": False}
    r = row.iloc[0]
    fdr_col = "FDR q-val" if "FDR q-val" in tab.columns else [c for c in tab.columns if "FDR" in c][0]
    return {"term": name, "nes": float(r["NES"]), "fdr": float(r[fdr_col]), "pval": float(r["NOM p-val"])}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("loading", flush=True)
    logx, zeb1 = load_log_symbols()
    clin = pd.read_csv(CLIN, sep="\t")
    clin["sample_id"] = clin["sample_id"].astype(str)
    common = [s for s in logx.columns if s in set(clin["sample_id"])]
    logx = logx[common]
    zeb1 = zeb1.loc[common]
    subtype = clin.set_index("sample_id").loc[common, "subtype"].astype(str)
    print(f"genes={logx.shape[0]} samples={logx.shape[1]} residualizing", flush=True)
    resid_g = residualize(logx, subtype)
    resid_z = residualize(pd.DataFrame([zeb1.to_numpy()], index=["ZEB1"], columns=zeb1.index), subtype).iloc[0]
    rnk = spearman_pos(resid_g, resid_z)
    rnk.to_csv(OUT / "rank_subtype_residual.tsv", sep="\t", header=False)
    print("GSEA", flush=True)
    tab = gp.prerank(
        rnk=rnk.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False),
        gene_sets=str(GMT),
        outdir=str(OUT / "gsea_subtype_residual"),
        permutation_num=1000,
        min_size=15,
        max_size=500,
        seed=SEED,
        threads=1,
        verbose=False,
    ).res2d
    tab.to_csv(OUT / "hallmark_subtype_residual.tsv", sep="\t", index=False)
    keys = [
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "HALLMARK_MYC_TARGETS_V1",
        "HALLMARK_MYC_TARGETS_V2",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ]
    fdr_col = "FDR q-val" if "FDR q-val" in tab.columns else [c for c in tab.columns if "FDR" in c][0]
    sig = tab.loc[tab[fdr_col] < 0.05, ["Term", "NES", fdr_col]].sort_values("NES", ascending=False)
    receipt = {
        "direction": "positive NES = enriched among genes higher with subtype-residual ZEB1-high",
        "n": int(logx.shape[1]),
        "n_genes": int(rnk.dropna().shape[0]),
        "key": {k: pick(tab, k) for k in keys},
        "FDR_lt_0.05": sig.to_dict("records"),
        "E2F_with_ZEB1_high": bool(
            pick(tab, "HALLMARK_E2F_TARGETS").get("nes", 0) > 0
            and pick(tab, "HALLMARK_E2F_TARGETS").get("fdr", 1) < 0.05
        ),
    }
    (OUT / "RESIDUAL_GSEA_RECEIPT.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps({"E2F": receipt["key"]["HALLMARK_E2F_TARGETS"], "G2M": receipt["key"]["HALLMARK_G2M_CHECKPOINT"], "MYC": receipt["key"]["HALLMARK_MYC_TARGETS_V1"], "n_sig": len(sig)}, indent=2))


if __name__ == "__main__":
    main()
