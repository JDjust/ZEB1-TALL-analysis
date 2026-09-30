#!/usr/bin/env python3
"""Score ZEB1 CollecTRI activity on the Gate 1 DN3 states already saved.

Uses the local human CollecTRI table. Targets are matched to mouse symbols by
the standard first-letter capitalization. Activity is a weighted mean of
gene-wise z-scores. The cell-cycle control drops S/G2M genes from the regulon
and also residualizes activity on S and G2M scores inside each library.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/dormant_lic/gate1")
TENX = OUT / "10x"
NET = Path("/data-b/liangfuhua/projects/bioinfo_direction_screening_20260903/config/collectri_human.tsv")
LIBS = {"TALL1": TENX / "TALL1_filt", "TALL2": TENX / "TALL2_filt"}


def log(msg: str) -> None:
    print(msg, flush=True)


def mouse_symbol(gene: str) -> str:
    g = str(gene)
    return g[:1].upper() + g[1:].lower()


def dense_genes(adata, genes):
    x = adata[:, genes].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=float)


def activity(adata, genes, weights) -> np.ndarray:
    x = dense_genes(adata, genes)
    z = (x - x.mean(0)) / (x.std(0) + 1e-8)
    return z @ weights / np.sum(np.abs(weights))


def load_regulon(var_names: set[str]) -> pd.DataFrame:
    net = pd.read_csv(NET, sep="\t")
    log(f"collectri columns={list(net.columns)} rows={len(net)}")
    src = "source_genesymbol" if "source_genesymbol" in net.columns else "source"
    tgt = "target_genesymbol" if "target_genesymbol" in net.columns else "target"
    stim = net["consensus_stimulation"].astype(str).str.lower().isin({"true", "1", "1.0"})
    inhib = net["consensus_inhibition"].astype(str).str.lower().isin({"true", "1", "1.0"})
    net = net.copy()
    net["weight"] = np.where(stim & ~inhib, 1.0, np.where(inhib & ~stim, -1.0, np.nan))
    zeb = net[net[src].astype(str).str.upper().eq("ZEB1") & net["weight"].notna()][[tgt, "weight"]].copy()
    zeb.columns = ["target_human", "weight"]
    zeb["target"] = zeb["target_human"].map(mouse_symbol)
    zeb = zeb[zeb["target"].isin(var_names)].drop_duplicates("target")
    log(f"ZEB1 targets in matrix={len(zeb)} of human ZEB1 rows={int(net[src].astype(str).str.upper().eq('ZEB1').sum())}")
    return zeb


def qc_norm(path: Path):
    adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[
        (adata.obs.n_genes_by_counts >= 200)
        & (adata.obs.n_genes_by_counts < 8000)
        & (adata.obs.pct_counts_mt < 20)
    ].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def contrast(df: pd.DataFrame, col: str, library: str, layer: str) -> dict:
    a = df.loc[df.state.eq("dDN3"), col].astype(float)
    b = df.loc[df.state.eq("pDN3"), col].astype(float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    return {
        "library": library,
        "layer": layer,
        "feature": col,
        "n_dDN3": int(a.size),
        "n_pDN3": int(b.size),
        "mean_dDN3": float(a.mean()) if len(a) else np.nan,
        "mean_pDN3": float(b.mean()) if len(b) else np.nan,
        "delta_d_minus_p": float(a.mean() - b.mean()) if len(a) and len(b) else np.nan,
    }


def residualize(df: pd.DataFrame, col: str) -> pd.Series:
    use = df[["S_score", "G2M_score", col]].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
    y = use[col].to_numpy()
    x = np.column_stack([np.ones(len(use)), use["S_score"].to_numpy(), use["G2M_score"].to_numpy()])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    out = pd.Series(np.nan, index=df.index)
    out.loc[use.index] = y - x @ beta
    return out


def main() -> None:
    saved = {}
    genes = None
    for name, path in LIBS.items():
        obs = pd.read_csv(OUT / f"{name}_dn3_obs.tsv.gz", sep="\t", index_col=0)
        adata = qc_norm(path)
        common = obs.index.intersection(adata.obs_names)
        log(f"{name} saved={len(obs)} matched={len(common)}")
        adata = adata[common].copy()
        adata.obs = obs.loc[common].copy()
        saved[name] = adata
        genes = set(adata.var_names) if genes is None else genes & set(adata.var_names)
    reg = load_regulon(genes)
    reg.to_csv(OUT / "zeb1_collectri_targets_in_matrix.tsv", sep="\t", index=False)
    cycle = set()
    # Cycle genes were stored only in uns of the first run; recover from the score columns' known lists.
    from gse260948_gate1 import S_GENES, G2M_GENES
    cycle = set(S_GENES) | set(G2M_GENES)
    rows = []
    for name, adata in saved.items():
        targets = reg["target"].tolist()
        weights = reg["weight"].astype(float).to_numpy()
        adata.obs["ZEB1_activity"] = activity(adata, targets, weights)
        keep = ~reg["target"].isin(cycle)
        adata.obs["ZEB1_activity_nocycle"] = activity(
            adata, reg.loc[keep, "target"].tolist(), reg.loc[keep, "weight"].astype(float).to_numpy()
        )
        adata.obs["ZEB1_activity_resid"] = residualize(adata.obs, "ZEB1_activity")
        adata.obs["ZEB1_activity_nocycle_resid"] = residualize(adata.obs, "ZEB1_activity_nocycle")
        for col, layer in [
            ("ZEB1_activity", "ZEB1"),
            ("ZEB1_activity_nocycle", "ZEB1_nocycle_targets"),
            ("ZEB1_activity_resid", "ZEB1_residualized_on_S_G2M"),
            ("ZEB1_activity_nocycle_resid", "ZEB1_nocycle_and_residualized"),
        ]:
            rows.append(contrast(adata.obs, col, name, layer))
        adata.obs.to_csv(OUT / f"{name}_dn3_obs.tsv.gz", sep="\t")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "zeb1_activity_contrasts.tsv", sep="\t", index=False)
    print(tab.to_string(index=False), flush=True)
    raw = tab[tab.layer.eq("ZEB1")]
    controlled = tab[tab.layer.eq("ZEB1_nocycle_and_residualized")]
    def both_neg(frame):
        signs = np.sign(frame["delta_d_minus_p"].to_numpy())
        return bool(len(signs) == 2 and set(signs) == {-1.0})
    decision = {
        "activity_both_negative": both_neg(raw),
        "activity_controlled_both_negative": both_neg(controlled),
        "n_targets": int(len(reg)),
        "n_targets_nocycle": int((~reg["target"].isin(cycle)).sum()),
    }
    decision["gate1_activity"] = "PASS" if decision["activity_both_negative"] and decision["activity_controlled_both_negative"] else "FAIL"
    pd.Series(decision).to_csv(OUT / "gate1_activity_decision.tsv", sep="\t", header=False)
    log(str(decision))


if __name__ == "__main__":
    main()
