#!/usr/bin/env python3
"""Module 2 rebuild: FACS-sorted postnatal thymus (GSE195812) + pediatric
thymus (GSE206710). ZEB1 remains continuous. Biological n for GSE195812 is
the FACS library (6 donors pooled). Biological n for GSE206710 is donor (n=3).
"""
from __future__ import annotations

import math
import shutil
import tarfile
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import entropy, spearmanr

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2")
M5 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5")
PROC, TAB, FIGP = ROOT / "processed", ROOT / "tables", ROOT / "figures"
for d in (PROC, TAB, FIGP):
    d.mkdir(parents=True, exist_ok=True)

RNG = 1
MAX_PER = 3500
STAGE_ORDER = ["DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP"]
STAGE_INDEX = {s: i for i, s in enumerate(STAGE_ORDER)}
GSE206_META = {
    "GSM6261174": {"sample": "TTA7", "donor": "D1", "fraction": "CD3neg", "rep": "r1", "age_days": 442},
    "GSM6261175": {"sample": "TTA8", "donor": "D1", "fraction": "CD3neg", "rep": "r2", "age_days": 442},
    "GSM6261176": {"sample": "TTA9", "donor": "D1", "fraction": "CD3pos", "rep": "r1", "age_days": 442},
    "GSM6261177": {"sample": "TTA10", "donor": "D1", "fraction": "CD3pos", "rep": "r2", "age_days": 442},
    "GSM6261178": {"sample": "TTA11", "donor": "D2", "fraction": "CD3neg", "rep": "r1", "age_days": 4867},
    "GSM6261179": {"sample": "TTA12", "donor": "D2", "fraction": "CD3pos", "rep": "r1", "age_days": 4867},
    "GSM6261180": {"sample": "TTA13", "donor": "D3", "fraction": "CD3neg", "rep": "r1", "age_days": 190},
    "GSM6261181": {"sample": "TTA14", "donor": "D3", "fraction": "CD3pos", "rep": "r1", "age_days": 190},
}
KEY_GENES = [
    "ZEB1", "ZEB2", "LMO2", "LYL1", "TCF7", "GATA3", "BCL11B", "MEF2C", "IL7R",
    "CD34", "TAL1", "DNTT", "RAG1", "RAG2", "CD1A", "CD1B", "CD3D", "CD3E",
    "CD4", "CD8A", "CD8B", "KIT", "SPI1", "PTCRA", "TRAC", "NOTCH1", "MYC",
    "LEF1", "HES1", "CD7", "CD5", "CD2", "LCK", "SATB1", "RUNX1", "HOXA9",
]
CENTROID_GENES = [
    "CD34", "KIT", "SPI1", "MEF2C", "LMO2", "LYL1", "IL7R", "HES1", "CD1A",
    "CD1B", "RAG1", "DNTT", "CD4", "CD8A", "CD3D", "CD3E", "BCL11B", "TCF7",
    "LEF1", "PTCRA", "TRAC", "GATA3", "TAL1", "NOTCH1", "CD7", "CD5",
]
LINEAGE = {
    "T": ["CD3D", "CD3E", "CD7", "CD2", "TRAC", "LCK", "CD5", "BCL11B", "CD1A"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "PRF1"],
    "B": ["MS4A1", "CD79A", "CD19", "CD74", "HLA-DRA"],
    "Myeloid": ["CD14", "LYZ", "S100A8", "S100A9", "CST3", "FCGR3A"],
    "Ery": ["HBB", "HBA1", "HBA2", "GYPA"],
    "DC": ["IRF8", "IRF7", "CLEC9A", "CD1C", "HLA-DRA"],
}
SIGS = {
    "ETP": ["CD34", "KIT", "IL7R", "LMO2", "LYL1", "MEF2C", "SPI1", "CD44", "HES1"],
    "Tdiff": ["CD1A", "CD3E", "CD4", "CD8A", "RAG1", "BCL11B", "TCF7", "LEF1", "DNTT"],
    "GATA3prog": ["GATA3", "TCF7", "BCL11B", "LEF1", "CD1A"],
    "LMO2prog": ["LMO2", "LYL1", "MEF2C", "SPI1", "CD34", "KIT"],
}


def log(m):
    print(m, flush=True)


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


def gene_vec(adata, gene, use_raw=True):
    src = adata.raw if use_raw and adata.raw is not None else adata
    if gene not in src.var_names:
        return np.full(adata.n_obs, np.nan)
    x = src[:, gene].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def score_list(adata, genes, name):
    present = [g for g in genes if g in (adata.raw.var_names if adata.raw is not None else adata.var_names)]
    if len(present) < 2:
        adata.obs[name] = np.nan
        log(f"score {name}: only {present}")
        return
    sc.tl.score_genes(adata, present, score_name=name, use_raw=True)


def extract_tar(tar_path: Path, dest: Path, done_marker: Path):
    if done_marker.exists():
        log(f"already extracted {tar_path.name}")
        return
    dest.mkdir(parents=True, exist_ok=True)
    log(f"extract {tar_path}")
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(dest)
    done_marker.write_text("ok\n")


def infer_stage(name: str) -> str:
    n = name.upper().replace("-", "_")
    if "DN1" in n:
        return "DN1"
    if "DN2" in n:
        return "DN2"
    if "DN3" in n:
        return "DN3"
    if "ISP" in n:
        return "ISP"
    if "CD3MIN" in n or "CD3NEG" in n or "CD3-" in name.upper():
        return "DP_CD3neg"
    if "CD3PLUS" in n or "CD3POS" in n or "CD3+" in name.upper():
        return "DP_CD3pos"
    if "CD8" in n:
        return "CD8SP"
    if "CD4" in n:
        return "CD4SP"
    return "unknown"


def load_gse195812() -> sc.AnnData:
    h5 = PROC / "GSE195812_qc.h5ad"
    if h5.exists():
        log(f"load cached {h5}")
        return sc.read_h5ad(h5)
    raw_dir = PROC / "gse195812_tenx"
    extract_tar(DATA / "GSE195812/suppl/GSE195812_RAW.tar", raw_dir, raw_dir / "_extracted")
    org = PROC / "gse195812_tenx_org"
    ads = []
    for mtx in sorted(raw_dir.rglob("*_matrix.mtx.gz")):
        prefix = mtx.name.replace("_matrix.mtx.gz", "")
        stage = infer_stage(prefix)
        dest = org / prefix
        dest.mkdir(parents=True, exist_ok=True)
        src = mtx.parent
        for kind, outn in [("matrix.mtx.gz", "matrix.mtx.gz"),
                           ("barcodes.tsv.gz", "barcodes.tsv.gz"),
                           ("features.tsv.gz", "features.tsv.gz")]:
            fp = src / f"{prefix}_{kind}"
            if fp.exists() and not (dest / outn).exists():
                shutil.copy2(fp, dest / outn)
        log(f"read {prefix} stage={stage}")
        a = sc.read_10x_mtx(dest, var_names="gene_symbols", make_unique=True)
        a.var_names_make_unique()
        gsm = prefix.split("_")[0]
        a.obs["sample"] = prefix
        a.obs["gsm"] = gsm
        a.obs["facs_stage"] = stage
        a.obs["atlas"] = "GSE195812"
        a.obs["donor"] = "pooled6"
        a.var["mt"] = a.var_names.str.upper().str.startswith(("MT-", "MT."))
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
        n0 = a.n_obs
        a = a[(a.obs["n_genes_by_counts"] >= 200) &
              (a.obs["n_genes_by_counts"] < 8000) &
              (a.obs["pct_counts_mt"] < 25)].copy()
        log(f"  QC {n0} -> {a.n_obs} median_genes={float(a.obs['n_genes_by_counts'].median()) if a.n_obs else 0:.0f}")
        if a.n_obs == 0:
            log(f"  skip empty {prefix}")
            continue
        if a.n_obs > MAX_PER:
            sc.pp.subsample(a, n_obs=MAX_PER, random_state=RNG)
            log(f"  subsampled to {MAX_PER}")
        ads.append(a)
        log(f"  keep {a.n_obs} cells")
    adata = ad.concat(ads, join="inner", index_unique="-")
    adata.obs_names_make_unique()
    log(f"GSE195812 concat {adata.n_obs} x {adata.n_vars}")
    return process_atlas(adata, batch_key="sample", cache=h5)


def load_gse206710() -> sc.AnnData:
    h5 = PROC / "GSE206710_qc.h5ad"
    if h5.exists():
        log(f"load cached {h5}")
        return sc.read_h5ad(h5)
    raw_dir = PROC / "gse206710_h5"
    extract_tar(DATA / "GSE206710/suppl/GSE206710_RAW.tar", raw_dir, raw_dir / "_extracted")
    ads = []
    for fp in sorted(raw_dir.rglob("*.h5")):
        gsm = fp.name.split("_")[0]
        meta = GSE206_META.get(gsm, {"sample": fp.stem, "donor": "NA", "fraction": "NA", "rep": "r1", "age_days": np.nan})
        log(f"read {fp.name} {meta}")
        a = sc.read_10x_h5(fp)
        a.var_names_make_unique()
        for k, v in meta.items():
            a.obs[k] = v
        a.obs["gsm"] = gsm
        a.obs["atlas"] = "GSE206710"
        a.obs["facs_stage"] = meta["fraction"]
        a.var["mt"] = a.var_names.str.upper().str.startswith(("MT-", "MT."))
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
        n0 = a.n_obs
        a = a[(a.obs["n_genes_by_counts"] >= 200) &
              (a.obs["n_genes_by_counts"] < 8000) &
              (a.obs["pct_counts_mt"] < 25)].copy()
        log(f"  QC {n0} -> {a.n_obs}")
        if a.n_obs == 0:
            continue
        if a.n_obs > MAX_PER:
            sc.pp.subsample(a, n_obs=MAX_PER, random_state=RNG)
            log(f"  subsampled to {MAX_PER}")
        ads.append(a)
        log(f"  keep {a.n_obs} cells")
    adata = ad.concat(ads, join="inner", index_unique="-")
    adata.obs_names_make_unique()
    log(f"GSE206710 concat {adata.n_obs} x {adata.n_vars}")
    return process_atlas(adata, batch_key="sample", cache=h5)


def process_atlas(adata: sc.AnnData, batch_key: str, cache: Path) -> sc.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=15)
    adata = adata[(adata.obs["n_genes_by_counts"] < 8000) & (adata.obs["pct_counts_mt"] < 20)].copy()
    log(f"after QC {adata.n_obs} x {adata.n_vars}")
    try:
        sc.pp.scrublet(adata, batch_key=batch_key)
        n_d = int(adata.obs.get("predicted_doublet", pd.Series(False, index=adata.obs_names)).sum())
        log(f"scrublet doublets {n_d}")
        adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
    except Exception as e:
        log(f"scrublet skipped: {e}")
        adata.obs["predicted_doublet"] = False

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key=batch_key)
    except Exception as e:
        log(f"HVG seurat_v3 failed {e}")
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key=batch_key)
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=50, svd_solver="arpack")
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    use_rep = "X_pca"
    try:
        import omicverse as ov
        ov.single.batch_correction(adata, batch_key=batch_key, methods="harmony", n_pcs=30)
        if "X_harmony" in adata.obsm:
            use_rep = "X_harmony"
        elif "X_pca_harmony" in adata.obsm:
            use_rep = "X_pca_harmony"
        log(f"harmony ok use_rep={use_rep} obsm={list(adata.obsm.keys())}")
    except Exception as e:
        log(f"harmony skipped {e}")
        try:
            sc.external.pp.harmony_integrate(adata, key=batch_key, max_iter_harmony=20)
            use_rep = "X_pca_harmony"
        except Exception as e2:
            log(f"scanpy harmony skipped {e2}")
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=min(30, adata.obsm[use_rep].shape[1]))
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata, resolution=0.6, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    except Exception as e:
        log(f"leiden failed {e}")
        adata.obs["leiden"] = "0"
    for nm, genes in LINEAGE.items():
        score_list(adata, genes, f"score_{nm}")
    for nm, genes in SIGS.items():
        score_list(adata, genes, f"sig_{nm}")
    lin_cols = [f"score_{k}" for k in LINEAGE if f"score_{k}" in adata.obs]
    adata.obs["lineage"] = adata.obs[lin_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    mx = adata.obs[lin_cols].max(axis=1)
    adata.obs.loc[mx < 0.05, "lineage"] = "Unknown"
    for g in KEY_GENES:
        adata.obs[f"g_{g}"] = gene_vec(adata, g)
    try:
        sc.tl.paga(adata, groups="leiden")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sc.pl.paga(adata, color="leiden", show=False, frameon=False)
        tag = "195812" if "GSE195812" in str(cache) else "206710"
        plt.savefig(FIGP / f"M2_r2_paga_{tag}.pdf", bbox_inches="tight")
        plt.savefig(FIGP / f"M2_r2_paga_{tag}.png", dpi=160, bbox_inches="tight")
        plt.close("all")
    except Exception as e:
        log(f"paga skip {e}")
    try:
        sc.tl.diffmap(adata)
        # root = earliest FACS stage among T-lineage
        tmask = adata.obs["lineage"].eq("T")
        if "facs_stage" in adata.obs and tmask.any():
            stages = [s for s in STAGE_ORDER if s in set(adata.obs.loc[tmask, "facs_stage"])]
            root_stage = stages[0] if stages else str(adata.obs.loc[tmask, "facs_stage"].mode().iloc[0])
            idx = np.flatnonzero((adata.obs["facs_stage"].astype(str) == root_stage) & tmask.to_numpy())
            if len(idx):
                adata.uns["iroot"] = int(idx[len(idx) // 2])
                sc.tl.dpt(adata)
                log(f"DPT root {root_stage} n={len(idx)}")
            else:
                adata.obs["dpt_pseudotime"] = np.nan
        else:
            adata.obs["dpt_pseudotime"] = np.nan
    except Exception as e:
        log(f"DPT failed {e}")
        adata.obs["dpt_pseudotime"] = np.nan
    adata.uns["use_rep"] = use_rep
    adata.write_h5ad(cache, compression="gzip")
    log(f"wrote {cache}")
    return adata


def annotate_206_stage(adata: sc.AnnData) -> sc.AnnData:
    """Marker-based developmental bins for pediatric thymus (no FACS DN/DP/SP)."""
    def col(name):
        k = f"g_{name}"
        if k in adata.obs:
            return pd.to_numeric(adata.obs[k], errors="coerce").fillna(0).to_numpy()
        return np.zeros(adata.n_obs)
    cd4, cd8, cd34, cd1a = col("CD4"), col("CD8A"), col("CD34"), col("CD1A")
    stage = np.full(adata.n_obs, "other", dtype=object)
    frac = adata.obs["fraction"].astype(str).to_numpy()
    m = frac == "CD3neg"
    if m.any():
        med34 = float(np.nanmedian(cd34[m]))
        med1a = float(np.nanmedian(cd1a[m]))
        stage[m] = "DN_late_ISP"
        stage[m & (cd34 >= med34)] = "DN_early"
        stage[m & (cd34 < med34) & (cd1a >= med1a)] = "DN_late_ISP"
        stage[m & (cd4 > 0.4) & (cd8 > 0.4)] = "DP_immature"
    m = frac == "CD3pos"
    if m.any():
        p40_8 = float(np.nanpercentile(cd8[m], 40))
        p40_4 = float(np.nanpercentile(cd4[m], 40))
        stage[m] = "DP_mature"
        stage[m & (cd4 > cd8) & (cd8 < p40_8)] = "CD4SP"
        stage[m & (cd8 > cd4) & (cd4 < p40_4)] = "CD8SP"
        stage[m & (cd4 > 0.25) & (cd8 > 0.25)] = "DP_mature"
    adata.obs["inferred_stage"] = stage
    return adata


def export_atlas(adata: sc.AnnData, tag: str, stage_col: str):
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    obs = pd.concat([adata.obs, umap], axis=1)
    rng = np.random.default_rng(RNG)
    parts = []
    for _, idx in obs.groupby("sample", observed=True).groups.items():
        idx = np.asarray(list(idx))
        if len(idx) > 2500:
            idx = rng.choice(idx, 2500, replace=False)
        parts.append(obs.loc[idx])
    plot = pd.concat(parts)
    plot.to_csv(TAB / f"M2_r2_{tag}_plot_cells.tsv.gz", sep="\t", index=False, compression="gzip")
    gcols = [c for c in obs.columns if c.startswith("g_")]
    sig = [c for c in obs.columns if c.startswith("sig_") or c.startswith("score_")]
    keep = ["sample", "gsm", "atlas", "donor", "facs_stage", "lineage", "leiden",
            "n_genes_by_counts", "pct_counts_mt", "dpt_pseudotime", "UMAP1", "UMAP2"]
    if "fraction" in obs.columns:
        keep.append("fraction")
    if "inferred_stage" in obs.columns:
        keep.append("inferred_stage")
    if "age_days" in obs.columns:
        keep.append("age_days")
    obs[[c for c in keep if c in obs.columns] + gcols + sig].to_csv(
        TAB / f"M2_r2_{tag}_all_obs.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    # library / donor means for key genes
    rows = []
    for keys, dd in obs.groupby(["sample", stage_col, "donor"], observed=True):
        rec = {"sample": keys[0], "stage": keys[1], "donor": keys[2], "n_cells": len(dd),
               "atlas": tag, "frac_T": float((dd["lineage"] == "T").mean())}
        for c in gcols + sig + ["dpt_pseudotime"]:
            if c in dd.columns:
                rec[c] = float(pd.to_numeric(dd[c], errors="coerce").mean())
        rows.append(rec)
    means = pd.DataFrame(rows)
    means.to_csv(TAB / f"M2_r2_{tag}_sample_means.tsv", sep="\t", index=False)
    # T-lineage only stage means
    t = obs[obs["lineage"].eq("T")].copy()
    if len(t) and stage_col in t.columns:
        sm = t.groupby(stage_col, observed=True)[gcols].mean()
        sm.to_csv(TAB / f"M2_r2_{tag}_stage_gene_means.tsv", sep="\t")
        ntab = t.groupby(stage_col, observed=True).size().rename("n_cells")
        ntab.to_csv(TAB / f"M2_r2_{tag}_stage_n.tsv", sep="\t")
    log(f"exported {tag} n={len(obs)}")
    return obs, means


def stage_trends(means: pd.DataFrame, tag: str):
    """Spearman of gene vs ordinal stage using sample means (not cells)."""
    d = means.copy()
    d["stage_i"] = d["stage"].map(STAGE_INDEX)
    d = d[d["stage_i"].notna()]
    gcols = [c for c in d.columns if c.startswith("g_")]
    rows = []
    for c in gcols:
        x = pd.to_numeric(d["stage_i"], errors="coerce")
        y = pd.to_numeric(d[c], errors="coerce")
        ok = x.notna() & y.notna()
        if ok.sum() < 5:
            continue
        rho, p = spearmanr(x[ok], y[ok])
        rows.append({"atlas": tag, "gene": c.replace("g_", ""), "rho_vs_stage": rho, "p": p,
                     "n_samples": int(ok.sum()),
                     "min": float(y[ok].min()), "max": float(y[ok].max()),
                     "range": float(y[ok].max() - y[ok].min())})
    out = pd.DataFrame(rows)
    if len(out):
        out["fdr"] = fdr_bh(out["p"].to_numpy())
        out.sort_values("p").to_csv(TAB / f"M2_r2_{tag}_gene_vs_stage.tsv", sep="\t", index=False)
    return out


def centroids_from_ref(obs: pd.DataFrame, stage_col: str):
    t = obs[obs["lineage"].eq("T")].copy()
    genes = [f"g_{g}" for g in CENTROID_GENES if f"g_{g}" in t.columns]
    cent = t.groupby(stage_col, observed=True)[genes].mean()
    cent = cent.reindex([s for s in STAGE_ORDER if s in cent.index])
    zeb = t.groupby(stage_col, observed=True)["g_ZEB1"].mean() if "g_ZEB1" in t.columns else None
    return cent, zeb, genes


def softmax_negdist(dist, tau=0.15):
    z = -dist / max(tau, 1e-6)
    z = z - np.nanmax(z, axis=1, keepdims=True)
    e = np.exp(z)
    e[~np.isfinite(e)] = 0
    s = e.sum(axis=1, keepdims=True)
    s[s == 0] = 1
    return e / s


def map_query(obs_q: pd.DataFrame, cent: pd.DataFrame, zeb_stage: pd.Series, genes: list[str], tag: str):
    present = [g for g in genes if g in obs_q.columns and g in cent.columns]
    if len(present) < 8:
        log(f"mapping {tag}: too few shared genes {present}")
        return pd.DataFrame()
    stages = list(cent.index)
    Xr = cent[present].to_numpy()
    Xr = Xr - np.nanmean(Xr, axis=0, keepdims=True)
    sd = np.nanstd(Xr, axis=0, keepdims=True)
    sd[sd == 0] = 1
    Xr = Xr / sd
    Xq = obs_q[present].to_numpy()
    Xq = (Xq - np.nanmean(cent[present].to_numpy(), axis=0, keepdims=True)) / sd
    Xq = np.nan_to_num(Xq, nan=0.0)
    Xr = np.nan_to_num(Xr, nan=0.0)
    # cosine distance
    nr = np.linalg.norm(Xr, axis=1, keepdims=True)
    nq = np.linalg.norm(Xq, axis=1, keepdims=True)
    nr[nr == 0] = 1
    nq[nq == 0] = 1
    sim = (Xq / nq) @ (Xr / nr).T
    dist = 1 - sim
    prob = softmax_negdist(dist, tau=0.2)
    nearest = [stages[i] for i in dist.argmin(axis=1)]
    H = entropy(np.clip(prob, 1e-12, 1), axis=1) / math.log(len(stages))
    out = obs_q.copy()
    out["nearest_stage"] = nearest
    out["map_entropy"] = H
    out["map_maxp"] = prob.max(axis=1)
    for i, s in enumerate(stages):
        out[f"p_{s}"] = prob[:, i]
    if zeb_stage is not None and "g_ZEB1" in out.columns:
        exp = out["nearest_stage"].map(zeb_stage.to_dict())
        out["ZEB1_expected"] = pd.to_numeric(exp, errors="coerce")
        out["ZEB1_residual"] = pd.to_numeric(out["g_ZEB1"], errors="coerce") - out["ZEB1_expected"]
    out.to_csv(TAB / f"M2_r2_map_{tag}_cells.tsv.gz", sep="\t", index=False, compression="gzip")
    return out


def map_gse227122(cent, zeb_stage, genes):
    h5 = M5 / "processed" / "GSE227122_annotated.h5ad"
    if not h5.exists():
        log("no M5 h5ad; skip T-ALL mapping")
        return
    log(f"map query {h5}")
    q = sc.read_h5ad(h5)
    mal = q.obs["malignant"].astype(str).eq("malignant") if "malignant" in q.obs else pd.Series(True, index=q.obs_names)
    q = q[mal].copy()
    for g in CENTROID_GENES + ["ZEB1", "ZEB2", "LMO2"]:
        col = f"g_{g}"
        if col not in q.obs.columns:
            q.obs[col] = gene_vec(q, g, use_raw=True)
    umap = pd.DataFrame(q.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=q.obs_names) if "X_umap" in q.obsm else None
    obs = q.obs.copy()
    if umap is not None:
        obs = pd.concat([obs, umap], axis=1)
    mapped = map_query(obs, cent, zeb_stage, genes, "GSE227122")
    if not len(mapped):
        return
    # patient-level (Dx)
    if "timepoint" in mapped.columns:
        d = mapped[mapped["timepoint"].astype(str).eq("Dx")].copy()
    else:
        d = mapped
    pcols = [c for c in d.columns if c.startswith("p_")]
    agg = {"n_cells": ("nearest_stage", "size"),
           "g_ZEB1": ("g_ZEB1", "mean"),
           "map_entropy": ("map_entropy", "mean"),
           "ZEB1_residual": ("ZEB1_residual", "mean") if "ZEB1_residual" in d.columns else ("g_ZEB1", "mean")}
    for c in pcols:
        agg[c] = (c, "mean")
    pat = d.groupby("patient", observed=True).agg(**{k: v for k, v in agg.items() if v[0] in d.columns}).reset_index()
    # nearest-stage majority
    maj = d.groupby("patient", observed=True)["nearest_stage"].agg(lambda x: Counter(x).most_common(1)[0][0])
    pat["majority_stage"] = pat["patient"].map(maj)
    pat.to_csv(TAB / "M2_r2_map_GSE227122_patient.tsv", sep="\t", index=False)
    # subtype if present
    if "patient" in d.columns:
        comp = d.groupby(["patient", "nearest_stage"], observed=True).size().reset_index(name="n")
        comp["frac"] = comp["n"] / comp.groupby("patient")["n"].transform("sum")
        comp.to_csv(TAB / "M2_r2_map_GSE227122_stage_comp.tsv", sep="\t", index=False)
    log(f"mapped patients {len(pat)}")


def main():
    log("=== GSE195812 FACS thymus ===")
    a195 = load_gse195812()
    obs195, means195 = export_atlas(a195, "GSE195812", "facs_stage")
    stage_trends(means195, "GSE195812")
    cent, zeb_stage, genes = centroids_from_ref(obs195, "facs_stage")
    cent.to_csv(TAB / "M2_r2_GSE195812_stage_centroids.tsv", sep="\t")
    if zeb_stage is not None:
        zeb_stage.to_csv(TAB / "M2_r2_GSE195812_ZEB1_by_stage.tsv", sep="\t")

    log("=== GSE206710 pediatric thymus ===")
    a206 = load_gse206710()
    a206 = annotate_206_stage(a206)
    obs206, means206 = export_atlas(a206, "GSE206710", "inferred_stage")
    # donor-level CD3neg vs CD3pos ZEB1 (true biological n=3)
    t206 = obs206[obs206["lineage"].eq("T")].copy() if "lineage" in obs206.columns else obs206
    dnr = t206.groupby(["donor", "fraction"], observed=True).agg(
        n=("g_ZEB1", "size"),
        ZEB1=("g_ZEB1", "mean"),
        ZEB2=("g_ZEB2", "mean") if "g_ZEB2" in t206.columns else ("g_ZEB1", "mean"),
        LMO2=("g_LMO2", "mean") if "g_LMO2" in t206.columns else ("g_ZEB1", "mean"),
        TCF7=("g_TCF7", "mean") if "g_TCF7" in t206.columns else ("g_ZEB1", "mean"),
        GATA3=("g_GATA3", "mean") if "g_GATA3" in t206.columns else ("g_ZEB1", "mean"),
    ).reset_index()
    dnr.to_csv(TAB / "M2_r2_GSE206710_donor_fraction.tsv", sep="\t", index=False)

    log("=== map T-ALL GSE227122 onto GSE195812 FACS centroids ===")
    map_gse227122(cent, zeb_stage, genes)

    # replication heatmap: ZEB1 by stage vs CD3 fraction
    note = TAB / "M2_r2_data_note.txt"
    note.write_text(
        "GSE195812: 8 FACS libraries, 6 healthy postnatal donors POOLED per library. "
        "Stage comparisons use the library as the sample (n=8 stages, n_donor=1 pooled). "
        "Cell-level UMAP/distributions are descriptive.\n"
        "GSE206710: 3 pediatric donors (190d, 442d, 4867d), CD3- vs CD3+ (donor1 has 2 reps). "
        "Biological n=3 for fraction tests.\n"
        "RDS All_thymus/Sort1/Sort2 are double-gzipped; RAW MTX used instead of the 13GB RDS.\n"
        "Slingshot/Monocle3/destiny/tradeSeq/SingleR not installed; DPT + FACS order used.\n"
        "T-ALL mapping: cosine to FACS stage centroids on 26 marker genes; patient n from GSE227122 Dx malignant cells.\n"
    )
    log("M2 prepare DONE")


if __name__ == "__main__":
    main()
