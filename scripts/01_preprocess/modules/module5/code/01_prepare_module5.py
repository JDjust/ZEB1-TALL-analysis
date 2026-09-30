#!/usr/bin/env python3
"""Module 5 prepare: GSE227122 scRNA (16 samples / 11 patients) + GSE248287 CITE-seq.

ZEB1 stays continuous. Tests that claim a state use patient as the biological n,
not cell. No MTX / LMO2 causal story is imposed.
"""
from __future__ import annotations

import gzip
import math
import re
import shutil
import tarfile
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import entropy, spearmanr, wilcoxon

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5")
PROC = ROOT / "processed"
TAB = ROOT / "tables"
FIGP = ROOT / "figures"
for d in (PROC, TAB, FIGP):
    d.mkdir(parents=True, exist_ok=True)

RAW_TAR = DATA / "GSE227122/suppl/GSE227122_RAW.tar"
TENX = PROC / "tenx"
H5AD = PROC / "GSE227122_annotated.h5ad"
MAX_PER_SAMPLE = 8000
RNG = 1

SIGS = {
    "Tcell_diff": ["CD1A", "CD1B", "CD3E", "CD3D", "CD4", "CD8A", "RAG1", "RAG2", "DNTT", "BCL11B", "TCF7", "LEF1", "CD5", "CD2", "LCK"],
    "ETP": ["CD34", "KIT", "IL7R", "LMO2", "LYL1", "HHEX", "WT1", "BAALC", "IGLL1", "SPINK2", "IGFBP7", "MEF2C", "SPI1", "CD44", "HES1"],
    "Stemness": ["PROM1", "BMI1", "KIT", "ABCG2", "SOX4", "MYCN", "ALDH1A1", "THY1", "NES", "MSI2", "HOXA9", "MEIS1"],
    "Prolif": ["MKI67", "PCNA", "TOP2A", "CDK1", "CCNB1", "CCNA2", "MCM2", "MCM7", "BIRC5", "UBE2C"],
    "Apoptosis": ["BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "CASP3", "CASP8", "FAS", "PMAIP1", "BBC3"],
    "TGFb": ["TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SERPINE1", "ID1", "SKIL"],
    "JAKSTAT": ["IL7R", "JAK1", "JAK3", "STAT5A", "STAT5B", "STAT1", "IL2RG", "SOCS1", "CISH"],
    "NOTCH": ["NOTCH1", "NOTCH3", "RBPJ", "HES1", "HEY1", "DTX1", "MYC", "NRARP"],
    "PI3K": ["AKT1", "MTOR", "RPS6KB1", "EIF4EBP1", "PIK3CA", "PTEN", "GSK3B", "TSC2"],
    "GATA3prog": ["GATA3", "TCF7", "BCL11B", "LEF1", "CD1A", "RAG1"],
    "LMO2prog": ["LMO2", "LYL1", "MEF2C", "SPI1", "CD34", "KIT"],
}
LINEAGE = {
    "T": ["CD3D", "CD3E", "CD7", "CD2", "TRAC", "LCK", "CD5"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "PRF1"],
    "B": ["MS4A1", "CD79A", "CD19", "CD74", "HLA-DRA"],
    "Myeloid": ["CD14", "LYZ", "S100A8", "S100A9", "CST3", "FCGR3A", "CSF1R"],
    "Ery": ["HBB", "HBA1", "HBA2", "GYPA"],
}
KEY_GENES = ["ZEB1", "ZEB2", "LMO2", "LYL1", "IL7R", "CD34", "CD7", "CD3D", "GATA3", "TCF7", "BCL11B", "MEF2C", "SPI1", "RAG1", "MKI67"]


def log(m):
    print(m, flush=True)


def parse_sample(name: str):
    m = re.match(r"(T\d+)(Dx|EOI|Rel)$", name)
    if not m:
        raise ValueError(name)
    return m.group(1), m.group(2)


def extract_tenx():
    marker = TENX / "T1Dx" / "matrix.mtx.gz"
    if marker.exists():
        log("tenx already extracted")
        return
    TENX.mkdir(parents=True, exist_ok=True)
    tmp = PROC / "_tar_extract"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    log(f"extract {RAW_TAR}")
    with tarfile.open(RAW_TAR, "r") as tf:
        tf.extractall(tmp)
    n = 0
    for f in tmp.rglob("*"):
        if not f.is_file():
            continue
        m = re.match(r"GSM\d+_(T\d+(?:Dx|EOI|Rel))_(barcodes|features|matrix)\.(.+)$", f.name)
        if not m:
            continue
        sample, kind, _rest = m.group(1), m.group(2), m.group(3)
        dest = TENX / sample
        dest.mkdir(parents=True, exist_ok=True)
        if kind == "matrix":
            out = dest / "matrix.mtx.gz"
        elif kind == "barcodes":
            out = dest / "barcodes.tsv.gz"
        else:
            out = dest / "features.tsv.gz"
        shutil.copy2(f, out)
        n += 1
    shutil.rmtree(tmp)
    log(f"organized {n} 10x files into {TENX}")


def load_gse227122() -> sc.AnnData:
    ads = []
    for sample_dir in sorted(p for p in TENX.iterdir() if p.is_dir()):
        mtx = sample_dir / "matrix.mtx.gz"
        if not mtx.exists():
            continue
        log(f"read {sample_dir.name}")
        a = sc.read_10x_mtx(sample_dir, var_names="gene_symbols", make_unique=True)
        a.var_names_make_unique()
        patient, tp = parse_sample(sample_dir.name)
        a.obs["sample"] = sample_dir.name
        a.obs["patient"] = patient
        a.obs["timepoint"] = tp
        a.obs["gsm"] = sample_dir.name
        if a.n_obs > MAX_PER_SAMPLE:
            sc.pp.subsample(a, n_obs=MAX_PER_SAMPLE, random_state=RNG)
            log(f"  subsampled to {MAX_PER_SAMPLE}")
        ads.append(a)
        log(f"  {a.n_obs} cells {a.n_vars} genes")
    adata = ad.concat(ads, join="inner", label="concat_batch", index_unique="-")
    adata.obs_names_make_unique()
    log(f"concat {adata.n_obs} x {adata.n_vars}")
    return adata


def score_list(adata, genes, name, use_raw=True):
    present = [g for g in genes if g in (adata.raw.var_names if use_raw and adata.raw is not None else adata.var_names)]
    if len(present) < 2:
        adata.obs[name] = np.nan
        log(f"score {name}: only {present}")
        return
    sc.tl.score_genes(adata, present, score_name=name, use_raw=use_raw)


def patient_entropy(labels, groups):
    out = {}
    for g in pd.unique(groups):
        sub = labels[groups == g]
        if len(sub) < 5:
            out[g] = np.nan
            continue
        c = np.array(list(Counter(sub).values()), dtype=float)
        out[g] = float(entropy(c / c.sum(), base=2))
    return out


def run_gse227122() -> sc.AnnData:
    if H5AD.exists():
        log(f"load cached {H5AD}")
        return sc.read_h5ad(H5AD)

    extract_tenx()
    adata = load_gse227122()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)

    qc = adata.obs.groupby("sample", observed=True).agg(
        n_cells=("sample", "size"),
        median_nFeature=("n_genes_by_counts", "median"),
        median_nCount=("total_counts", "median"),
        median_mito=("pct_counts_mt", "median"),
        patient=("patient", "first"),
        timepoint=("timepoint", "first"),
    ).reset_index()
    qc.to_csv(TAB / "M5_5.1_sample_qc.tsv", sep="\t", index=False)

    # keep a QC snapshot before filtering for 5.2-5.4
    adata.obs[["sample", "patient", "timepoint", "n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(
        TAB / "M5_qc_cells_prefilter.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=20)
    adata = adata[(adata.obs["n_genes_by_counts"] < 8000) & (adata.obs["pct_counts_mt"] < 20)].copy()
    log(f"after QC filter {adata.n_obs} x {adata.n_vars}")

    try:
        sc.pp.scrublet(adata, batch_key="sample")
        n_d = int(adata.obs.get("predicted_doublet", pd.Series(False, index=adata.obs_names)).sum())
        log(f"scrublet doublets {n_d}")
        adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
    except Exception as e:
        log(f"scrublet skipped: {e}")
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_score"] = np.nan

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", layer="counts", batch_key="sample")
    except Exception as e:
        log(f"seurat_v3 HVG failed ({e}); seurat flavor")
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="sample")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns["pca"]
    adata.varm["PCs"] = np.zeros((adata.n_vars, adata_hvg.obsm["X_pca"].shape[1]))

    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    adata.obsm["X_umap_unint"] = np.array(adata.obsm["X_umap"], copy=True)

    use_rep = "X_pca"
    try:
        import omicverse as ov
        ov.single.batch_correction(adata, batch_key="sample", methods="harmony", n_pcs=30)
        if "X_harmony" in adata.obsm:
            use_rep = "X_harmony"
            log("omicverse Harmony returned coordinates; convergence not independently verified")
        elif "X_pca_harmony" in adata.obsm:
            use_rep = "X_pca_harmony"
            log("omicverse pca_harmony OK")
        else:
            log("omicverse batch_correction ran; obsm=" + ",".join(adata.obsm.keys()))
            if "X_emb" in adata.obsm:
                use_rep = "X_emb"
    except Exception as e:
        log(f"harmony skipped: {e}")
        try:
            sc.external.pp.harmony_integrate(adata, key="sample", max_iter_harmony=30)
            use_rep = "X_pca_harmony"
            log("scanpy harmony OK")
        except Exception as e2:
            log(f"scanpy harmony skipped: {e2}")

    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=min(30, adata.obsm[use_rep].shape[1]))
    sc.tl.umap(adata)
    try:
        sc.tl.leiden(adata, resolution=0.6, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
        sc.tl.leiden(adata, resolution=1.0, key_added="leiden_hi", flavor="igraph", n_iterations=2, directed=False)
        log("leiden igraph OK")
    except Exception as e:
        log(f"leiden igraph failed: {e}")
        conn = adata.obsp["connectivities"]
        import igraph as ig
        sources, targets = conn.nonzero()
        g = ig.Graph(n=conn.shape[0], edges=list(zip(sources.tolist(), targets.tolist())), directed=False)
        g.es["weight"] = conn[sources, targets].A1
        part = g.community_multilevel(weights="weight")
        adata.obs["leiden"] = pd.Categorical([str(x) for x in part.membership])
        adata.obs["leiden_hi"] = adata.obs["leiden"]
        log("igraph community_multilevel OK")
    log("clusters " + str(adata.obs["leiden"].value_counts().to_dict()))

    for nm, genes in LINEAGE.items():
        score_list(adata, genes, f"score_{nm}")
    for nm, genes in SIGS.items():
        score_list(adata, genes, f"sig_{nm}")

    lin_cols = [f"score_{k}" for k in LINEAGE]
    adata.obs["lineage"] = adata.obs[lin_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    mx = adata.obs[lin_cols].max(axis=1)
    adata.obs.loc[mx < 0.05, "lineage"] = "Unknown"

    # malignancy: T-lineage clusters that are patient-private (low patient entropy)
    ent = patient_entropy(adata.obs["patient"].to_numpy(), adata.obs["leiden"].to_numpy())
    adata.obs["cluster_patient_entropy"] = adata.obs["leiden"].map(ent).astype(float)
    cl_t = adata.obs.groupby("leiden", observed=True)["score_T"].median()
    cl_my = adata.obs.groupby("leiden", observed=True)["score_Myeloid"].median()
    cl_b = adata.obs.groupby("leiden", observed=True)["score_B"].median()
    mal_cl = []
    for cl, e in ent.items():
        t = float(cl_t.get(cl, np.nan))
        my = float(cl_my.get(cl, np.nan))
        b = float(cl_b.get(cl, np.nan))
        if (t > my) and (t > b) and (np.isfinite(e) and e < 1.6 or t > 0.4):
            mal_cl.append(str(cl))
    adata.obs["malignant"] = np.where(
        adata.obs["leiden"].astype(str).isin(mal_cl) & adata.obs["lineage"].isin(["T", "Unknown"]),
        "malignant",
        "nonmalignant",
    )
    # force T-high cells in private clusters
    log("malignant clusters " + ",".join(mal_cl) + f"  frac={adata.obs.malignant.eq('malignant').mean():.3f}")

    # CNV-like proxy: 9p21 (CDKN2A/B often deleted) and T-ALL blast vs myeloid
    score_list(adata, ["CDKN2A", "CDKN2B", "MTAP", "IFNA1"], "cnv_9p21")
    adata.obs["cnv_proxy"] = adata.obs["score_T"] - adata.obs["score_Myeloid"]

    # CytoTRACE substitute: gene detection among malignant, rank inverted so high = less differentiated
    adata.obs["cytotrace_proxy"] = adata.obs["n_genes_by_counts"].astype(float)
    mal = adata.obs["malignant"].eq("malignant")
    if mal.any():
        x = adata.obs.loc[mal, "cytotrace_proxy"]
        adata.obs.loc[mal, "cytotrace_proxy"] = (x - x.mean()) / (x.std(ddof=0) + 1e-8)

    # PAGA + DPT on all cells, root = malignant cluster with max ETP
    sc.tl.paga(adata, groups="leiden")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sc.pl.paga(adata, color="leiden", show=False, frameon=False, title="PAGA leiden")
        plt.savefig(FIGP / "M5_5.23_paga.pdf", bbox_inches="tight")
        plt.savefig(FIGP / "M5_5.23_paga.png", dpi=180, bbox_inches="tight")
        plt.close("all")
    except Exception as e:
        log(f"paga plot skip {e}")

    sc.tl.diffmap(adata)
    etp_by_cl = adata.obs.loc[mal].groupby("leiden", observed=True)["sig_ETP"].median() if mal.any() else pd.Series(dtype=float)
    if len(etp_by_cl):
        root_cl = str(etp_by_cl.idxmax())
        idx = np.flatnonzero((adata.obs["leiden"].astype(str) == root_cl).to_numpy())
        adata.uns["iroot"] = int(idx[len(idx) // 2]) if len(idx) else 0
        adata.uns["dpt_root_cluster"] = root_cl
        log(f"DPT root cluster {root_cl} n={len(idx)}")
        try:
            sc.tl.dpt(adata)
        except Exception as e:
            log(f"dpt failed {e}")
            adata.obs["dpt_pseudotime"] = np.nan
    else:
        adata.obs["dpt_pseudotime"] = np.nan

    # gene values from raw
    raw = adata.raw.to_adata()
    for g in KEY_GENES:
        if g in raw.var_names:
            vec = raw[:, g].X
            if sparse.issparse(vec):
                vec = vec.toarray().ravel()
            else:
                vec = np.asarray(vec).ravel()
            adata.obs[f"g_{g}"] = vec
        else:
            adata.obs[f"g_{g}"] = np.nan
            log(f"missing gene {g}")

    adata.write_h5ad(H5AD, compression="gzip")
    log(f"wrote {H5AD}")
    return adata


def subsample_plot_obs(adata, n_per=2500):
    parts = []
    rng = np.random.default_rng(RNG)
    for _s, idx in adata.obs.groupby("sample", observed=True).groups.items():
        idx = np.asarray(list(idx))
        if len(idx) > n_per:
            idx = rng.choice(idx, n_per, replace=False)
        parts.append(adata.obs.loc[idx])
    return pd.concat(parts, axis=0)


def export_tables(adata):
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    u2 = pd.DataFrame(adata.obsm.get("X_umap_unint", adata.obsm["X_umap"]), columns=["UMAP1_unint", "UMAP2_unint"], index=adata.obs_names)
    obs = pd.concat([adata.obs, umap, u2], axis=1)
    keep = [c for c in obs.columns if c not in ("concat_batch",)]
    plot = subsample_plot_obs(adata)
    plot = plot.join(umap).join(u2)
    plot.to_csv(TAB / "M5_plot_cells.tsv.gz", sep="\t", index=False, compression="gzip")
    obs[keep].to_csv(TAB / "M5_all_obs.tsv.gz", sep="\t", index=False, compression="gzip")

    # composition
    comp = obs.groupby(["sample", "patient", "timepoint", "lineage"], observed=True).size().reset_index(name="n")
    comp["frac"] = comp["n"] / comp.groupby("sample")["n"].transform("sum")
    comp.to_csv(TAB / "M5_5.13_composition.tsv", sep="\t", index=False)

    malc = obs.groupby(["sample", "patient", "timepoint", "malignant"], observed=True).size().reset_index(name="n")
    malc["frac"] = malc["n"] / malc.groupby("sample")["n"].transform("sum")
    malc.to_csv(TAB / "M5_malignant_composition.tsv", sep="\t", index=False)

    # marker means by lineage
    gcols = [c for c in obs.columns if c.startswith("g_")]
    if gcols:
        mm = obs.groupby("lineage", observed=True)[gcols].mean().T
        mm.to_csv(TAB / "M5_5.9_marker_means.tsv", sep="\t")

    # patient-level stats on Dx malignant cells (biological n)
    dx = obs[(obs["timepoint"] == "Dx") & (obs["malignant"] == "malignant")].copy()
    sig_cols = [c for c in obs.columns if c.startswith("sig_")]
    gene_cols = [c for c in gcols]
    pat = dx.groupby("patient", observed=True).agg(
        n_cells=("sample", "size"),
        **{c: (c, "mean") for c in sig_cols + gene_cols + ["cytotrace_proxy", "dpt_pseudotime", "cnv_proxy"]}
    ).reset_index()
    # ZEB1-low fraction: within-patient bottom quartile would be tautological.
    # Use ETP-high fraction among malignant (score > cohort median).
    if "sig_ETP" in dx.columns and dx["sig_ETP"].notna().any():
        thr = float(dx["sig_ETP"].median())
        etp_hi = dx.assign(hi=dx["sig_ETP"] > thr).groupby("patient", observed=True)["hi"].mean().rename("frac_ETP_high")
        pat = pat.merge(etp_hi, on="patient", how="left")
    if "g_ZEB1" in dx.columns:
        zmed = float(dx["g_ZEB1"].median())
        zlo = dx.assign(lo=dx["g_ZEB1"] < zmed).groupby("patient", observed=True)["lo"].mean().rename("frac_ZEB1_below_global_median")
        pat = pat.merge(zlo, on="patient", how="left")
    pat.to_csv(TAB / "M5_5.36_patient_pseudobulk.tsv", sep="\t", index=False)

    # Spearman at patient n
    rows = []
    if "g_ZEB1" in pat.columns:
        for c in [x for x in pat.columns if x not in ("patient", "n_cells", "g_ZEB1")]:
            x = pd.to_numeric(pat["g_ZEB1"], errors="coerce")
            y = pd.to_numeric(pat[c], errors="coerce")
            ok = x.notna() & y.notna()
            if ok.sum() < 6:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            rows.append({"feature": c, "rho": rho, "p": p, "n": int(ok.sum())})
    pdf = pd.DataFrame(rows)
    if len(pdf):
        pvals = pdf["p"].to_numpy()
        n = len(pvals)
        order = np.argsort(pvals)
        ranked = np.empty(n, dtype=float)
        prev = 1.0
        for i, k in enumerate(order[::-1], start=0):
            rank = n - i
            val = min(prev, pvals[k] * n / rank)
            ranked[k] = val
            prev = val
        pdf["fdr"] = ranked
        pdf.sort_values("p").to_csv(TAB / "M5_5.41_patient_ZEB1_associations.tsv", sep="\t", index=False)

    # paired Dx vs EOI, patient means of malignant cells
    both = obs[obs["malignant"] == "malignant"].copy()
    pair_rows = []
    for patient, dd in both.groupby("patient", observed=True):
        tps = set(dd["timepoint"])
        if not {"Dx", "EOI"} <= tps:
            continue
        rec = {"patient": patient}
        for col in ["g_ZEB1", "g_ZEB2", "g_LMO2", "g_IL7R", "sig_ETP", "sig_Stemness", "sig_Tcell_diff", "dpt_pseudotime"]:
            if col not in dd.columns:
                continue
            rec[f"{col}_Dx"] = dd.loc[dd.timepoint == "Dx", col].mean()
            rec[f"{col}_EOI"] = dd.loc[dd.timepoint == "EOI", col].mean()
        rec["n_Dx"] = int((dd.timepoint == "Dx").sum())
        rec["n_EOI"] = int((dd.timepoint == "EOI").sum())
        pair_rows.append(rec)
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(TAB / "M5_5.37_Dx_EOI_paired.tsv", sep="\t", index=False)
    if len(pairs) >= 3 and "g_ZEB1_Dx" in pairs.columns:
        with open(TAB / "M5_5.37_paired_tests.tsv", "w") as fh:
            fh.write("feature\tn\tW\tp\n")
            for col in ["g_ZEB1", "g_LMO2", "sig_ETP", "sig_Stemness"]:
                a = pd.to_numeric(pairs[f"{col}_Dx"], errors="coerce")
                b = pd.to_numeric(pairs[f"{col}_EOI"], errors="coerce")
                ok = a.notna() & b.notna()
                if ok.sum() < 3:
                    continue
                try:
                    w, p = wilcoxon(a[ok], b[ok], zero_method="wilcox", alternative="two-sided")
                except Exception:
                    w, p = np.nan, np.nan
                fh.write(f"{col}\t{int(ok.sum())}\t{w}\t{p}\n")

    # relapse observation
    rel = obs[obs["timepoint"] == "Rel"]
    rel.groupby(["patient", "malignant", "lineage"], observed=True).size().reset_index(name="n").to_csv(
        TAB / "M5_5.39_relapse_counts.tsv", sep="\t", index=False
    )

    # cluster summary
    cl = obs.groupby("leiden", observed=True).agg(
        n=("leiden", "size"),
        n_patients=("patient", "nunique"),
        entropy=("cluster_patient_entropy", "first"),
        frac_malignant=("malignant", lambda s: np.mean(s == "malignant")),
        ZEB1=("g_ZEB1", "mean") if "g_ZEB1" in obs.columns else ("leiden", "size"),
        ETP=("sig_ETP", "mean") if "sig_ETP" in obs.columns else ("leiden", "size"),
        lineage_mode=("lineage", lambda s: s.value_counts().index[0]),
    ).reset_index()
    cl.to_csv(TAB / "M5_cluster_summary.tsv", sep="\t", index=False)
    log("exported tables")


def run_gse248287():
    """5.42 second scRNA/CITE-seq validation at patient level."""
    suppl = DATA / "GSE248287/suppl"
    meta_f = suppl / "GSE248287_10x_metadata_all.txt.gz"
    mtx_f = suppl / "GSE248287_10x_RNA_matrix.mtx.gz"
    feat_f = suppl / "GSE248287_10x_RNA_features.tsv.gz"
    bc_f = suppl / "GSE248287_10x_barcodes.tsv.gz"
    if not (meta_f.exists() and mtx_f.exists()):
        log("GSE248287 files missing")
        return
    meta = pd.read_csv(meta_f, sep="\t")
    log(f"GSE248287 meta {meta.shape} cols={list(meta.columns)[:25]}")
    meta.to_csv(TAB / "M5_5.42_metadata_preview.tsv", sep="\t", index=False)
    # identify likely columns
    cols = {c.lower(): c for c in meta.columns}

    def pick(*cands):
        for c in cands:
            if c.lower() in cols:
                return cols[c.lower()]
        for k, v in cols.items():
            for c in cands:
                if c.lower() in k:
                    return v
        return None

    c_type = pick("celltype", "cell_type", "annotation", "cluster")
    c_pat = pick("patient", "donor", "sample", "orig.ident", "orig_ident")
    c_diag = pick("diagnosis", "disease", "status", "tissue", "group")
    log(f"248287 mapped type={c_type} patient={c_pat} diag={c_diag}")

    # load matrix only if we can map barcodes; subsample later
    try:
        from scipy.io import mmread
        log("read GSE248287 mtx")
        mat = mmread(gzip.open(mtx_f, "rb")).tocsr()
        if mat.shape[0] < mat.shape[1]:
            # genes x cells typical
            pass
        genes = pd.read_csv(feat_f, sep="\t", header=None)
        bcs = pd.read_csv(bc_f, sep="\t", header=None)[0].astype(str)
        if genes.shape[1] >= 2:
            gnames = genes[1].astype(str)
        else:
            gnames = genes[0].astype(str)
        # mtx is genes x cells
        if mat.shape[0] != len(gnames) and mat.shape[1] == len(gnames):
            mat = mat.T
        adata = ad.AnnData(X=mat.T if mat.shape[0] == len(gnames) else mat)
        if adata.n_vars != len(gnames):
            adata = ad.AnnData(X=mat.T)
        adata.var_names = pd.Index(gnames).make_unique()
        n_cells = adata.n_obs
        if len(bcs) == n_cells:
            adata.obs_names = pd.Index(bcs).make_unique()
        log(f"248287 adata {adata.n_obs} x {adata.n_vars}")
    except Exception as e:
        log(f"248287 matrix load failed: {e}")
        return

    # merge meta
    m = meta.copy()
    bc_col = pick("barcode", "cell", "cell_id", "barcode_id")
    if bc_col is None:
        # first column
        bc_col = m.columns[0]
    m["_bc"] = m[bc_col].astype(str).str.replace("-1$", "", regex=True)
    adata.obs["_bc"] = adata.obs_names.astype(str).str.replace("-1$", "", regex=True)
    # try several join keys
    joined = False
    for left, right in [
        (adata.obs_names.astype(str), m[bc_col].astype(str)),
        (adata.obs["_bc"], m["_bc"]),
    ]:
        imap = dict(zip(right, range(len(right))))
        hit = sum(x in imap for x in left[: min(5000, len(left))])
        log(f"meta match probe {hit}/5000 using {bc_col}")
        if hit > 100:
            idx = [imap.get(x, -1) for x in left]
            keep = np.array(idx) >= 0
            adata = adata[keep].copy()
            idx = np.array(idx)[keep]
            for col in m.columns:
                adata.obs[col] = m.iloc[idx][col].to_numpy()
            joined = True
            break
    if not joined:
        log("could not join GSE248287 metadata; skip expression scoring")
        m.head(20).to_csv(TAB / "M5_5.42_unjoined_meta_head.tsv", sep="\t", index=False)
        return

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for g in ["ZEB1", "ZEB2", "LMO2", "IL7R"]:
        if g in adata.var_names:
            vec = adata[:, g].X
            adata.obs[f"g_{g}"] = vec.toarray().ravel() if sparse.issparse(vec) else np.asarray(vec).ravel()
    for nm, genes in SIGS.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"sig_{nm}", use_raw=False)

    pat_col = c_pat if c_pat in adata.obs.columns else None
    if pat_col is None:
        log("no patient column in joined 248287")
        adata.obs.head(0)
        return
    # prefer T-ALL malignant / tumor if a type column exists
    sub = adata.obs.copy()
    if c_type and c_type in sub.columns:
        t = sub[c_type].astype(str).str.lower()
        mask = t.str.contains("t-all|blast|malignant|leukemia|t cell|tcell|cd7")
        if mask.sum() >= 500:
            sub = sub[mask]
            log(f"248287 subset {c_type} n={len(sub)}")
    gcols = [c for c in sub.columns if c.startswith("g_") or c.startswith("sig_")]
    pb = sub.groupby(pat_col, observed=True)[gcols].mean()
    pb["n_cells"] = sub.groupby(pat_col, observed=True).size()
    pb = pb.reset_index()
    pb.to_csv(TAB / "M5_5.42_patient_pseudobulk.tsv", sep="\t", index=False)
    if "g_ZEB1" in pb.columns and len(pb) >= 6:
        rows = []
        for c in gcols:
            if c == "g_ZEB1":
                continue
            x = pd.to_numeric(pb["g_ZEB1"], errors="coerce")
            y = pd.to_numeric(pb[c], errors="coerce")
            ok = x.notna() & y.notna()
            if ok.sum() < 6:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            rows.append({"feature": c, "rho": rho, "p": p, "n": int(ok.sum())})
        pd.DataFrame(rows).to_csv(TAB / "M5_5.42_ZEB1_patient_assoc.tsv", sep="\t", index=False)
    log("GSE248287 patient table done")


def main():
    sc.settings.verbosity = 2
    sc.settings.figdir = str(FIGP)
    adata = run_gse227122()
    export_tables(adata)
    try:
        run_gse248287()
    except Exception as e:
        log(f"GSE248287 failed: {e}")
    log("MODULE5 PREPARE DONE")
    log(f"cells={adata.n_obs} genes={adata.n_vars} samples={adata.obs['sample'].nunique()} patients={adata.obs['patient'].nunique()}")


if __name__ == "__main__":
    main()
