#!/usr/bin/env python3
"""5.42: GSE248287 pediatric T-ALL Dx samples from RAW.tar (combined mtx/meta are truncated)."""
from __future__ import annotations

import re
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
RAW = DATA / "GSE248287/suppl/GSE248287_RAW.tar"
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5")
TENX = ROOT / "processed" / "gse248287_dx"
TAB = ROOT / "tables"
TENX.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

SIGS = {
    "ETP": ["CD34", "KIT", "IL7R", "LMO2", "LYL1", "HHEX", "WT1", "BAALC", "IGLL1", "SPINK2", "MEF2C", "SPI1"],
    "Stemness": ["PROM1", "BMI1", "KIT", "SOX4", "THY1", "MSI2", "HOXA9", "MEIS1"],
    "Tcell_diff": ["CD1A", "CD3E", "CD3D", "RAG1", "BCL11B", "TCF7", "LEF1", "CD8A"],
    "GATA3prog": ["GATA3", "TCF7", "BCL11B", "LEF1", "CD1A"],
    "LMO2prog": ["LMO2", "LYL1", "MEF2C", "SPI1", "CD34"],
}


def log(m):
    print(m, flush=True)


def extract_dx():
    want_re = re.compile(r"PedTALL_P\d+_Dx_(RNA_matrix|RNA_features|barcodes|metadata)\.")
    marker = TENX / "P1_Dx" / "matrix.mtx.gz"
    if marker.exists():
        log("248287 Dx already extracted")
        return
    log(f"extract Dx 10x from {RAW}")
    with tarfile.open(RAW, "r") as tf:
        for m in tf.getmembers():
            base = Path(m.name).name
            if not want_re.search(base):
                continue
            mm = re.search(r"(PedTALL_P\d+_Dx)_(RNA_matrix|RNA_features|barcodes|metadata)", base)
            if not mm:
                continue
            sample, kind = mm.group(1), mm.group(2)
            dest = TENX / sample.replace("PedTALL_", "")
            dest.mkdir(parents=True, exist_ok=True)
            if kind == "RNA_matrix":
                out = dest / "matrix.mtx.gz"
            elif kind == "RNA_features":
                out = dest / "features.tsv.gz"
            elif kind == "barcodes":
                out = dest / "barcodes.tsv.gz"
            else:
                out = dest / "metadata.tsv.gz"
            src = tf.extractfile(m)
            if src is None:
                continue
            out.write_bytes(src.read())
            log(f"  {out}")


def load_one(p: Path) -> ad.AnnData:
    feat = pd.read_csv(p / "features.tsv.gz", sep="\t", header=None)
    bcs = pd.read_csv(p / "barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    from scipy.io import mmread
    import gzip
    mat = mmread(gzip.open(p / "matrix.mtx.gz", "rb")).tocsr()
    gnames = (feat[1] if feat.shape[1] >= 2 else feat[0]).astype(str).tolist()
    bcs = bcs.astype(str).tolist()
    if mat.shape == (len(gnames), len(bcs)):
        x = mat.T
    elif mat.shape == (len(bcs), len(gnames)):
        x = mat
    else:
        raise ValueError(f"{p.name} mtx {mat.shape} genes {len(gnames)} bcs {len(bcs)}")
    a = ad.AnnData(X=x.tocsr())
    a.var_names = gnames
    a.obs_names = bcs
    a.var_names_make_unique()
    a.obs_names_make_unique()
    a.obs["sample"] = p.name
    a.obs["patient"] = p.name.split("_")[0]
    a.obs["timepoint"] = "Dx"
    if a.n_obs > 4000:
        sc.pp.subsample(a, n_obs=4000, random_state=1)
    meta = p / "metadata.tsv.gz"
    if meta.exists():
        md = pd.read_csv(meta, sep="\t")
        bc = md.columns[0]
        md["_bc"] = md[bc].astype(str)
        a.obs["_bc"] = a.obs_names.astype(str)
        for c in md.columns:
            cl = str(c).lower()
            if cl in ("celltypes_all", "celltypes_initial", "type") or "celltype" in cl:
                imap = dict(zip(md["_bc"], md[c]))
                a.obs[str(c)] = a.obs["_bc"].map(imap)
    return a


def main():
    extract_dx()
    ads = []
    for p in sorted(TENX.iterdir()):
        if (p / "matrix.mtx.gz").exists():
            log(f"read {p.name}")
            ads.append(load_one(p))
    if not ads:
        raise SystemExit("no 248287 Dx matrices")
    adata = ad.concat(ads, join="inner", index_unique="-")
    log(f"concat {adata.n_obs} x {adata.n_vars} patients={adata.obs.patient.nunique()}")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for g in ["ZEB1", "ZEB2", "LMO2", "IL7R", "GATA3", "TCF7"]:
        if g in adata.var_names:
            vec = adata[:, g].X
            adata.obs[f"g_{g}"] = vec.toarray().ravel() if sparse.issparse(vec) else np.asarray(vec).ravel()
    for nm, genes in SIGS.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"sig_{nm}", use_raw=False)
    gcols = [c for c in adata.obs.columns if c.startswith("g_") or c.startswith("sig_")]
    pb = adata.obs.groupby("patient", observed=True)[gcols].mean()
    pb["n_cells"] = adata.obs.groupby("patient", observed=True).size()
    pb = pb.reset_index()
    pb.to_csv(TAB / "M5_5.42_patient_pseudobulk.tsv", sep="\t", index=False)
    rows = []
    if "g_ZEB1" in pb.columns:
        for c in gcols:
            if c == "g_ZEB1":
                continue
            x = pd.to_numeric(pb["g_ZEB1"], errors="coerce")
            y = pd.to_numeric(pb[c], errors="coerce")
            ok = x.notna() & y.notna()
            if int(ok.sum()) < 6:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            rows.append({"feature": c, "rho": float(rho), "p": float(p), "n": int(ok.sum())})
    pd.DataFrame(rows).to_csv(TAB / "M5_5.42_ZEB1_patient_assoc.tsv", sep="\t", index=False)
    log("GSE248287 Dx patient tables written")


if __name__ == "__main__":
    main()
