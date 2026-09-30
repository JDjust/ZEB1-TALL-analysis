#!/usr/bin/env python3
"""Module 10: GSE287751 Lmo2 KO scRNA (HTO-hashed library).
Biological n is HTO cluster / assigned sample, not cell.
Bulk DN1 KO already showed Lmo2 down, Zeb1 not down.
"""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10")
PROC, TAB = ROOT / "processed", ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

KEYS = ["Zeb1", "Zeb2", "Lmo2", "Tal1", "Lyl1", "Tcf7", "Gata3", "Bcl11b",
        "Il7r", "Cd34", "Dntt", "Rag1", "Kit", "Spi1", "Mef2c", "Cd44"]


def log(m):
    print(m, flush=True)


def extract():
    tar = DATA / "GSE287751/suppl/GSE287751_RAW.tar"
    dest = PROC / "gse287751_scrna"
    dest.mkdir(parents=True, exist_ok=True)
    want = ("cDNA_barcodes", "cDNA_features", "cDNA_matrix", "hto_cluster_umap", "HTO_ft")
    with tarfile.open(tar, "r") as tf:
        for m in tf:
            if not m.isfile():
                continue
            base = Path(m.name).name
            if not any(k.lower() in base.lower() for k in want):
                continue
            out = dest / base
            if out.exists() and out.stat().st_size > 100:
                log(f"have {base}")
                continue
            src = tf.extractfile(m)
            out.write_bytes(src.read())
            log(f"extracted {base} {out.stat().st_size}")
    return dest


def read_xlsx(path: Path) -> pd.DataFrame:
    try:
        return pd.read_excel(path)
    except Exception as e:
        log(f"read_excel failed: {e}")
    try:
        import zipfile
        from xml.etree import ElementTree as ET
        ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        with zipfile.ZipFile(path) as z:
            ss = []
            if "xl/sharedStrings.xml" in z.namelist():
                root = ET.fromstring(z.read("xl/sharedStrings.xml"))
                for si in root.findall("m:si", ns):
                    ss.append("".join(t.text or "" for t in si.iter("{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")))
            sheet = next(n for n in z.namelist() if n.startswith("xl/worksheets/sheet") and n.endswith(".xml"))
            root = ET.fromstring(z.read(sheet))
            rows = []
            for row in root.findall("m:sheetData/m:row", ns):
                vals = []
                for c in row.findall("m:c", ns):
                    t = c.attrib.get("t")
                    v = c.find("m:v", ns)
                    if v is None or v.text is None:
                        vals.append("")
                    elif t == "s":
                        vals.append(ss[int(v.text)] if int(v.text) < len(ss) else v.text)
                    else:
                        vals.append(v.text)
                rows.append(vals)
        if not rows:
            return pd.DataFrame()
        hdr = [str(x) if x else f"c{i}" for i, x in enumerate(rows[0])]
        n = len(hdr)
        body = [r + [""] * (n - len(r)) for r in rows[1:]]
        return pd.DataFrame(body, columns=hdr)
    except Exception as e2:
        log(f"zip xlsx failed: {e2}")
        return pd.DataFrame()


def load_hto(path: Path) -> pd.DataFrame:
    """First line = HTO feature barcodes; later lines = cell barcode + 0/1."""
    rows = []
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split()
        for line in fh:
            p = line.strip().split()
            if len(p) < 2:
                continue
            rec = {"barcode": p[0]}
            nums = [int(x) for x in p[1:]]
            rec["n_pos"] = int(sum(x > 0 for x in nums))
            rec["hto_i"] = int(int(np.argmax(nums)) if nums else -1)
            rec["hto_tag"] = header[rec["hto_i"]] if 0 <= rec["hto_i"] < len(header) else "none"
            rows.append(rec)
    d = pd.DataFrame(rows)
    d["singlet"] = d["n_pos"].eq(1)
    return d, header


def main():
    dest = extract()
    mtx = dest / "GSM9286835_cDNA_matrix.mtx.gz"
    feat = dest / "GSM9286835_cDNA_features.tsv.gz"
    bcs = dest / "GSM9286835_cDNA_barcodes.tsv.gz"
    xlsx = list(dest.glob("*hto_cluster_umap.xlsx")) + list(dest.glob("*assignment*.xlsx"))
    if not mtx.exists():
        raise SystemExit("no mtx")
    genes = pd.read_csv(feat, sep="\t", header=None)
    barcodes = pd.read_csv(bcs, sep="\t", header=None)[0].astype(str)
    log(f"features {len(genes)} barcodes {len(barcodes)} reading mtx")
    mat = mmread(gzip.open(mtx, "rb")).tocsr()
    gnames = (genes[1] if genes.shape[1] >= 2 else genes[0]).astype(str).tolist()
    if mat.shape == (len(gnames), len(barcodes)):
        x = mat.T
    elif mat.shape == (len(barcodes), len(gnames)):
        x = mat
    else:
        raise SystemExit(f"shape {mat.shape} genes {len(gnames)} bcs {len(barcodes)}")
    adata = ad.AnnData(X=x.tocsr())
    adata.var_names = pd.Index(gnames).astype(str)
    adata.obs_names = barcodes.astype(str)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    log(f"AnnData {adata.n_obs} x {adata.n_vars}")
    adata.obs["_bc"] = adata.obs_names.astype(str).str.replace(r"-\d+$", "", regex=True)
    hto_p = dest / "GSM9286836_HTO_ft_count_bin.txt.gz"
    if hto_p.exists():
        hto, header = load_hto(hto_p)
        log(f"HTO cells {len(hto)} tags {header} singlets {int(hto.singlet.sum())}")
        hto.to_csv(TAB / "M10_r2_HTO_assignment_raw.tsv", sep="\t", index=False)
        imap = dict(zip(hto["barcode"].astype(str), hto["hto_tag"]))
        imap_s = dict(zip(hto["barcode"].astype(str), hto["singlet"]))
        adata.obs["hto_tag"] = adata.obs["_bc"].map(imap)
        adata.obs["hto_singlet"] = adata.obs["_bc"].map(imap_s)
        if adata.obs["hto_tag"].isna().mean() > 0.5:
            adata.obs["hto_tag"] = adata.obs_names.astype(str).map(imap)
    if xlsx:
        assign = read_xlsx(xlsx[0])
        log(f"xlsx {xlsx[0].name} {assign.shape} cols={list(assign.columns)[:20]}")
        if len(assign.columns):
            assign.to_csv(TAB / "M10_r2_HTO_xlsx.tsv", sep="\t", index=False)
            bc_col = assign.columns[0]
            assign["_bc"] = assign[bc_col].astype(str)
            for c in assign.columns:
                cl = str(c).lower()
                if any(k in cl for k in ("hto", "hash", "cluster", "sample", "geno", "group",
                                         "umap", "assign", "id", "cond")):
                    adata.obs[str(c)] = adata.obs["_bc"].map(dict(zip(assign["_bc"], assign[c])))
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for g in KEYS:
        if g in adata.var_names:
            vec = adata[:, g].X
            adata.obs[f"g_{g}"] = vec.toarray().ravel() if sparse.issparse(vec) else np.asarray(vec).ravel()
    # find grouping column
    gcol = None
    for c in adata.obs.columns:
        nun = adata.obs[c].nunique(dropna=True)
        cl = str(c).lower()
        if 2 <= nun <= 20 and any(k in cl for k in ("hto", "hash", "cluster", "sample", "geno", "group", "assign")):
            gcol = c
            break
    if gcol is None:
        for c in adata.obs.columns:
            nun = adata.obs[c].nunique(dropna=True)
            if 2 <= nun <= 16 and c not in ("_bc",):
                gcol = c
                break
    log(f"group col={gcol}")
    if "hto_tag" in adata.obs.columns and adata.obs["hto_tag"].notna().sum() >= 50:
        gcol = "hto_tag"
        if "hto_singlet" in adata.obs.columns:
            keep = adata.obs["hto_singlet"].fillna(False).astype(bool)
            if keep.mean() > 0.3:
                adata = adata[keep].copy()
                log(f"kept singlets {adata.n_obs}")
    log(f"group col={gcol}")
    gcols = [c for c in adata.obs.columns if c.startswith("g_")]
    if gcol:
        pb = adata.obs.groupby(gcol, observed=True)[gcols].mean()
        pb["n_cells"] = adata.obs.groupby(gcol, observed=True).size()
        pb = pb.reset_index()
        pb.to_csv(TAB / "M10_r2_HTO_group_means.tsv", sep="\t", index=False)
        # heuristic KO vs WT
        lab = pb[gcol].astype(str)
        pb["geno"] = np.where(lab.str.contains("KO|ko|Lmo2", regex=True), "Lmo2_KO",
                              np.where(lab.str.contains("WT|Ctrl|control|Control", regex=True), "control", "other"))
        pb.to_csv(TAB / "M10_r2_HTO_group_means.tsv", sep="\t", index=False)
        if "g_Zeb1" in adata.obs.columns and "g_Lmo2" in adata.obs.columns:
            rows = []
            for name, sub in adata.obs.groupby(gcol, observed=True):
                rec = {"group": name, "n": len(sub),
                       "Zeb1_mean": float(sub["g_Zeb1"].mean()),
                       "Lmo2_mean": float(sub["g_Lmo2"].mean())}
                rows.append(rec)
            pd.DataFrame(rows).to_csv(TAB / "M10_r2_HTO_Zeb1_Lmo2.tsv", sep="\t", index=False)
        # cell-level MW only as description; n is groups
        note = [f"group_col={gcol}", f"n_groups={adata.obs[gcol].nunique()}",
                f"n_cells={adata.n_obs}"]
        (TAB / "M10_r2_scrna_note.txt").write_text("\n".join(note) + "\n", encoding="utf-8")
    else:
        (TAB / "M10_r2_scrna_note.txt").write_text(
            "HTO assignment could not be joined; xlsx columns did not map to barcodes.\n",
            encoding="utf-8",
        )
        adata.obs[gcols].mean().to_frame("mean").to_csv(TAB / "M10_r2_library_gene_means.tsv", sep="\t")
    # subsample for a UMAP table if grouping worked
    if gcol and adata.n_obs > 4000:
        sc.pp.subsample(adata, n_obs=4000, random_state=1)
    if gcol and adata.n_vars > 50:
        sc.pp.highly_variable_genes(adata, n_top_genes=1500, subset=False)
        sc.pp.pca(adata, n_comps=20)
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)
        plot = adata.obs[[gcol] + [c for c in gcols if c in adata.obs.columns]].copy()
        plot["UMAP1"] = adata.obsm["X_umap"][:, 0]
        plot["UMAP2"] = adata.obsm["X_umap"][:, 1]
        plot.to_csv(TAB / "M10_r2_scrna_plot_cells.tsv.gz", sep="\t", index=False)
    log("M10 scRNA PREPARE DONE")


if __name__ == "__main__":
    main()
