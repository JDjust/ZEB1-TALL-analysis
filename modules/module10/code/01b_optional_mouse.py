#!/usr/bin/env python3
"""Optional M10: GSE287751 Lmo2 KO bulk (DN1 n=2 vs 2) and GSE260948
Tal1/Lmo2 mouse thymus scRNA. ZEB1 stays continuous. Do not claim Lmo2-|Zeb1.
"""
from __future__ import annotations

import gzip
import math
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

LOCAL = Path(r"d:\_bioinformation\local_data_staging\revision_tier1")
ROOT = Path(r"d:\_bioinformation\modules\module10")
PROC, TAB, FIG = ROOT / "processed", ROOT / "tables", ROOT / "figures"
for d in (PROC, TAB, FIG, ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

ENS = {
    "ENSMUSG00000024238": "Zeb1", "ENSMUSG00000026872": "Zeb2",
    "ENSMUSG00000032698": "Lmo2", "ENSMUSG00000028717": "Tal1",
    "ENSMUSG00000034041": "Lyl1", "ENSMUSG00000000782": "Tcf7",
    "ENSMUSG00000015619": "Gata3", "ENSMUSG00000048251": "Bcl11b",
    "ENSMUSG00000003882": "Il7r", "ENSMUSG00000016494": "Cd34",
    "ENSMUSG00000028934": "Dntt", "ENSMUSG00000039748": "Rag1",
    "ENSMUSG00000005087": "Cd44", "ENSMUSG00000005672": "Kit",
    "ENSMUSG00000031004": "Mki67", "ENSMUSG00000026923": "Notch1",
    "ENSMUSG00000022346": "Myc", "ENSMUSG00000002111": "Spi1",
    "ENSMUSG00000005583": "Mef2c",
}
KEYS = list(ENS.values())
BULK_META = {
    "GSM8750806": {"stage": "DN1", "geno": "control", "rep": "r1"},
    "GSM8750807": {"stage": "DN1", "geno": "control", "rep": "r2"},
    "GSM8750808": {"stage": "DN1", "geno": "Lmo2_KO", "rep": "r1"},
    "GSM8750809": {"stage": "DN1", "geno": "Lmo2_KO", "rep": "r2"},
    "GSM8750810": {"stage": "DN2", "geno": "Lmo2_KO", "rep": "r1"},
    "GSM8750811": {"stage": "DN2", "geno": "Lmo2_KO", "rep": "r2"},
}
SC_META = {
    "GSM8128573": {"state": "WT", "genotype": "WT", "age": "6w"},
    "GSM8128574": {"state": "preleukemic", "genotype": "Tal1_Lmo2", "age": "6w"},
    "GSM8128575": {"state": "TALL1", "genotype": "Tal1_Lmo2", "age": "3mo"},
    "GSM8128576": {"state": "TALL2", "genotype": "Tal1_Lmo2", "age": "3mo"},
}
MAX_PER = 2500
MIN_UMI = 500


def log(m):
    print(m, flush=True)


def hedges_g(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    sp = math.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if not sp:
        return np.nan
    return float(((a.mean() - b.mean()) / sp) * (1 - 3 / (4 * (na + nb) - 9)))


def extract_member(tar_path: Path, dest_dir: Path, pred):
    dest_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r") as tf:
        for m in tf.getmembers():
            name = Path(m.name).name
            if not pred(name):
                continue
            out = dest_dir / name
            if out.exists() and out.stat().st_size > 100:
                continue
            out.write_bytes(tf.extractfile(m).read())
            log(f"extracted {name}")


def bulk_gse287751():
    tar = LOCAL / "GSE287751/suppl/GSE287751_RAW.tar"
    dest = PROC / "gse287751_rsem"
    extract_member(tar, dest, lambda n: n.endswith("genes.results.gz"))
    recs = []
    for fp in sorted(dest.glob("GSM875*.genes.results.gz")):
        gsm = fp.name.split("_")[0]
        d = pd.read_csv(fp, sep="\t")
        d["ens"] = d["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        d["symbol"] = d["ens"].map(ENS)
        sub = d[d["symbol"].notna()][["ens", "symbol", "TPM", "expected_count"]].copy()
        meta = BULK_META[gsm]
        sub["gsm"] = gsm
        for k, v in meta.items():
            sub[k] = v
        recs.append(sub)
        log(f"{fp.name} key_genes={len(sub)}")
    long = pd.concat(recs, ignore_index=True)
    long.to_csv(TAB / "M10_r2_GSE287751_key_long.tsv", sep="\t", index=False)
    wide = long.pivot_table(index="symbol", columns="gsm", values="TPM")
    wide.to_csv(TAB / "M10_r2_GSE287751_key_TPM.tsv", sep="\t")
    rows = []
    dn1 = long[long.stage.eq("DN1")]
    for gene, g in dn1.groupby("symbol"):
        a = g.loc[g.geno == "Lmo2_KO", "TPM"].to_numpy()
        b = g.loc[g.geno == "control", "TPM"].to_numpy()
        rows.append({
            "gene": gene, "stage": "DN1",
            "TPM_KO_r1_r2": ",".join(f"{x:.3f}" for x in a),
            "TPM_ctrl_r1_r2": ",".join(f"{x:.3f}" for x in b),
            "mean_KO": float(np.mean(a)), "mean_control": float(np.mean(b)),
            "log2FC_KO_vs_ctrl": math.log2((np.mean(a) + 0.1) / (np.mean(b) + 0.1)),
            "hedges_g": hedges_g(a, b), "n_KO": len(a), "n_ctrl": len(b),
        })
    stats = pd.DataFrame(rows).sort_values("gene")
    stats.to_csv(TAB / "M10_r2_GSE287751_DN1_KO_vs_ctrl.tsv", sep="\t", index=False)
    log(stats.to_string(index=False))
    (TAB / "M10_r2_GSE287751_note.txt").write_text(
        "GSE287751 bulk RSEM: DN1 control n=2, DN1 Lmo2 KO n=2, DN2 Lmo2 KO n=2 "
        "(no DN2 control). scRNA 235 MB mtx not analyzed. n=2 cannot support causality.\n"
    )
    return stats


def _stream_keep_barcodes(mtx_gz: Path, n_bc: int) -> np.ndarray:
    colsum = np.zeros(n_bc, dtype=np.float64)
    with gzip.open(mtx_gz, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit() and not hasattr(_stream_keep_barcodes, "_hdr"):
                # first non-comment is dims
                _stream_keep_barcodes._hdr = True
                continue
            if len(parts) < 3:
                continue
            b = int(parts[1]) - 1
            if 0 <= b < n_bc:
                colsum[b] += float(parts[2])
    _stream_keep_barcodes._hdr = False
    keep = np.flatnonzero((colsum >= MIN_UMI) & (colsum < 40000))
    log(f"  UMI filter {n_bc} -> {len(keep)} (min_UMI={MIN_UMI})")
    if len(keep) > MAX_PER:
        rng = np.random.default_rng(1)
        keep = np.sort(rng.choice(keep, MAX_PER, replace=False))
        log(f"  subsample {MAX_PER}")
    return keep


def _write_filtered_10x(mtx_gz, bc_gz, ft_gz, out_dir: Path, keep_idx: np.ndarray):
    out_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(bc_gz, "rt") as fh:
        bcs = [ln.rstrip() for ln in fh]
    keep_set = set(int(i) for i in keep_idx)
    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as fo:
        for i in keep_idx:
            fo.write(bcs[int(i)] + "\n")
    shutil.copyfile(ft_gz, out_dir / "features.tsv.gz")
    old_to_new = {int(old): new + 1 for new, old in enumerate(keep_idx)}
    entries = []
    n_genes = None
    with gzip.open(mtx_gz, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            g, b, n = line.split()[:3]
            if n_genes is None:
                n_genes = int(g)
                continue
            bi = int(b) - 1
            if bi in old_to_new:
                entries.append((int(g), old_to_new[bi], n))
    with gzip.open(out_dir / "matrix.mtx.gz", "wt") as fo:
        fo.write("%%MatrixMarket matrix coordinate real general\n")
        fo.write(f"{n_genes} {len(keep_idx)} {len(entries)}\n")
        for g, b, n in entries:
            fo.write(f"{g} {b} {n}\n")
    log(f"  wrote filtered 10x n={len(keep_idx)} nnz={len(entries)}")


def scrna_gse260948():
    import scanpy as sc
    tar = LOCAL / "GSE260948/suppl/GSE260948_RAW.tar"
    dest = PROC / "gse260948_10x"
    extract_member(tar, dest, lambda n: n.endswith((".mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")))
    ads = []
    for gsm, meta in SC_META.items():
        mtx = next(dest.glob(f"{gsm}_*_matrix.mtx.gz"))
        bcs = next(dest.glob(f"{gsm}_*_barcodes.tsv.gz"))
        fts = next(dest.glob(f"{gsm}_*_features.tsv.gz"))
        n_bc = sum(1 for _ in gzip.open(bcs, "rt"))
        log(f"{gsm} {meta} barcodes={n_bc}")
        lib = dest / f"{gsm}_filt"
        if n_bc > 30000:
            keep = _stream_keep_barcodes(mtx, n_bc)
            _write_filtered_10x(mtx, bcs, fts, lib, keep)
        else:
            lib.mkdir(exist_ok=True)
            for src, name in ((mtx, "matrix.mtx.gz"), (bcs, "barcodes.tsv.gz"), (fts, "features.tsv.gz")):
                tgt = lib / name
                if not tgt.exists():
                    shutil.copyfile(src, tgt)
        a = sc.read_10x_mtx(lib, var_names="gene_symbols", cache=False)
        a.var_names_make_unique()
        for k, v in meta.items():
            a.obs[k] = v
        a.obs["gsm"] = gsm
        a.var["mt"] = a.var_names.str.upper().str.startswith(("MT-", "MT."))
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        n0 = a.n_obs
        a = a[(a.obs["n_genes_by_counts"] >= 200) &
              (a.obs["n_genes_by_counts"] < 8000) &
              (a.obs["pct_counts_mt"] < 20)].copy()
        log(f"  QC {n0} -> {a.n_obs}")
        if a.n_obs > MAX_PER:
            sc.pp.subsample(a, n_obs=MAX_PER, random_state=1)
        ads.append(a)
    import anndata as ad
    adata = ad.concat(ads, join="inner", index_unique="-")
    adata.obs_names_make_unique()
    log(f"concat {adata.n_obs} x {adata.n_vars}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for g in KEYS:
        if g in adata.var_names:
            x = adata[:, g].X
            adata.obs[f"g_{g}"] = np.asarray(x.todense()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()
    obs = adata.obs.copy()
    gene_cols = [c for c in obs.columns if c.startswith("g_")]
    obs.to_csv(TAB / "M10_r2_GSE260948_obs.tsv.gz", sep="\t")
    libm = obs.groupby(["gsm", "state", "genotype"], observed=True)[gene_cols].mean().reset_index()
    ntab = obs.groupby("gsm").size().rename("n")
    libm = libm.merge(ntab, on="gsm")
    libm.to_csv(TAB / "M10_r2_GSE260948_library_means.tsv", sep="\t", index=False)
    log(libm.to_string(index=False))
    (TAB / "M10_r2_GSE260948_note.txt").write_text(
        "GSE260948: WT / Tal1-Lmo2 preleukemic / T-ALL#1 / T-ALL#2. "
        f"TALL2 raw had ~1.8e6 barcodes (unfiltered); UMI>={MIN_UMI} then subsample {MAX_PER}. "
        "Biological n: 1 WT, 1 PL, 2 T-ALL libraries. No trajectory.\n"
    )


def main():
    log("=== GSE287751 Lmo2 KO bulk ===")
    bulk_gse287751()
    log("=== GSE260948 Tal1/Lmo2 thymus scRNA ===")
    scrna_gse260948()
    log("OPTIONAL M10 DONE")


if __name__ == "__main__":
    main()
