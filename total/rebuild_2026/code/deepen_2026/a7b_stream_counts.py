"""A7-B: stream GSE248287 MTX into 15 Dx-malignant patient sums. No classifier."""
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a7b_malignant")

print("load metadata", flush=True)
md = pd.read_csv(P / "GSE248287_10x_metadata_all.txt.gz", sep="\t", low_memory=False)
keep = (
    (md["Type"] == "T-ALL")
    & (md["Timepoint"] == "Dx")
    & (md["Celltypes_all"] == "Malignant")
)
patients = sorted(md.loc[keep, "orig.ident.sample"].unique())
assert len(patients) == 15, patients
p_ix = {p: i for i, p in enumerate(patients)}
col_to_p = np.full(len(md), -1, dtype=np.int32)
col_to_p[np.flatnonzero(keep.to_numpy())] = (
    md.loc[keep, "orig.ident.sample"].map(p_ix).to_numpy(dtype=np.int32)
)
n_keep = int((col_to_p >= 0).sum())
print("dx_malignant_cols", n_keep, "patients", patients, flush=True)

with gzip.open(P / "GSE248287_10x_RNA_features.tsv.gz", "rt") as f:
    genes = [ln.strip() for ln in f]
n_genes = len(genes)
counts = np.zeros((n_genes, 15), dtype=np.int64)

print("stream MTX", flush=True)
n_used = 0
with gzip.open(P / "GSE248287_10x_RNA_matrix.mtx.gz", "rt") as f:
    hdr = f.readline()
    assert hdr.startswith("%%MatrixMarket")
    dims = f.readline()
    while dims.startswith("%"):
        dims = f.readline()
    ng, nc, nnz = [int(x) for x in dims.split()]
    assert ng == n_genes and nc == len(md), (ng, nc, n_genes, len(md))
    print("mtx", ng, nc, nnz, flush=True)
    for k, line in enumerate(f, 1):
        i_s, j_s, x_s = line.split()
        p = col_to_p[int(j_s) - 1]
        if p >= 0:
            counts[int(i_s) - 1, p] += int(x_s)
            n_used += 1
        if k % 20_000_000 == 0:
            print("nnz_read", k, "used", n_used, flush=True)

print("nnz_used", n_used, "libsizes", counts.sum(axis=0).tolist(), flush=True)
lib = pd.DataFrame({
    "patient": patients,
    "n_cells": [int((md.loc[keep, "orig.ident.sample"] == p).sum()) for p in patients],
    "n_umi": counts.sum(axis=0),
})
lib.to_csv(P / "a7b_patient_lib.tsv", sep="\t", index=False)
out = pd.DataFrame(counts, index=genes, columns=patients)
out.index.name = "symbol"
out.to_csv(P / "a7b_patient_counts.tsv.gz", sep="\t")
print("wrote counts", out.shape, flush=True)
