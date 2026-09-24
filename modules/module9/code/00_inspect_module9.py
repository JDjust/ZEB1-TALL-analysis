#!/usr/bin/env python3
"""Inspect M9 datasets, R packages, and gene annotation on hulu."""
from __future__ import annotations

import gzip
import json
import subprocess
import tarfile
from pathlib import Path

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
AN = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
ACC = [
    "GSE110635", "GSE110632", "GSE62143", "GSE62141",
    "GSE186943", "GSE144035", "GSE154675", "GSE188225",
    "GSE287751", "GSE110636",
]


def head_gz(p: Path, n=8):
    try:
        with gzip.open(p, "rt", errors="replace") as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                print(f"    {line.rstrip()[:200]}")
    except Exception as e:
        print("    HEAD_FAIL", e)


def sample_titles(matrix: Path, n=40):
    if not matrix.exists():
        print("  no matrix", matrix)
        return
    opener = gzip.open if str(matrix).endswith(".gz") else open
    print("  MATRIX", matrix.name, matrix.stat().st_size)
    with opener(matrix, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession") or line.startswith("!Sample_characteristics") or line.startswith("!Sample_source"):
                print("   ", line.rstrip()[:400])


print("=== DATASETS ===")
for acc in ACC:
    d = DATA / acc
    print(f"\n## {acc} exists={d.exists()}")
    if not d.exists():
        continue
    files = sorted(p for p in d.rglob("*") if p.is_file())
    print(f"  n_files={len(files)} bytes={sum(p.stat().st_size for p in files)}")
    for p in files:
        print(f"  {p.relative_to(d)} {p.stat().st_size}")
    for m in sorted((d / "matrix").glob("*")) if (d / "matrix").exists() else []:
        sample_titles(m)
    for extra in [
        d / "suppl" / f"{acc}_FPKM_and_read_count_matrix_data_file.csv.gz",
        d / "suppl" / f"{acc}_All_counts_data.txt.gz",
        d / "suppl" / f"{acc}_Lmo2_transgenic_thymocyte_CellPop_GenewiseCounts.txt.gz",
    ]:
        if extra.exists():
            print("  COUNTS_HEAD", extra.name)
            head_gz(extra, 6)
    tar = d / "suppl" / f"{acc}_RAW.tar"
    if tar.exists() and tar.stat().st_size < 50_000_000:
        print("  TAR_MEMBERS")
        with tarfile.open(tar, "r") as tf:
            for m in tf.getmembers()[:15]:
                print(f"    {m.name} {m.size}")

print("\n=== MODULE PROCESSED ===")
for m in ["module3", "module4"]:
    p = AN / m / "processed"
    print(m, p.exists())
    if p.exists():
        for f in sorted(p.glob("*"))[:30]:
            print(" ", f.name, f.stat().st_size)

print("\n=== R PKGS ===")
rs = "/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript"
code = r"""
pkgs <- c("DESeq2","limma","edgeR","fgsea","ggplot2","data.table","ggrepel",
          "ComplexHeatmap","circlize","patchwork","RobustRankAggreg","RRHO2","RRHO",
          "biomaRt","org.Hs.eg.db","org.Mm.eg.db","AnnotationDbi","PCAtools",
          "decoupleR","msigdbr","apeglm","ashr","pheatmap","igraph")
for (p in pkgs) cat(p, requireNamespace(p, quietly=TRUE), "\n")
"""
r = subprocess.run([rs, "-e", code], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-800:] if r.stderr else "")

print("\n=== PY PKGS ===")
for pkg in ["pandas", "numpy", "pyBigWig", "scipy", "statsmodels"]:
    try:
        m = __import__(pkg)
        print(pkg, getattr(m, "__version__", "ok"))
    except Exception as e:
        print(pkg, "NO", e)
