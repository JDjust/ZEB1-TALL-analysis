#!/usr/bin/env python3
"""Inspect scRNA packages and GSE227122 / GSE248287 layout on hulu."""
from __future__ import annotations

import gzip
import importlib
import sys
import tarfile
from pathlib import Path

print("python", sys.version.replace("\n", " "))
for m in [
    "scanpy", "anndata", "harmonypy", "omicverse", "scrublet",
    "infercnvpy", "scvi", "bbknn", "leidenalg", "igraph",
    "matplotlib", "seaborn", "sklearn", "scipy", "pandas", "numpy",
]:
    try:
        mod = importlib.import_module(m)
        print("OK", m, getattr(mod, "__version__", "?"))
    except Exception as e:
        print("NO", m, type(e).__name__)

raw = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE227122/suppl/GSE227122_RAW.tar")
print("RAW", raw.exists(), raw.stat().st_size if raw.exists() else 0)
if raw.exists():
    with tarfile.open(raw, "r") as tf:
        names = tf.getnames()
    print("tar_n", len(names))
    print("tar_head", names[:6])

p248 = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE248287/suppl")
print("248287")
for p in sorted(p248.glob("*")):
    print(" ", p.name, p.stat().st_size)

meta = p248 / "GSE248287_10x_metadata_all.txt.gz"
if meta.exists():
    with gzip.open(meta, "rt") as fh:
        header = fh.readline().rstrip("\n")
        row1 = fh.readline().rstrip("\n")
    print("meta_header", header[:500])
    print("meta_row1", row1[:400])
    n = 0
    with gzip.open(meta, "rt") as fh:
        for n, _ in enumerate(fh, 0):
            pass
    print("meta_lines_incl_header", n)
