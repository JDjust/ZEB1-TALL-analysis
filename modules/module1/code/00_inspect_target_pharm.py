#!/usr/bin/env python3
"""Inspect TARGET lineage, pharmacotype lineage, GPL570 probes, R packages."""
from __future__ import annotations
import json
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/inspect")
OUT.mkdir(parents=True, exist_ok=True)

# TARGET clinical
clin = DATA / "TARGET-ALL-P2/Clinical_Supplement/na"
for p in sorted(clin.glob("*.xlsx")):
    try:
        xl = pd.ExcelFile(p)
        print("FILE", p.name, "sheets", xl.sheet_names[:8])
        df = pd.read_excel(p, sheet_name=0)
        print("  shape", df.shape, "cols", list(df.columns)[:25])
        for c in df.columns:
            cl = str(c).lower()
            if any(k in cl for k in ["lineage", "immunopheno", "all type", "disease", "diagnosis", "t-all", "tall", "subtype", "phenotype", "who"]):
                print("  COL", c, "nunique", df[c].nunique(dropna=False))
                print("   ", df[c].astype(str).value_counts(dropna=False).head(12).to_dict())
    except Exception as e:
        print("ERR", p.name, e)

# RNA file names
rna_dir = DATA / "TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts"
files = sorted(rna_dir.glob("*"))
print("RNA n", len(files))
print("RNA example", files[0].name if files else None)
if files:
    print(files[0].read_text(errors="replace").splitlines()[:6])

# manifest
man = DATA / "TARGET-ALL-P2/gdc_open_manifest.json"
if man.exists():
    obj = json.loads(man.read_text())
    print("manifest type", type(obj).__name__, "len", len(obj) if hasattr(obj, "__len__") else None)
    if isinstance(obj, list):
        print("manifest0 keys", list(obj[0].keys())[:20] if obj else None)
        print("manifest0", {k: obj[0].get(k) for k in list(obj[0])[:15]})
    elif isinstance(obj, dict):
        print("dict keys", list(obj.keys())[:20])

# GPL570 probes
annot = DATA / "GSE13159/GPL570.annot.gz"
print("annot exists", annot.exists(), annot.stat().st_size if annot.exists() else 0)

# pharmacotype lc50
lc = DATA / "StJude_ALL_Pharmacotype/imputed_lc50_data.csv"
df = pd.read_csv(lc, nrows=5)
print("lc50 cols", list(df.columns))
df = pd.read_csv(lc)
print("lc50 shape", df.shape)
for c in df.columns:
    cl = str(c).lower()
    if any(k in cl for k in ["lineage", "immuno", "subtype", "diagnosis", "phenotype", "tall", "ball", "all"]):
        print("LC COL", c, df[c].astype(str).value_counts(dropna=False).head(15).to_dict())
print("lc50 all cols", list(df.columns))

# xlsx supplement
x = DATA / "StJude_ALL_Pharmacotype/41591_2022_2112_MOESM3_ESM.xlsx"
xl = pd.ExcelFile(x)
print("pharm xlsx sheets", xl.sheet_names)
for s in xl.sheet_names[:8]:
    d = pd.read_excel(x, sheet_name=s, nrows=3)
    print(" sheet", s, "cols", list(d.columns)[:20])

# FPKM zip
zpath = DATA / "StJude_ALL_Pharmacotype/pharmacotyping_ped_rnaseq_fpkm.zip"
with zipfile.ZipFile(zpath) as z:
    print("fpkm zip", z.namelist()[:20])

print("DONE")
