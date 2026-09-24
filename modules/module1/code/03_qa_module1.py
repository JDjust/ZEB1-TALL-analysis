#!/usr/bin/env python3
import gzip
from pathlib import Path
import pandas as pd

PROC = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")

mile = pd.read_csv(PROC / "gse13159_keygenes.tsv", sep="\t")
print("MILE ZEB1 range", float(mile.ZEB1.min()), float(mile.ZEB1.max()), "median", float(mile.ZEB1.median()))
print(mile.groupby("disease").ZEB1.agg(["count","median","mean"]))

probes = pd.read_csv(PROC / "gse13159_keygene_probes.tsv", sep="\t", index_col=0)
print("probe shape", probes.shape)
print(probes.describe().T)

sm = DATA / "GSE13159/GSE13159_series_matrix.txt.gz"
want = {"212764_at", "203603_s_at", "204249_s_at", "1007_s_at"}
with gzip.open(sm, "rt", encoding="utf-8", errors="replace") as f:
    in_table = False
    found = 0
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            _ = next(f)
            continue
        if not in_table:
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        pid = line.split("\t", 1)[0].strip().strip('"')
        if pid in want:
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:8]]
            print("RAW", pid, vals)
            found += 1
            if found == 4:
                break

tgt = pd.read_csv(PROC / "target_keygenes.tsv", sep="\t")
print("TARGET n", len(tgt), "usi unique", tgt.usi.nunique())
print(tgt.groupby("disease")[["ZEB1","ZEB2","LMO2"]].median())
print(tgt.sample_type.value_counts())
prim = tgt[tgt.sample_type.astype(str).str.contains("Primary", na=False)]
print("primary n", len(prim), dict(prim.disease.value_counts()), "unique usi", prim.usi.nunique())

ph = pd.read_csv(PROC / "pharmacotype_keygenes.tsv", sep="\t")
print("PHARM ZEB1 range", float(ph.ZEB1.min()), float(ph.ZEB1.max()))
print(ph.groupby("disease")[["ZEB1","ZEB2","LMO2"]].median())
print("pharm disease", dict(ph.disease.value_counts(dropna=False)))
