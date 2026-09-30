#!/usr/bin/env python3
import pandas as pd

p = "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/processed/gse287751_scrna/author_cell_metadata.tsv"
full = pd.read_csv(p, sep="\t")
print("columns", list(full.columns))
print("n", len(full))
print(full.head(2).to_string())
for c in full.columns:
    nunq = full[c].nunique(dropna=False)
    if nunq <= 20:
        print(c, full[c].astype(str).value_counts(dropna=False).to_dict())
