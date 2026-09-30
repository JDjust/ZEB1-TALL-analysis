"""A7-B0: GSE248287 metadata pre-gate. No scores. No classifier."""
from pathlib import Path
import gzip
import pandas as pd

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a7b_malignant")
fp = P / "GSE248287_10x_metadata_all.txt.gz"
out = P / "pregate"
out.mkdir(exist_ok=True)

df = pd.read_csv(fp, sep="\t", low_memory=False)
print("n_rows", len(df), "n_cols", df.shape[1], flush=True)
cols = pd.DataFrame({
    "column": df.columns,
    "dtype": [str(t) for t in df.dtypes],
    "n_unique": [df[c].nunique(dropna=False) for c in df.columns],
    "n_na": [int(df[c].isna().sum()) for c in df.columns],
})
cols.to_csv(out / "a7b0_columns.tsv", sep="\t", index=False)

keys = []
for c in df.columns:
    cl = c.lower()
    if any(k in cl for k in (
        "patient", "sample", "donor", "time", "diagnos", "remiss", "treat",
        "leuk", "malign", "tumor", "blast", "tall", "t-all", "cell", "type",
        "annot", "ident", "cluster", "group", "subtype", "bm", "pb", "blood",
        "marrow", "tissue", "barcode", "hto", "hash",
    )):
        keys.append(c)
(out / "a7b0_priority_cols.txt").write_text("\n".join(keys), encoding="utf-8")

for c in df.columns:
    nu = df[c].nunique(dropna=False)
    if nu <= 80:
        tab = df[c].astype(str).value_counts(dropna=False)
        tab.rename_axis("level").reset_index(name="n").to_csv(
            out / f"a7b0_{''.join(ch if ch.isalnum() else '_' for ch in c)}.tsv",
            sep="\t", index=False,
        )

# first 5 rows, all columns, for field inspection
df.head(5).to_csv(out / "a7b0_head5.tsv", sep="\t", index=False)
print("cols", list(df.columns), flush=True)
print("priority", keys, flush=True)
