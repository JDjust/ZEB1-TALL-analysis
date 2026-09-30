"""A8: merge GSE280250 HTSeq counts. No scores."""
from pathlib import Path
import gzip
import pandas as pd

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a8_adult")
raw = P / "raw_counts"
files = sorted(raw.glob("GSM*.txt.gz"))
assert len(files) == 79, len(files)

cols = {}
genes0 = None
for fp in files:
    gsm = fp.name.split("_", 1)[0]
    g, c = [], []
    with gzip.open(fp, "rt") as fh:
        for ln in fh:
            a, b = ln.rstrip("\n").split("\t")
            g.append(a)
            c.append(int(b))
    if genes0 is None:
        genes0 = g
    else:
        assert g == genes0
    cols[gsm] = c

mat = pd.DataFrame(cols, index=genes0)
mat.index.name = "symbol"
mat.to_csv(P / "a8_counts.tsv.gz", sep="\t")
print("counts", mat.shape, "lib", mat.sum().describe().to_dict(), flush=True)

# series matrix phenotype: GEO repeats !Sample_characteristics_ch1
acc = title = None
age = sex = disease = None
with gzip.open(P / "GSE280250_series_matrix.txt.gz", "rt") as fh:
    for ln in fh:
        if not ln.startswith("!Sample_"):
            continue
        parts = ln.rstrip("\n").split("\t")
        key = parts[0]
        vals = [x.strip().strip('"') for x in parts[1:]]
        if key == "!Sample_geo_accession":
            acc = vals
        elif key == "!Sample_title":
            title = vals
        elif key == "!Sample_characteristics_ch1" and vals:
            tag = vals[0].split(":", 1)[0].strip().lower()
            payload = [x.split(":", 1)[1].strip() if ":" in x else x for x in vals]
            if tag == "age":
                age = [int(x) for x in payload]
            elif tag == "sex":
                sex = payload
            elif tag == "disease":
                disease = payload
assert acc and title and age and sex and disease
ph = pd.DataFrame({
    "geo_accession": acc, "title": title,
    "age": age, "sex": sex, "disease": disease,
})
ph.to_csv(P / "a8_phenotype.tsv", sep="\t", index=False)
print(ph.disease.value_counts().to_string(), flush=True)
print("age_min", ph.age.min(), "age_max", ph.age.max(), "n_lt18", int((ph.age < 18).sum()), flush=True)
