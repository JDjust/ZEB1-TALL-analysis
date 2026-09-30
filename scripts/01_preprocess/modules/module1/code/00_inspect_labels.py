#!/usr/bin/env python3
"""Inspect labels and gene presence for Module 1 datasets."""
from __future__ import annotations
import gzip
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/inspect")
OUT.mkdir(parents=True, exist_ok=True)

# ---- GSE13159 metadata ----
brief = DATA / "GSE13159" / "GSE13159_gsm_brief.txt"
samples = []
cur = {}
with brief.open(encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gsm": line.split("=", 1)[-1].strip()}
        elif "=" in line and line.startswith("!"):
            k, v = line.split("=", 1)
            k = k.strip()[1:]
            v = v.strip()
            if k == "Sample_characteristics_ch1":
                cur.setdefault("chars", []).append(v)
            else:
                cur[k] = v
    if cur:
        samples.append(cur)

char_keys = Counter()
diag = Counter()
titles = []
for s in samples:
    for c in s.get("chars", []):
        char_keys[c.split(":")[0].strip().lower()] += 1
        cl = c.lower()
        if any(x in cl for x in ["leukemia", "diagnosis", "disease", "phenotype", "cell type", "tissue"]):
            diag[c] += 1
    titles.append(s.get("Sample_title", ""))

title_prefix = Counter()
for t in titles:
    title_prefix[t.split()[0] if t else "NA"] += 1

with (OUT / "gse13159_char_keys.txt").open("w") as f:
    for k, n in char_keys.most_common():
        f.write(f"{n}\t{k}\n")
with (OUT / "gse13159_diag.txt").open("w") as f:
    for k, n in diag.most_common(80):
        f.write(f"{n}\t{k}\n")
with (OUT / "gse13159_title_prefix.txt").open("w") as f:
    for k, n in title_prefix.most_common(40):
        f.write(f"{n}\t{k}\n")

print("GSE13159 n_samples", len(samples))
print("char keys", char_keys.most_common(20))
print("diag top", diag.most_common(15))
print("title prefix", title_prefix.most_common(15))

# series matrix sample titles / characteristics
sm = DATA / "GSE13159" / "GSE13159_series_matrix.txt.gz"
meta_rows = {}
with gzip.open(sm, "rt", encoding="utf-8", errors="replace") as f:
    for i, line in enumerate(f):
        if line.startswith("!"):
            key = line.split("\t", 1)[0]
            if key in {
                "!Sample_title",
                "!Sample_geo_accession",
                "!Sample_source_name_ch1",
                "!Sample_characteristics_ch1",
                "!Sample_characteristics_ch2",
            } or "characteristics" in key.lower() or key in {"!Sample_title"}:
                ntab = line.count("\t")
                meta_rows.setdefault(key, []).append((i, ntab, line[:300]))
        if line.startswith("!series_matrix_table_begin"):
            print("table begin at line", i)
            header = next(f)
            cols = header.rstrip("\n").split("\t")
            print("matrix n_cols", len(cols), "first", cols[:6])
            break
        if i > 200:
            break
print("meta keys", {k: len(v) for k, v in meta_rows.items()})
with (OUT / "gse13159_matrix_meta_keys.txt").open("w") as f:
    for k, vs in meta_rows.items():
        f.write(f"{k}\tn={len(vs)} first={vs[0][2][:200]}\n")

# ---- TARGET clinical ----
clin_dir = DATA / "TARGET-ALL-P2" / "Clinical_Supplement"
clin_files = list(clin_dir.rglob("*"))
print("TARGET clinical files", len(clin_files))
for p in clin_files[:30]:
    if p.is_file():
        print(" CLIN", p.relative_to(clin_dir), p.stat().st_size)

# ---- pharmacotype ----
pharm = DATA / "StJude_ALL_Pharmacotype"
print("pharmacotype files", list(pharm.iterdir()))

print("DONE inspect part1")
