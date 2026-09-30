#!/usr/bin/env python3
"""Inspect TARGET manifest and one RNA file without pandas/openpyxl."""
from __future__ import annotations
import json
from pathlib import Path

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/TARGET-ALL-P2")
man = json.loads((DATA / "gdc_open_manifest.json").read_text())
print("type", type(man).__name__)
if isinstance(man, dict):
    print("keys", list(man.keys())[:30])
    for k, v in list(man.items())[:2]:
        print(" item", k, type(v).__name__, str(v)[:400])
elif isinstance(man, list):
    print("len", len(man))
    print("keys0", list(man[0].keys()) if man else None)
    print("sample0", json.dumps(man[0], indent=2)[:1500])

rna = DATA / "Gene_Expression_Quantification/STAR_-_Counts"
files = sorted(p for p in rna.iterdir() if p.is_file())
print("n_rna_files", len(files))
print("name0", files[0].name)
print("head0:")
print("\n".join(files[0].read_text().splitlines()[:12]))
