#!/usr/bin/env python3
"""ZEB1 log2(1+CPM) by mapped_state from existing 1908 counts. Read-only."""
from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path

SRC = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/public_expansion/runs/1908")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/q1_bioinfo_gse248287")
GENES = ["ZEB1", "BCL11B", "MKI67", "SPINK2", "CD34", "CD1A"]


def main() -> None:
    rows = []
    for meta_path in sorted(SRC.glob("GSE248287_P*_pseudobulk_metadata.tsv"), key=lambda p: int(p.name.split("_")[1][1:])):
        patient = meta_path.name.split("_")[1]
        meta = list(csv.DictReader(meta_path.open(encoding="utf-8"), delimiter="\t"))
        lib = {r["mapped_state"]: float(r["lib_size"]) for r in meta}
        ncell = {r["mapped_state"]: int(float(r["n_cells"])) for r in meta}
        counts_path = SRC / f"GSE248287_{patient}_pseudobulk_counts.tsv.gz"
        with gzip.open(counts_path, "rt", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col_state = {}
            for i, h in enumerate(header[1:], start=1):
                st = h.split("|")[-1]
                col_state[i] = st
            found = {st: {} for st in col_state.values()}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                g = parts[0]
                if g not in GENES:
                    continue
                for i, st in col_state.items():
                    found[st][g] = float(parts[i])
        for st, genes in found.items():
            if st not in lib or lib[st] <= 0:
                continue
            row = {"patient": patient, "mapped_state": st, "n_cells": ncell.get(st), "library_umi": lib[st]}
            for g in GENES:
                umi = genes.get(g, 0.0)
                row[g] = math.log2(1.0 + 1e6 * umi / lib[st])
            rows.append(row)
    keys = list(rows[0].keys())
    with (OUT / "zeb1_by_mapped_state.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("states", len(rows), "patients", len({r['patient'] for r in rows}))


if __name__ == "__main__":
    main()
