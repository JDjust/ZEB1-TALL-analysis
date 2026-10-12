#!/usr/bin/env python3
"""Read-only extract of key genes from existing GSE248287 Dx malignant pseudobulk (run 1908)."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

SRC = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/public_expansion/runs/1908")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/q1_bioinfo_gse248287")
GENES = [
    "ZEB1", "ZEB2", "BCL11B", "MKI67", "E2F1", "CDK1", "CCNE1", "MYC", "PCNA", "TOP2A",
    "SPINK2", "CD34", "CEBPA", "LYL1", "MEF2C", "LMO2", "KIT", "IL7R",
    "CD1A", "CD3E", "CD3D", "CD8A", "RAG1", "TCF7", "LEF1", "GATA3",
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    meta_rows = []
    for meta_path in sorted(SRC.glob("GSE248287_P*_pseudobulk_metadata.tsv"), key=lambda p: int(p.name.split("_")[1][1:])):
        with meta_path.open(encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for rec in r:
                meta_rows.append(rec)
                if rec.get("mapped_state") == "ALL" and rec.get("timepoint") in {"D0", "Dx", "0"}:
                    pass
        patient = meta_path.name.split("_")[1]
        counts_path = SRC / f"GSE248287_{patient}_pseudobulk_counts.tsv.gz"
        with gzip.open(counts_path, "rt", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            all_cols = [i for i, h in enumerate(header) if h.endswith("|ALL")]
            if not all_cols:
                continue
            col = all_cols[0]
            found = {}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                g = parts[0]
                if g in GENES:
                    found[g] = float(parts[col])
        # library size from ALL metadata
        n_cells = lib = None
        with meta_path.open(encoding="utf-8") as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                if rec.get("mapped_state") == "ALL":
                    n_cells = int(float(rec["n_cells"]))
                    lib = float(rec["lib_size"])
                    sample = rec.get("sample")
                    break
        if not lib:
            continue
        row = {"patient": patient, "sample": sample, "n_malignant": n_cells, "library_umi": lib}
        for g in GENES:
            umi = found.get(g)
            row[g + "_umi"] = umi
            row[g] = None if umi is None else (0.0 if lib <= 0 else __import__("math").log2(1.0 + 1e6 * umi / lib))
        rows.append(row)
        print(patient, "n", n_cells, "ZEB1", row.get("ZEB1"), flush=True)
    # write
    if not rows:
        raise SystemExit("no patients")
    keys = list(rows[0].keys())
    with (OUT / "patient_malignant_keygenes.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with (OUT / "pseudobulk_metadata_all.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(meta_rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(meta_rows)
    (OUT / "EXTRACT_RECEIPT.json").write_text(
        json.dumps(
            {
                "source": str(SRC),
                "n_patients": len(rows),
                "genes": GENES,
                "unit": "patient Dx malignant ALL pseudobulk",
                "normalization": "log2(1+CPM) from existing 1908 summed UMI",
                "did_not_modify": str(SRC),
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
