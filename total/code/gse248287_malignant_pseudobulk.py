#!/usr/bin/env python3
"""GSE248287 diagnosis-only malignant-cell pseudobulk from author annotations.

Metadata lacks a barcode column. The GEO matrix/barcode/metadata files share row
order; require exact per-cell nCount_RNA agreement before using that alignment.
The unit for correlations is the patient, never the cell.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr

SETS = {
    "ETP": ["CD34", "KIT", "IL7R", "LMO2", "LYL1", "HHEX", "WT1", "BAALC", "IGLL1", "SPINK2", "MEF2C", "SPI1"],
    "Stemness": ["PROM1", "BMI1", "KIT", "SOX4", "THY1", "MSI2", "HOXA9", "MEIS1"],
    "Tcell_diff": ["CD1A", "CD3E", "CD3D", "RAG1", "BCL11B", "TCF7", "LEF1", "CD8A"],
    "GATA3prog": ["GATA3", "TCF7", "BCL11B", "LEF1", "CD1A"],
    "LMO2prog": ["LMO2", "LYL1", "MEF2C", "SPI1", "CD34"],
}
GENES = ["ZEB1", "ZEB2", "LMO2", "IL7R", "GATA3", "TCF7", "CD34", "CD3E"]


def load_sample(folder: Path):
    meta = pd.read_csv(folder / "metadata.tsv.gz", sep="\t")
    features = pd.read_csv(folder / "features.tsv.gz", sep="\t", header=None)[0].astype(str)
    barcodes = pd.read_csv(folder / "barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    with gzip.open(folder / "matrix.mtx.gz", "rb") as fh:
        counts = mmread(fh).tocsr()
    if counts.shape != (len(features), len(barcodes)) or len(meta) != len(barcodes):
        raise ValueError(f"{folder.name}: dimensions do not match")
    if not np.array_equal(np.asarray(counts.sum(axis=0)).ravel().astype(int),
                          meta["nCount_RNA"].to_numpy(dtype=int)):
        raise ValueError(f"{folder.name}: metadata row order does not match matrix columns")
    # The authors' matrices include Dx/Dxh and technical library suffixes
    # such as P8_Dx_2; the metadata independently confirms diagnosis status.
    prefix = "PedTALL_" + folder.name
    if not barcodes.str.match(r"^" + prefix + r"(?:h|_[0-9]+)?_[ACGT]+-[0-9]+$").all():
        raise ValueError(f"{folder.name}: unexpected barcode prefix")
    if not meta["Timepoint"].eq("Dx").all():
        raise ValueError(f"{folder.name}: non-diagnosis cells")
    if features.duplicated().any():
        raise ValueError(f"{folder.name}: duplicate gene symbols; resolve before aggregation")
    return meta, features, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min-cells", type=int, default=30)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    rows, qc = [], []
    for folder in sorted(args.input.glob("P*_Dx"), key=lambda p: int(p.name.split("_")[0][1:])):
        meta, features, counts = load_sample(folder)
        malignant = meta["Celltypes_all"].eq("Malignant").to_numpy()
        n = int(malignant.sum())
        qc.append({"patient": folder.name.split("_")[0], "n_all": len(meta),
                   "n_malignant": n, "n_other": len(meta) - n,
                   "malignant_fraction": n / len(meta), "metadata_alignment": "exact_count_match",
                   "included": n >= args.min_cells})
        if n < args.min_cells:
            continue
        pooled = np.asarray(counts[:, malignant].sum(axis=1)).ravel().astype(float)
        library = float(pooled.sum())
        if library <= 0:
            raise ValueError(f"{folder.name}: zero malignant library")
        index = pd.Index(features)
        cpm = np.log2(1 + 1e6 * pooled / library)
        row = {"patient": folder.name.split("_")[0], "n_malignant": n,
               "library_umi": int(library), "method": "log2(1+CPM) from summed raw malignant UMI"}
        for gene in GENES:
            row[gene] = float(cpm[index.get_loc(gene)]) if gene in index else np.nan
        for name, genes in SETS.items():
            present = [g for g in genes if g in index]
            row[name] = float(np.mean([cpm[index.get_loc(g)] for g in present])) if present else np.nan
            row[name + "_genes"] = ",".join(present)
            row[name + "_n_genes"] = len(present)
        rows.append(row)
    patients = pd.DataFrame(rows)
    qc_table = pd.DataFrame(qc)
    if len(patients) < 6:
        raise ValueError("fewer than six patients pass malignant-cell threshold")
    assoc = []
    features = [g for g in GENES if g != "ZEB1"] + list(SETS)
    for feature in features:
        x, y = patients["ZEB1"], patients[feature]
        ok = x.notna() & y.notna()
        rho, p = spearmanr(x[ok], y[ok]) if ok.sum() >= 6 else (np.nan, np.nan)
        assoc.append({"feature": feature, "n_patients": int(ok.sum()),
                      "rho": float(rho), "p": float(p),
                      "unit": "patient", "cell_selection": "author Celltypes_all=Malignant; Dx"})
    assoc = pd.DataFrame(assoc)
    # Benjamini-Hochberg across all feature correlations in this validation family.
    valid = assoc["p"].notna()
    p = assoc.loc[valid, "p"].to_numpy()
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.clip(q, 0, 1)
    assoc.loc[valid, "fdr_bh"] = out
    patients.to_csv(args.output / "patient_malignant_pseudobulk.tsv", sep="\t", index=False)
    qc_table.to_csv(args.output / "patient_cell_qc.tsv", sep="\t", index=False)
    assoc.to_csv(args.output / "zeb1_patient_associations.tsv", sep="\t", index=False)
    (args.output / "run.json").write_text(json.dumps({
        "input": str(args.input), "min_malignant_cells": args.min_cells,
        "n_samples": len(qc), "n_patients_included": len(patients),
        "cell_label": "author Celltypes_all=Malignant", "timepoint": "Dx",
        "normalization": "sum raw UMI per patient then log2(1+CPM)",
        "multiple_testing": "BH across 12 candidate features",
        "gene_sets": SETS,
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
