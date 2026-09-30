#!/usr/bin/env python3
"""Patient-level raw-UMI aggregation for a fixed four-gene scRNA check.

GSE227122 uses the prior project's heuristic candidate-malignant label.
GSE248287 uses the original authors' malignant annotation. They are distinct
evidence tiers and are never pooled as cells.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from gse248287_malignant_pseudobulk import load_sample

GENES = ("ZEB1", "GATA3", "TCF7", "BCL11B")


def row(cohort: str, patient: str, counts: np.ndarray, index: pd.Index,
        n_cells: int, label: str) -> dict:
    library = float(np.sum(counts))
    if n_cells < 30 or library <= 0:
        raise ValueError(f"{cohort}/{patient}: insufficient cells or library")
    result = dict(cohort=cohort, patient=patient, n_cells=n_cells,
                  library_umi=int(library), malignant_label_source=label,
                  expression_unit="log2(1+CPM), summed raw UMI")
    for gene in GENES:
        if gene not in index:
            raise ValueError(f"{cohort}: {gene} absent")
        value = float(counts[index.get_loc(gene)])
        result[f"{gene}_count"] = int(value)
        result[gene] = float(np.log2(1 + value * 1e6 / library))
    return result


def gse227122(path: Path) -> list[dict]:
    a = ad.read_h5ad(path)
    if a.var_names.has_duplicates or "counts" not in a.layers:
        raise ValueError("GSE227122 raw count schema changed")
    obs = a.obs
    mask = obs.timepoint.eq("Dx") & obs.malignant.eq("malignant")
    counts = a.layers["counts"]
    # obs.total_counts was computed before gene filtering. The saved count
    # layer excludes a small fraction of genes, so require close agreement
    # in the selected cells rather than exact identity.
    layer_totals = np.asarray(counts.sum(axis=1)).ravel()
    recorded_totals = obs.total_counts.to_numpy(dtype=float)
    rel = np.abs(layer_totals - recorded_totals) / np.maximum(recorded_totals, 1)
    if np.quantile(rel[mask.to_numpy()], .99) >= .01 or rel[mask.to_numpy()].max() >= .05:
        raise ValueError("GSE227122 selected-cell count totals diverge from metadata")
    index = pd.Index(a.var_names)
    rows = []
    for patient in sorted(obs.loc[mask, "patient"].unique()):
        selected = (mask & obs.patient.eq(patient)).to_numpy()
        pooled = np.asarray(counts[selected, :].sum(axis=0)).ravel()
        rows.append(row("GSE227122", str(patient), pooled, index, int(selected.sum()),
                        "project heuristic candidate-malignant"))
    if len(rows) != 10:
        raise ValueError(f"GSE227122 expected 10 patients; got {len(rows)}")
    return rows


def gse248287(path: Path) -> list[dict]:
    rows = []
    for folder in sorted(path.glob("P*_Dx"), key=lambda p: int(p.name.split("_")[0][1:])):
        meta, features, counts = load_sample(folder)
        selected = meta.Celltypes_all.eq("Malignant").to_numpy()
        pooled = np.asarray(counts[:, selected].sum(axis=1)).ravel()
        rows.append(row("GSE248287", folder.name.split("_")[0], pooled,
                        pd.Index(features), int(selected.sum()), "author Celltypes_all=Malignant"))
    if len(rows) != 15:
        raise ValueError(f"GSE248287 expected 15 patients; got {len(rows)}")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse227122", type=Path, required=True)
    parser.add_argument("--gse248287", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(gse227122(args.gse227122) + gse248287(args.gse248287))
    if result[["cohort", "patient"]].duplicated().any():
        raise ValueError("duplicate patient")
    result.to_csv(args.output, sep="\t", index=False)
    print(f"saved {args.output}: {result.groupby('cohort').size().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
