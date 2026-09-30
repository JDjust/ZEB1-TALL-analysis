"""Descriptive condition pseudobulk using author genotype/time/priming labels."""
import argparse
import gzip
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument("--raw", type=Path, required=True)
parser.add_argument("--metadata", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
meta = pd.read_csv(args.metadata, sep="\t")
barcodes = pd.read_csv(args.raw / "GSM9286835_cDNA_barcodes.tsv.gz", header=None)[0]
features = pd.read_csv(args.raw / "GSM9286835_cDNA_features.tsv.gz", sep="\t", header=None)
assert barcodes.is_unique and meta["index"].is_unique
positions = pd.Index(barcodes).get_indexer(meta["index"])
assert (positions >= 0).all(), "Author barcode missing from raw matrix"
with gzip.open(args.raw / "GSM9286835_cDNA_matrix.mtx.gz", "rb") as stream:
    matrix = mmread(stream).tocsc()
assert matrix.shape == (len(features), len(barcodes))
genes = ["Zeb1", "Lmo2", "Tcf7", "Gata3", "Bcl11b", "Il7r", "Lyl1", "Tal1", "Cd34", "Kit", "Spi1"]
gene_indices = {}
for gene in genes:
    hits = features.index[features[1].eq(gene)]
    assert len(hits) == 1, f"Ambiguous gene: {gene}"
    gene_indices[gene] = hits[0]
records = []
for condition in sorted(meta.hto.unique()):
    selected = positions[meta.hto.eq(condition)]
    counts = np.asarray(matrix[:, selected].sum(axis=1)).ravel()
    total = float(counts.sum())
    timepoint, priming, genotype = condition.split("_")
    for gene, index in gene_indices.items():
        records.append(dict(condition=condition, timepoint=timepoint, priming=priming,
                            genotype=genotype, gene=gene, n_cells=len(selected),
                            library_umi=int(total), gene_umi=int(counts[index]),
                            log2_cpm=float(np.log2(1 + counts[index] * 1e6 / total))))
result = pd.DataFrame(records)
result.to_csv(args.output / "author_condition_pseudobulk.tsv", sep="\t", index=False)
pivot = result.pivot(index=["timepoint", "priming", "gene"], columns="genotype", values="log2_cpm")
assert set(pivot.columns) == {"EV", "KO"} and pivot.notna().all().all()
pivot["delta_KO_minus_EV"] = pivot.KO - pivot.EV
pivot.reset_index().to_csv(args.output / "author_condition_contrasts.tsv", sep="\t", index=False)
(args.output / "analysis_scope.json").write_text(json.dumps({"unit": "author-labelled experimental condition",
    "n_conditions": 8, "n_author_retained_cells": len(meta), "n_matched_cells": len(positions),
    "genes": genes, "statistics": "Descriptive log2(1+CPM) differences; no cell-level P or FDR",
    "limitation": "Four context pairs are not independent biological replicates of one perturbation"}, indent=2))
print(pivot.loc[(slice(None), slice(None), ["Zeb1", "Lmo2"]), :].to_string(), flush=True)
