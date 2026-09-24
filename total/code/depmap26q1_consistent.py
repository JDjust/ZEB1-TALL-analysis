#!/usr/bin/env python3
"""Rebuild Figure 7 C/D/E from one DepMap Public 26Q1 release.

The model table uses 'B-Cell Acute Lymphoblastic Leukemia' for B-ALL. Earlier
scripts looked for a different name and silently mixed releases in the figure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

T_LABEL = "T-Lymphoblastic Leukemia/Lymphoma"
B_LABEL = "B-Cell Acute Lymphoblastic Leukemia"


def selected_table(path: Path, genes: list[str], *, first_col_index: bool = True):
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    chosen = {g: next((c for c in cols if c.startswith(g + " (")), None) for g in genes}
    if any(v is None for v in chosen.values()):
        raise ValueError(f"missing genes in {path.name}: {[g for g, v in chosen.items() if v is None]}")
    keep = [cols[0]] + list(chosen.values())
    x = pd.read_csv(path, usecols=keep)
    x = x.rename(columns={cols[0]: "ModelID", **{v: g for g, v in chosen.items()}})
    if x["ModelID"].duplicated().any():
        raise ValueError(f"duplicate ModelID in {path.name}")
    return x


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--genes", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    genes = [x.strip() for x in args.genes.read_text().splitlines() if x.strip()]
    if len(genes) != len(set(genes)) or "ZEB1" not in genes:
        raise ValueError("invalid target gene list")
    paths = {name: args.data / name for name in [
        "Model.csv", "CRISPRGeneEffect.csv", "CRISPRGeneDependency.csv",
        "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"]}
    model = pd.read_csv(paths["Model.csv"])[["ModelID", "CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"]]
    if model.ModelID.duplicated().any():
        raise ValueError("duplicate ModelID in Model.csv")
    model["lineage"] = np.select([
        model.OncotreePrimaryDisease.eq(T_LABEL),
        model.OncotreePrimaryDisease.eq(B_LABEL),
        model.OncotreeLineage.eq("Myeloid")],
        ["T-ALL", "B-ALL", "myeloid"], default="other")
    ge = selected_table(paths["CRISPRGeneEffect.csv"], genes)
    ge = ge.rename(columns={g: g + "_GE" for g in genes})
    dep = selected_table(paths["CRISPRGeneDependency.csv"], ["ZEB1"])
    dep = dep.rename(columns={"ZEB1": "ZEB1_dependency"})
    expr_cols = pd.read_csv(paths["OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"], nrows=0).columns
    zeb_col = next((c for c in expr_cols if c.startswith("ZEB1 (")), None)
    if zeb_col is None:
        raise ValueError("ZEB1 expression absent")
    expr = pd.read_csv(paths["OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"],
                       usecols=["ModelID", "IsDefaultEntryForModel", zeb_col])
    default = expr["IsDefaultEntryForModel"].astype(str).str.lower().isin(["true", "1", "yes"])
    expr = expr.loc[default, ["ModelID", zeb_col]].rename(columns={zeb_col: "ZEB1_RNA"})
    if expr.ModelID.duplicated().any():
        raise ValueError("multiple default expression entries for one ModelID")
    merged = model.merge(ge, on="ModelID", how="left", validate="one_to_one")
    merged = merged.merge(dep, on="ModelID", how="left", validate="one_to_one")
    merged = merged.merge(expr, on="ModelID", how="left", validate="one_to_one")
    tall = merged[merged.lineage.eq("T-ALL")].copy()
    if len(tall) != 28 or tall.ZEB1_GE.notna().sum() != 8 or tall.ZEB1_RNA.notna().sum() != 18:
        raise ValueError("26Q1 T-ALL cohort changed; re-check release")
    selected = merged[merged.lineage.isin(["T-ALL", "B-ALL", "myeloid"])].copy()
    long = selected.melt(id_vars=["ModelID", "CellLineName", "lineage", "ZEB1_RNA",
                                  "ZEB1_dependency"],
                         value_vars=[g + "_GE" for g in genes], var_name="gene", value_name="gene_effect")
    long["gene"] = long["gene"].str.removesuffix("_GE")
    summary = (long.groupby(["gene", "lineage"], observed=True)["gene_effect"]
               .agg(n="count", mean="mean", median="median").reset_index())
    pivot = summary.pivot(index="gene", columns="lineage", values=["n", "mean"])
    pivot.columns = [f"{stat}_{lineage.replace('-', '')}" for stat, lineage in pivot.columns]
    pivot = pivot.reset_index()
    for field in ["n_TALL", "n_BALL", "n_myeloid", "mean_TALL", "mean_BALL", "mean_myeloid"]:
        if field not in pivot.columns:
            pivot[field] = np.nan
    pivot["delta_T_minus_B"] = pivot["mean_TALL"] - pivot["mean_BALL"]
    pivot["delta_T_minus_M"] = pivot["mean_TALL"] - pivot["mean_myeloid"]
    pivot.to_csv(args.output / "target_gene_effect_summary.tsv", sep="\t", index=False)
    tall[["ModelID", "CellLineName", "ZEB1_RNA", "ZEB1_GE", "ZEB1_dependency"]].to_csv(
        args.output / "tall_zeb1_dependency.tsv", sep="\t", index=False)
    long.to_csv(args.output / "selected_model_gene_effects.tsv", sep="\t", index=False)
    provenance = {"release": "DepMap Public 26Q1", "n_models_by_lineage": model.lineage.value_counts().to_dict(),
                  "n_tall_crispr": int(tall.ZEB1_GE.notna().sum()), "n_genes": len(genes),
                  "gene_list_sha256": sha256(args.genes),
                  "input_sha256": {name: sha256(path) for name, path in paths.items()}}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps({"n_tall_crispr": provenance["n_tall_crispr"],
                      "n_b_all_crispr_zeb1": int(merged.loc[merged.lineage.eq("B-ALL"), "ZEB1_GE"].notna().sum())}))


if __name__ == "__main__":
    main()
