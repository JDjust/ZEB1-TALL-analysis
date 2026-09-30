#!/usr/bin/env python3
"""Separate 4 ZEB1-dependent vs 4 non-dependent T-ALL models in DepMap 26Q1.

Pre-specified axes are scored as one number per model. A genome-wide scan is
saved only as a calibration: with 4 vs 4, perfect separation is expected for
hundreds of genes by chance.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1")
OUT = DATA / "_zeb1_forensic"
OUT.mkdir(exist_ok=True)

# Dependency probability from the locked Fig S9 table, not re-derived here.
GROUP = {
    "ACH-000981": ("DND-41", "dependent"),
    "ACH-000197": ("TALL-1", "dependent"),
    "ACH-000995": ("JURKAT", "dependent"),
    "ACH-000953": ("SUP-T1", "dependent"),
    "ACH-000937": ("PF-382", "nondependent"),
    "ACH-000519": ("PEER", "nondependent"),
    "ACH-000101": ("KE-37", "nondependent"),
    "ACH-001737": ("CCRF-HSB-2-DM", "nondependent"),
}

# One score per axis. Membership is a literature panel, not a fitted signature.
AXES = {
    "NOTCH_targets": [
        "HES1", "HES4", "HEY1", "DTX1", "NRARP", "MYC", "PTCRA", "NOTCH3", "IL7R", "SHQ1",
    ],
    "NOTCH_receptors_ligands": ["NOTCH1", "NOTCH2", "NOTCH3", "JAG1", "DLL1", "DLL4", "FBXW7"],
    "TAL1_LMO": ["TAL1", "LMO1", "LMO2", "LYL1", "GATA3", "MYB", "RUNX1", "ERG"],
    "TLX_NKX": ["TLX1", "TLX3", "NKX2-1", "NKX2-5", "RANBP17"],
    "ETP_immature": ["LYL1", "LMO2", "HHEX", "MEF2C", "SPI1", "CD34", "KIT", "BCL2", "CD44"],
    "cortical_mature_T": [
        "BCL11B", "TCF7", "LEF1", "CD1A", "CD4", "CD8A", "CD3D", "CD3E", "RAG1", "PTCRA", "LCK", "ZAP70",
    ],
    "MYC_targets_core": ["MYC", "MYCN", "ODC1", "LDHA", "PKM", "NCL", "CAD", "NOP56", "WDR43", "NPM1"],
    "PTEN_PI3K_AKT_RNA": ["PTEN", "PIK3CA", "PIK3CD", "AKT1", "MTOR", "RICTOR"],
    "TCR_signaling": ["CD3D", "CD3E", "CD3G", "LCK", "ZAP70", "LAT", "ITK", "FYN", "TRAC", "CD28"],
    "chromatin_regulators": [
        "EZH2", "SUZ12", "EED", "KMT2A", "KDM6A", "HDAC1", "CREBBP", "EP300", "PHF6", "SMARCA4",
    ],
    "ZEB1_EMT_regulon": ["ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH1", "EPCAM", "CTBP1", "CTBP2"],
    "folate_MTX": [
        "DHFR", "TYMS", "FPGS", "GGH", "SLC19A1", "MTHFR", "ATIC", "GART", "SHMT1", "SHMT2", "MTHFD1", "MTHFD2",
    ],
}

MARKERS = sorted({
    "ZEB1", "ZEB2", "NOTCH1", "NOTCH3", "FBXW7", "HES1", "HES4", "DTX1", "NRARP", "MYC", "PTCRA",
    "TAL1", "LMO1", "LMO2", "LYL1", "TLX1", "TLX3", "NKX2-1", "BCL11B", "TCF7", "LEF1", "RAG1",
    "CD1A", "CD34", "HHEX", "MEF2C", "PTEN", "PIK3CA", "AKT1", "MYB", "GATA3", "IL7R", "LCK",
    "CD3D", "TRAC", "DHFR", "TYMS", "FPGS", "SLC19A1", "EZH2", "PHF6", "CREBBP", "EP300",
    "CDKN2A", "CDKN2B", "RB1", "TP53", "WT1", "JAK1", "STAT5B", "USP7",
})

CN_GENES = [
    "NOTCH1", "FBXW7", "PTEN", "CDKN2A", "CDKN2B", "MYC", "MYB", "TAL1", "LMO2", "TLX1", "TLX3",
    "PHF6", "LEF1", "BCL11B", "RB1", "TP53", "EZH2", "CREBBP", "EP300", "ZEB1", "WT1", "IL7R",
]


def gene_symbol(col: str) -> str:
    return str(col).split(" (")[0].strip()


def load_models(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    if "ModelID" in df.columns:
        raw_id = "ModelID"
    else:
        raw_id = df.columns[0]
    keep = df[raw_id].astype(str).isin(GROUP)
    out = df.loc[keep].copy()
    rename = {c: gene_symbol(c) for c in out.columns if c != raw_id}
    out = out.rename(columns=rename)
    out[raw_id] = out[raw_id].astype(str)
    meta_cols = ["SequencingID", "ModelConditionID", "IsDefaultEntryForModel", "IsDefaultEntryForMC"]
    out = out.drop(columns=[c for c in meta_cols if c in out.columns])
    out = out.groupby(raw_id, as_index=True).mean(numeric_only=True)
    out = out.loc[:, ~out.columns.duplicated()]
    out.index.name = "ModelID"
    return out


def auc_high_in_dep(values: pd.Series, dep_ids, non_ids) -> float:
    d = values.reindex(dep_ids).to_numpy(dtype=float)
    n = values.reindex(non_ids).to_numpy(dtype=float)
    if np.isnan(d).any() or np.isnan(n).any():
        return np.nan
    wins = ties = 0
    for a in d:
        for b in n:
            if a > b:
                wins += 1
            elif a == b:
                ties += 1
    return (wins + 0.5 * ties) / (len(d) * len(n))


def axis_score(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in mat.columns]
    if len(present) < 3:
        return pd.Series(np.nan, index=mat.index)
    sub = mat[present]
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def scan(mat: pd.DataFrame, dep_ids, non_ids, kind: str) -> pd.DataFrame:
    rows = []
    for gene in mat.columns:
        vals = mat[gene]
        if vals.notna().sum() < 8:
            continue
        a = auc_high_in_dep(vals, dep_ids, non_ids)
        if np.isnan(a):
            continue
        if a >= 0.875 or a <= 0.125:
            rows.append((gene, a, float(vals.reindex(dep_ids).mean() - vals.reindex(non_ids).mean())))
    out = pd.DataFrame(rows, columns=["gene", "auc_dep_higher", "mean_diff_dep_minus_non"])
    out.insert(0, "matrix", kind)
    return out.sort_values("auc_dep_higher", ascending=False)


def main() -> None:
    dep_ids = [i for i, (_, g) in GROUP.items() if g == "dependent"]
    non_ids = [i for i, (_, g) in GROUP.items() if g == "nondependent"]

    expr = load_models(DATA / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    crispr = load_models(DATA / "CRISPRGeneEffect.csv")
    cn = load_models(DATA / "OmicsCNGeneWGS.csv")
    order = list(GROUP)
    expr = expr.reindex(order)
    crispr = crispr.reindex(order)
    cn = cn.reindex(order)
    print("expr models", int(expr.notna().any(axis=1).sum()), "genes", expr.shape[1])
    print("crispr models", int(crispr.notna().any(axis=1).sum()), "genes", crispr.shape[1])
    print("cn models", int(cn.notna().any(axis=1).sum()), "genes", cn.shape[1])

    meta = pd.DataFrame(
        [{"ModelID": i, "CellLineName": n, "group": g} for i, (n, g) in GROUP.items()]
    ).set_index("ModelID")

    axis_rows = []
    for name, genes in AXES.items():
        for kind, mat in (("expression", expr), ("CRISPR_gene_effect", crispr)):
            score = axis_score(mat, genes)
            present = [g for g in genes if g in mat.columns]
            a = auc_high_in_dep(score, dep_ids, non_ids)
            axis_rows.append({
                "axis": name,
                "matrix": kind,
                "n_genes_present": len(present),
                "genes_missing": ",".join(g for g in genes if g not in mat.columns),
                "auc_dep_higher": a,
                **{GROUP[i][0]: float(score.loc[i]) for i in GROUP},
            })
    axes = pd.DataFrame(axis_rows)
    axes.to_csv(OUT / "axis_scores.tsv", sep="\t", index=False)

    marker_e = expr.reindex(columns=[g for g in MARKERS if g in expr.columns])
    marker_c = crispr.reindex(columns=[g for g in MARKERS if g in crispr.columns])
    marker_e.insert(0, "CellLineName", [GROUP[i][0] for i in marker_e.index])
    marker_e.insert(1, "group", [GROUP[i][1] for i in marker_e.index])
    marker_c.insert(0, "CellLineName", [GROUP[i][0] for i in marker_c.index])
    marker_c.insert(1, "group", [GROUP[i][1] for i in marker_c.index])
    marker_e.to_csv(OUT / "marker_expression.tsv", sep="\t")
    marker_c.to_csv(OUT / "marker_crispr.tsv", sep="\t")

    cn_sub = cn.reindex(columns=[g for g in CN_GENES if g in cn.columns])
    cn_sub.insert(0, "CellLineName", [GROUP[i][0] for i in cn_sub.index])
    cn_sub.insert(1, "group", [GROUP[i][1] for i in cn_sub.index])
    cn_sub.to_csv(OUT / "marker_copynumber.tsv", sep="\t")

    expr_scan = scan(expr, dep_ids, non_ids, "expression")
    crispr_scan = scan(crispr, dep_ids, non_ids, "CRISPR_gene_effect")
    both = pd.concat([expr_scan, crispr_scan], ignore_index=True)
    both.to_csv(OUT / "near_perfect_separators.tsv", sep="\t", index=False)

    n_genes_e = int(expr.notna().all().sum())
    n_genes_c = int(crispr.notna().all().sum())
    # One-sided perfect split among 8 distinct values: 1/C(8,4) = 1/70.
    summary = pd.DataFrame([
        {"matrix": "expression", "n_complete_genes": n_genes_e,
         "n_auc_1": int((expr_scan.auc_dep_higher == 1).sum()),
         "n_auc_0": int((expr_scan.auc_dep_higher == 0).sum()),
         "expected_one_direction": n_genes_e / 70},
        {"matrix": "CRISPR_gene_effect", "n_complete_genes": n_genes_c,
         "n_auc_1": int((crispr_scan.auc_dep_higher == 1).sum()),
         "n_auc_0": int((crispr_scan.auc_dep_higher == 0).sum()),
         "expected_one_direction": n_genes_c / 70},
    ])
    summary.to_csv(OUT / "separator_null.tsv", sep="\t", index=False)
    print(axes[["axis", "matrix", "n_genes_present", "auc_dep_higher"]].to_string(index=False))
    print(summary.to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
