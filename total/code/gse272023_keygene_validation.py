#!/usr/bin/env python3
"""Independent NOPHO GSE272023 patient-level ZEB1/key-gene check."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, pearsonr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "validation"
GENES = ["ZEB1", "GATA3", "TCF7", "BCL11B", "IL7R", "LMO2", "ZEB2"]


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) /
                               np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p)); result[order] = np.minimum(q, 1)
    return result


def main():
    mat = pd.read_csv(DATA / "GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz", sep="\t")
    meta = pd.read_csv(DATA / "gse272023_sample_map.tsv", sep="\t")
    mapping = pd.read_csv(DATA / "gse272023_keygene_id_map.tsv", sep="\t")
    if mat.columns[0] != "ID" or mat.ID.duplicated().any():
        raise ValueError("GSE272023 matrix gene schema changed")
    if set(mapping.gene_symbol) != set(GENES) or mapping.ensembl_gene.duplicated().any():
        raise ValueError("key-gene mapping invalid")
    tumors = meta[meta.tissue_type.eq("Tumor")].copy()
    if len(tumors) != 108 or tumors["sample"].duplicated().any():
        raise ValueError("GSE272023 tumor sample metadata changed")
    if not set(tumors["sample"]).issubset(mat.columns):
        raise ValueError("tumor sample IDs not in VST matrix")
    sub = mat.set_index("ID").loc[mapping.ensembl_gene, tumors["sample"]].T
    sub.columns = mapping.gene_symbol.tolist()
    if sub.isna().any().any() or not np.isfinite(sub.to_numpy(dtype=float)).all():
        raise ValueError("nonfinite selected expression")
    d = tumors.set_index("sample")[["gsm", "cimp", "tissue_type"]].join(sub)
    if d.cimp.value_counts().to_dict() != {"CIMP high": 65, "CIMP low": 43}:
        raise ValueError("CIMP group sizes changed")
    # Rank residualization estimates a conditional association; CIMP class is
    # coarse and is not a complete molecular subtype adjustment.
    cov = np.column_stack([np.ones(len(d)), d.cimp.eq("CIMP high").to_numpy(dtype=float)])
    def residual_rank(values):
        y = rankdata(values)
        return y - cov @ np.linalg.lstsq(cov, y, rcond=None)[0]
    zres = residual_rank(d.ZEB1)
    records = []
    for gene in GENES[1:]:
        rho, p = spearmanr(d.ZEB1, d[gene])
        r_cond, p_cond = pearsonr(zres, residual_rank(d[gene]))
        records.append(dict(gene=gene, n_patients=len(d), rho=float(rho),
                            p_exploratory=float(p), cimp_partial_rank_r=float(r_cond),
                            cimp_partial_p_exploratory=float(p_cond),
                            expression="deposited VST counts", unit="patient"))
    out = pd.DataFrame(records)
    out["fdr_bh_6_genes"] = bh(out.p_exploratory)
    out["cimp_partial_fdr_bh_6_genes"] = bh(out.cimp_partial_p_exploratory)
    d.reset_index().to_csv(DATA / "gse272023_keygene_patients.tsv", sep="\t", index=False)
    out.to_csv(DATA / "gse272023_keygene_associations.tsv", sep="\t", index=False)
    (DATA / "gse272023_keygene_run.json").write_text(json.dumps({
        "accession": "GSE272023", "n_tumor_patients": len(d),
        "n_cimp_high": 65, "n_cimp_low": 43,
        "input_matrix": "GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz",
        "unit": "deposited VST expression per patient",
        "gene_mapping": "gse272023_keygene_id_map.tsv",
        "caution": "CIMP rank adjustment is not molecular subtype or batch adjustment",
    }, indent=2), encoding="utf-8")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
