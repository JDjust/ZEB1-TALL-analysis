#!/usr/bin/env python3
"""Build Figure 4B from patient-level ZEB1 correlations in three bulk cohorts."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "data" / "validation"
MODULE = ROOT.parent / "modules" / "module4" / "tables"
OUT = ROOT / "data" / "SourceData_Fig4B_crosscohort.tsv"
GENES = ["GATA3", "TCF7", "BCL11B", "IL7R", "LMO2", "ZEB2"]


def main():
    target = pd.read_csv(MODULE / "M4_keygene_stats.tsv", sep="\t")
    pharma = pd.read_csv(VAL / "pharmacotype_zeb1_gene_stats.tsv", sep="\t")
    gse = pd.read_csv(VAL / "gse272023_keygene_associations.tsv", sep="\t")
    rows = []
    for cohort, frame, n, family in [
        ("TARGET-ALL-P2", target, 265, "module4 target gene-stat family"),
        ("StJude_Pharmacotype", pharma, 116, "module4 pharmacotype gene-stat family"),
    ]:
        d = frame[frame.gene.isin(GENES)]
        if len(d) != len(GENES) or not d.n.eq(n).all():
            raise ValueError(f"{cohort}: key gene/sample mismatch")
        for _, r in d.iterrows():
            rows.append(dict(gene=r.gene, cohort=cohort, n_patients=n,
                             rho=float(r.rho_zeb1), p=float(r.p_rho),
                             fdr=float(r.fdr_rho), fdr_family=family,
                             method="unadjusted patient Spearman"))
    if set(gse.gene) != set(GENES) or not gse.n_patients.eq(108).all():
        raise ValueError("GSE272023 key gene/sample mismatch")
    for _, r in gse.iterrows():
        rows.append(dict(gene=r.gene, cohort="GSE272023", n_patients=108,
                         rho=float(r.rho), p=float(r.p_exploratory),
                         fdr=float(r.fdr_bh_6_genes),
                         fdr_family="six selected genes", method="unadjusted patient Spearman"))
    out = pd.DataFrame(rows)
    if len(out) != 18 or out[["rho", "p", "fdr"]].isna().any().any():
        raise ValueError("Figure 4B source table incomplete")
    out.to_csv(OUT, sep="\t", index=False)
    print(f"saved {OUT}: {len(out)} rows")


if __name__ == "__main__":
    main()
