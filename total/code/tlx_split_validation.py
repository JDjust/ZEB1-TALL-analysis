#!/usr/bin/env python3
"""Exploratory TLX1/TLX3 split in two independent public T-ALL cohorts.

Input is the unchanged Hulu module3/processed/independent_keygenes.tsv export.
These subtype calls are study supplied; the analysis is post hoc and descriptive.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/validation/independent_keygenes.tsv"
OUT = ROOT / "data/validation/tlx1_tlx3_split.tsv"


def effect(a, b):
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na - 1) * np.var(a, ddof=1) +
                      (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    g = (np.mean(a) - np.mean(b)) / pooled * (1 - 3 / (4 * (na + nb) - 9))
    return float(g)


def main():
    d = pd.read_csv(INPUT, sep="\t")
    d = d[d.disease.eq("T-ALL")].copy()
    if d.gsm.duplicated().any():
        raise ValueError("sample IDs are duplicated")
    rows = []
    for cohort, labels in {
        "GSE26713": {"TLX1": "TLX1", "TLX3": "TLX3"},
        "GSE62156": {"TLX1": "tlx1 subgroup", "TLX3": "tlx3 subgroup"},
    }.items():
        c = d[d.cohort.eq(cohort)].copy()
        tlx = c[c.level1.eq("TLX")]
        rest = c[~c.level1.eq("TLX")]
        assert len(tlx) == sum(tlx.level2.eq(v).sum() for v in labels.values())
        for subtype, label in labels.items():
            a = tlx.loc[tlx.level2.eq(label), "ZEB1"].dropna().to_numpy()
            b = rest.ZEB1.dropna().to_numpy()
            if len(a) < 3 or len(b) < 3:
                raise ValueError("insufficient patient units")
            rows.append({"cohort": cohort, "subtype": subtype,
                         "comparator": "all non-TLX T-ALL in cohort",
                         "n_subtype": len(a), "n_comparator": len(b),
                         "mean_subtype": np.mean(a), "mean_comparator": np.mean(b),
                         "hedges_g": effect(a, b),
                         "wilcoxon_p_exploratory": mannwhitneyu(a, b, alternative="two-sided").pvalue,
                         "patient_unit": "GSM", "inference": "post-hoc, two cohorts only"})
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    print(OUT)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
