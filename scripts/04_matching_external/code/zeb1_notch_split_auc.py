#!/usr/bin/env python3
"""Split the 10-gene NOTCH score and rescore the same 8 DepMap models."""
from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1")
OUT = DATA / "_zeb1_forensic"

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
SETS = {
    "canonical_NOTCH_TF": ["HES1", "HES4", "HEY1", "DTX1", "NRARP"],
    "developmental_NOTCH_output": ["PTCRA", "NOTCH3", "IL7R"],
    "MYC_only": ["MYC"],
    "original_10": ["HES1", "HES4", "HEY1", "DTX1", "NRARP", "MYC", "PTCRA", "NOTCH3", "IL7R", "SHQ1"],
    "original_minus_MYC_PTCRA_IL7R": ["HES1", "HES4", "HEY1", "DTX1", "NRARP", "NOTCH3", "SHQ1"],
}


def gene_symbol(col):
    return str(col).split(" (")[0].strip()


def auc(values, dep_ids, non_ids):
    d = values.reindex(dep_ids).to_numpy(float)
    n = values.reindex(non_ids).to_numpy(float)
    wins = ties = 0
    for a in d:
        for b in n:
            if a > b:
                wins += 1
            elif a == b:
                ties += 1
    return (wins + 0.5 * ties) / (len(d) * len(n))


def main():
    df = pd.read_csv(DATA / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", low_memory=False)
    df = df[df["ModelID"].astype(str).isin(GROUP)].copy()
    df = df.rename(columns={c: gene_symbol(c) for c in df.columns if c != "ModelID"})
    df["ModelID"] = df["ModelID"].astype(str)
    mat = df.groupby("ModelID").mean(numeric_only=True)
    order = list(GROUP)
    mat = mat.reindex(order)
    dep = [i for i, (_, g) in GROUP.items() if g == "dependent"]
    non = [i for i, (_, g) in GROUP.items() if g == "nondependent"]
    rows = []
    wide = {}
    for name, genes in SETS.items():
        present = [g for g in genes if g in mat.columns]
        sub = mat[present]
        if name == "MYC_only":
            score = sub.iloc[:, 0]
        else:
            z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
            score = z.mean(axis=1)
        wide[name] = score
        rows.append({
            "set": name,
            "n_genes": len(present),
            "missing": ",".join(g for g in genes if g not in mat.columns),
            "auc_dep_higher": auc(score, dep, non),
            **{GROUP[i][0]: float(score.loc[i]) for i in order},
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "notch_split_auc.tsv", sep="\t", index=False)
    print(out[["set", "n_genes", "auc_dep_higher"]].to_string(index=False))
    print(out.drop(columns=["missing"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
