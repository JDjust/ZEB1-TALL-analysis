#!/usr/bin/env python3
"""Pan-cancer ZEB1 CRISPR co-dependency, plus a pre-specified NOTCH-machinery panel."""
from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/CRISPRGeneEffect.csv")
OUT = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/_zeb1_forensic")

PANEL = [
    "RBPJ", "MAML1", "MAML2", "MAML3", "ADAM10", "PSEN1", "NCSTN", "APH1A",
    "HES1", "HES4", "HEY1", "DTX1", "USP7", "BRD4", "EP300", "CREBBP",
    "NOTCH1", "FBXW7", "MYC", "BCL11B", "TCF7",
]
# Class labels are descriptive only. They are not used to rank co-dependencies.
CLASS = {
    "ACH-000981": "A_NOTCH_high_dep",
    "ACH-000197": "A_NOTCH_high_dep",
    "ACH-000953": "A_NOTCH_high_dep",
    "ACH-000995": "B_Jurkat",
    "ACH-000937": "C_nondependent",
    "ACH-000519": "C_nondependent",
    "ACH-000101": "C_nondependent",
    "ACH-001737": "C_nondependent",
}


def symbol(col):
    return str(col).split(" (")[0].strip()


def main():
    df = pd.read_csv(DATA, low_memory=False)
    df = df.rename(columns={c: symbol(c) for c in df.columns})
    df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.loc[:, ~df.columns.duplicated()]
    if "ZEB1" not in df.columns:
        raise SystemExit("ZEB1 missing")
    z = df["ZEB1"]
    use = df.columns[df.columns != "ZEB1"]
    # pairwise pearson with ZEB1; positive rho = co-dependent (both more negative together)
    rows = []
    zv = z.to_numpy(float)
    for gene in use:
        y = df[gene].to_numpy(float)
        m = np.isfinite(zv) & np.isfinite(y)
        n = int(m.sum())
        if n < 100:
            continue
        a = zv[m]
        b = y[m]
        a = a - a.mean()
        b = b - b.mean()
        denom = np.sqrt((a * a).sum() * (b * b).sum())
        rho = float((a * b).sum() / denom) if denom else np.nan
        rows.append((gene, n, rho))
    res = pd.DataFrame(rows, columns=["gene", "n", "pearson_with_ZEB1"])
    res = res.sort_values("pearson_with_ZEB1", ascending=False)
    res.head(40).to_csv(OUT / "zeb1_codependency_top.tsv", sep="\t", index=False)
    res.tail(15).to_csv(OUT / "zeb1_codependency_bottom.tsv", sep="\t", index=False)
    panel = res[res.gene.isin(PANEL)].sort_values("pearson_with_ZEB1", ascending=False)
    panel.to_csv(OUT / "zeb1_codependency_panel.tsv", sep="\t", index=False)
    print("n_lines", int(z.notna().sum()), "n_genes", len(res))
    print(res.head(20).to_string(index=False))
    print("--- panel ---")
    print(panel.to_string(index=False))
    # values of the panel in the eight lines, not used to choose the panel
    sub = df.loc[df.index.isin(CLASS), ["ZEB1"] + [g for g in PANEL if g in df.columns]]
    sub.insert(0, "class", [CLASS.get(i, "") for i in sub.index])
    sub.to_csv(OUT / "zeb1_codependency_eight.tsv", sep="\t")
    print(sub.round(3).to_string())


if __name__ == "__main__":
    main()
