#!/usr/bin/env python3
"""T-ALL versus B-ALL gene-effect uncertainty in locked DepMap 26Q1 models."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    data = pd.read_csv(args.input, sep="\t")
    rng = np.random.default_rng(20260923)
    rows = []
    for gene, frame in data.groupby("gene"):
        t = frame.loc[frame.lineage.eq("T-ALL"), "gene_effect"].dropna().to_numpy()
        b = frame.loc[frame.lineage.eq("B-ALL"), "gene_effect"].dropna().to_numpy()
        if len(t) < 4 or len(b) < 4:
            continue
        boot = np.array([rng.choice(t, len(t), replace=True).mean() -
                         rng.choice(b, len(b), replace=True).mean() for _ in range(5000)])
        rows.append({"gene": gene, "n_TALL": len(t), "n_BALL": len(b),
                     "mean_delta_T_minus_B": float(t.mean() - b.mean()),
                     "bootstrap_ci95_lo": float(np.quantile(boot, 0.025)),
                     "bootstrap_ci95_hi": float(np.quantile(boot, 0.975)),
                     "mann_whitney_p": float(mannwhitneyu(t, b, alternative="two-sided").pvalue)})
    out = pd.DataFrame(rows).sort_values("mann_whitney_p").reset_index(drop=True)
    n = len(out)
    q = np.minimum.accumulate((out.mann_whitney_p.to_numpy() * n /
                               np.arange(1, n + 1))[::-1])[::-1]
    out["fdr_bh_32_genes"] = np.clip(q, 0, 1)
    out.to_csv(args.output, sep="\t", index=False)
    print(out[out.gene.isin(["ZEB1", "GATA3", "BCL11B"])].to_string(index=False))


if __name__ == "__main__":
    main()
