#!/usr/bin/env python3
"""Check a frozen GATA3/TCF7/BCL11B signal in two patient scRNA cohorts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "validation" / "three_gene_patient"
INPUT = DATA / "patient_raw_umi_four_genes.tsv"
GENES = ["GATA3", "TCF7", "BCL11B"]
SEED = 20260924


def perm_p(x: np.ndarray, y: np.ndarray, n_perm: int, rng) -> float:
    xr, yr = rankdata(x), rankdata(y)
    xr = (xr - xr.mean()) / np.linalg.norm(xr - xr.mean())
    yr = (yr - yr.mean()) / np.linalg.norm(yr - yr.mean())
    observed = abs(float(xr @ yr))
    count = 0
    for _ in range(n_perm):
        count += abs(float(xr @ rng.permutation(yr))) >= observed - 1e-12
    return (count + 1) / (n_perm + 1)


def boot_ci(x: np.ndarray, y: np.ndarray, n_boot: int, rng) -> tuple[float, float]:
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(x[idx])) < 3 or len(np.unique(y[idx])) < 3:
            continue
        vals.append(float(spearmanr(x[idx], y[idx]).statistic))
    return tuple(np.quantile(vals, [0.025, 0.975])) if vals else (np.nan, np.nan)


def main():
    d = pd.read_csv(INPUT, sep="\t")
    if set(d.cohort) != {"GSE227122", "GSE248287"} or d[["cohort", "patient"]].duplicated().any():
        raise ValueError("patient cohort schema changed")
    if d.groupby("cohort").size().to_dict() != {"GSE227122": 10, "GSE248287": 15}:
        raise ValueError("patient counts changed")
    rng = np.random.default_rng(SEED)
    out, results = [], []
    for cohort, g in d.groupby("cohort", sort=True):
        g = g.copy().sort_values("patient")
        z = g[GENES].to_numpy(dtype=float)
        if not np.isfinite(z).all() or not np.isfinite(g.ZEB1).all():
            raise ValueError(f"{cohort}: nonfinite expression")
        std = z.std(axis=0, ddof=1)
        if (std <= 0).any():
            raise ValueError(f"{cohort}: constant signature gene")
        g["three_gene_score"] = ((z - z.mean(axis=0)) / std).mean(axis=1)
        g["score_definition"] = "mean of within-cohort patient z scores: GATA3,TCF7,BCL11B"
        out.append(g)
        for feature in ["three_gene_score"] + GENES:
            x = g.ZEB1.to_numpy(dtype=float)
            y = g[feature].to_numpy(dtype=float)
            rho = float(spearmanr(x, y).statistic)
            loo = [float(spearmanr(np.delete(x, i), np.delete(y, i)).statistic)
                   for i in range(len(x))]
            lo, hi = boot_ci(x, y, 10000, rng)
            results.append(dict(cohort=cohort, feature=feature, n_patients=len(g), rho=rho,
                                permutation_p=perm_p(x, y, 100000, rng),
                                bootstrap_ci95_lo=lo, bootstrap_ci95_hi=hi,
                                loo_rho_min=min(loo), loo_rho_max=max(loo),
                                malignant_label_source=g.malignant_label_source.iloc[0],
                                unit="patient", analysis="exploratory fixed-gene check"))
    patients = pd.concat(out, ignore_index=True)
    assoc = pd.DataFrame(results)
    for cohort, idx in assoc.groupby("cohort").groups.items():
        p = assoc.loc[idx, "permutation_p"].to_numpy()
        order = np.argsort(p)
        q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
        arr = np.empty(len(p)); arr[order] = np.clip(q, 0, 1)
        assoc.loc[idx, "fdr_bh_4_features"] = arr
    DATA.mkdir(parents=True, exist_ok=True)
    patients.to_csv(DATA / "patient_three_gene_score.tsv", sep="\t", index=False)
    assoc.to_csv(DATA / "three_gene_associations.tsv", sep="\t", index=False)
    (DATA / "run.json").write_text(json.dumps({"genes": GENES, "n_permutations": 100000,
            "n_bootstrap": 10000, "seed": SEED, "correlation_unit": "patient",
            "multiple_testing": "BH within cohort across signature and three individual genes",
            "cohorts_not_pooled": True}, indent=2), encoding="utf-8")
    print(assoc.to_string(index=False))


if __name__ == "__main__":
    main()
