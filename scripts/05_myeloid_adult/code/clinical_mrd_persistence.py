#!/usr/bin/env python3
"""Exploratory day-15-positive MRD persistence analysis from public source table.

The 0.01% threshold was confirmed in Lee et al., PMID 36604538.
The source day-15 column includes day 19 for Total XV.
All predictors are measured at diagnosis or day 15, before the late endpoint.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

THRESHOLD = 0.01  # percent; source-defined MRD negativity boundary


def loo_predict(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.full(len(y), np.nan)
    for train, test in LeaveOneOut().split(x):
        model = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
        model.fit(x[train], y[train])
        out[test[0]] = model.predict_proba(x[test])[0, 1]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    raw = pd.read_csv(args.input, sep="\t")
    required = ["sample", "ZEB1_log", "Day 15 MRD (%)", "Day 42 or 46 MRD (%)"]
    if any(c not in raw.columns for c in required) or raw["sample"].duplicated().any():
        raise ValueError("source-table schema or unique sample ID failed")
    early = "Day 15 MRD (%)"
    late = "Day 42 or 46 MRD (%)"
    eligible = raw[raw[early].notna() & raw[late].notna() & raw[early].ge(THRESHOLD)].copy()
    eligible["late_positive"] = eligible[late].ge(THRESHOLD).astype(int)
    eligible["log10_early_mrd_pct"] = np.log10(eligible[early])
    if len(eligible) != 66 or int(eligible["late_positive"].sum()) != 18:
        raise ValueError("source cohort changed; reassess design before using this script")
    if eligible[required].isna().any().any():
        raise ValueError("missing primary input in eligible cohort")
    y = eligible["late_positive"].to_numpy(dtype=int)
    x0 = eligible[["log10_early_mrd_pct"]].to_numpy(dtype=float)
    x1 = eligible[["log10_early_mrd_pct", "ZEB1_log"]].to_numpy(dtype=float)
    p0 = loo_predict(x0, y)
    p1 = loo_predict(x1, y)
    eligible["loo_prob_early_mrd"] = p0
    eligible["loo_prob_early_mrd_plus_zeb1"] = p1
    # Descriptive maximum-likelihood odds ratio, per 1 SD diagnosis ZEB1.
    z = (eligible["ZEB1_log"] - eligible["ZEB1_log"].mean()) / eligible["ZEB1_log"].std(ddof=1)
    design = sm.add_constant(pd.DataFrame({"log10_early_mrd_pct": eligible["log10_early_mrd_pct"],
                                            "ZEB1_per_SD": z}, index=eligible.index))
    fit = sm.Logit(y, design).fit(disp=False)
    ci = fit.conf_int().loc["ZEB1_per_SD"]
    rho_early, p_early = spearmanr(eligible["ZEB1_log"], eligible[early])
    rho_late, p_late = spearmanr(eligible["ZEB1_log"], eligible[late])
    summary = {
        "n_source": len(raw), "n_eligible": len(eligible), "n_late_positive": int(y.sum()),
        "threshold_percent": THRESHOLD,
        "unit": "patient", "endpoint": "day 42 or 46 MRD >= 0.01% among day 15 MRD >= 0.01%",
        "model": "day-15 log10 MRD percent, with or without diagnosis ZEB1_log",
        "loo_auc_baseline": float(roc_auc_score(y, p0)),
        "loo_auc_plus_zeb1": float(roc_auc_score(y, p1)),
        "loo_auc_delta": float(roc_auc_score(y, p1) - roc_auc_score(y, p0)),
        "loo_brier_baseline": float(brier_score_loss(y, p0)),
        "loo_brier_plus_zeb1": float(brier_score_loss(y, p1)),
        "zeb1_or_per_sd": float(np.exp(fit.params["ZEB1_per_SD"])),
        "zeb1_or_ci95_lo": float(np.exp(ci.iloc[0])),
        "zeb1_or_ci95_hi": float(np.exp(ci.iloc[1])),
        "zeb1_wald_p": float(fit.pvalues["ZEB1_per_SD"]),
        "zeb1_spearman_early_rho": float(rho_early),
        "zeb1_spearman_early_p": float(p_early),
        "zeb1_spearman_late_rho": float(rho_late),
        "zeb1_spearman_late_p": float(p_late),
        "validation_status": "internal leave-one-out only; no independent clinical validation",
    }
    eligible.to_csv(args.output / "paired_patients_predictions.tsv", sep="\t", index=False)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
