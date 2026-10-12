#!/usr/bin/env python3
"""Nested Cox C-index / LRT: does ZEB1 add discrimination beyond ETP+subtype+MRD?

Does not overwrite SU/MC/D or rewrite existing Cox HR numbers.
Uses cox_analysis_table.tsv (same table as run_cox.py).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import chi2

OUT = Path(r"d:\_bioinformation\ZEB1\02_analyses\03_q1_bioinfo_paper")
TABLE = OUT / "cox_analysis_table.tsv"
SEED = 42
N_BOOT = 200

CLIN_MRD = "age_years + C(sex) + log_wbc + C(protocol) + C(etp3) + C(subtype_g) + log_mrd"
CLIN_NO_MRD = "age_years + C(sex) + log_wbc + C(protocol) + C(etp3) + C(subtype_g)"
NONETP_MRD = "age_years + C(sex) + log_wbc + C(protocol) + C(subtype_g) + log_mrd"
NUM_COLS = ["ZEB1", "age_years", "log_wbc", "log_mrd"]


def fit_pair(df: pd.DataFrame, duration: str, event: str, reduced: str, full: str) -> dict:
    needed = [c for c in NUM_COLS if c in (reduced + " " + full)]
    sub = df.dropna(subset=[duration, event] + needed).copy()
    cph_r = CoxPHFitter()
    cph_f = CoxPHFitter()
    cph_r.fit(sub, duration_col=duration, event_col=event, formula=reduced, show_progress=False)
    cph_f.fit(sub, duration_col=duration, event_col=event, formula=full, show_progress=False)
    ll_r = float(cph_r.log_likelihood_)
    ll_f = float(cph_f.log_likelihood_)
    lr = 2.0 * (ll_f - ll_r)
    df_diff = int(cph_f.params_.shape[0] - cph_r.params_.shape[0])
    p_lrt = float(chi2.sf(max(lr, 0.0), df_diff)) if df_diff > 0 else None
    c_r = float(cph_r.concordance_index_)
    c_f = float(cph_f.concordance_index_)
    zeb = None
    if "ZEB1" in cph_f.params_.index:
        s = cph_f.summary.loc["ZEB1"]
        zeb = {
            "hr": float(s["exp(coef)"]),
            "p_wald": float(s["p"]),
            "lo": float(s["exp(coef) lower 95%"]),
            "hi": float(s["exp(coef) upper 95%"]),
        }
    return {
        "n": int(len(sub)),
        "n_event": int(sub[event].sum()),
        "reduced": reduced,
        "full": full,
        "C_reduced": c_r,
        "C_full": c_f,
        "delta_C": c_f - c_r,
        "AIC_reduced": float(cph_r.AIC_partial_),
        "AIC_full": float(cph_f.AIC_partial_),
        "LL_reduced": ll_r,
        "LL_full": ll_f,
        "LRT_stat": lr,
        "LRT_df": df_diff,
        "LRT_p": p_lrt,
        "ZEB1": zeb,
    }


def bootstrap_delta(
    df: pd.DataFrame,
    duration: str,
    event: str,
    reduced: str,
    full: str,
    n_boot: int,
    seed: int,
) -> dict:
    needed = [c for c in NUM_COLS if c in (reduced + " " + full)]
    sub = df.dropna(subset=[duration, event] + needed).reset_index(drop=True)
    rng = np.random.default_rng(seed)
    deltas = []
    n_fail = 0
    n = len(sub)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boot = sub.iloc[idx].copy()
        try:
            cph_r = CoxPHFitter()
            cph_f = CoxPHFitter()
            cph_r.fit(boot, duration_col=duration, event_col=event, formula=reduced, show_progress=False)
            cph_f.fit(boot, duration_col=duration, event_col=event, formula=full, show_progress=False)
            deltas.append(float(cph_f.concordance_index_) - float(cph_r.concordance_index_))
        except Exception:
            n_fail += 1
    arr = np.array(deltas, dtype=np.float64)
    if len(arr) < 20:
        return {"n_boot_ok": int(len(arr)), "n_fail": n_fail, "error": "too few successful boots"}
    return {
        "n_boot": n_boot,
        "n_boot_ok": int(len(arr)),
        "n_fail": n_fail,
        "delta_C_mean": float(arr.mean()),
        "delta_C_median": float(np.median(arr)),
        "ci95_lo": float(np.quantile(arr, 0.025)),
        "ci95_hi": float(np.quantile(arr, 0.975)),
        "frac_delta_gt_0": float(np.mean(arr > 0)),
        "frac_delta_ge_0_01": float(np.mean(arr >= 0.01)),
    }


def main() -> None:
    df = pd.read_csv(TABLE, sep="\t")
    pairs = {
        "OS_clinical_mrd_plus_ZEB1": ("os_time", "os_status", CLIN_MRD, CLIN_MRD + " + ZEB1", df),
        "OS_clinical_nomrd_plus_ZEB1": ("os_time", "os_status", CLIN_NO_MRD, CLIN_NO_MRD + " + ZEB1", df),
        "OS_mrd_increment_noZEB1": ("os_time", "os_status", CLIN_NO_MRD, CLIN_MRD, df),
        "EFS_clinical_mrd_plus_ZEB1": ("efs_time", "efs_status", CLIN_MRD, CLIN_MRD + " + ZEB1", df),
        "EFS_clinical_nomrd_plus_ZEB1": ("efs_time", "efs_status", CLIN_NO_MRD, CLIN_NO_MRD + " + ZEB1", df),
        "EFS_mrd_increment_noZEB1": ("efs_time", "efs_status", CLIN_NO_MRD, CLIN_MRD, df),
        "OS_nonetp_mrd_plus_ZEB1": (
            "os_time",
            "os_status",
            NONETP_MRD,
            NONETP_MRD + " + ZEB1",
            df[df["etp"] == "Non-ETP"].copy(),
        ),
        "EFS_nonetp_mrd_plus_ZEB1": (
            "efs_time",
            "efs_status",
            NONETP_MRD,
            NONETP_MRD + " + ZEB1",
            df[df["etp"] == "Non-ETP"].copy(),
        ),
    }

    results = {"n_table": int(len(df)), "n_boot": N_BOOT, "seed": SEED, "pairs": {}}
    for name, (dur, ev, red, full, use) in pairs.items():
        rec = fit_pair(use, dur, ev, red, full)
        print(
            f"{name} n={rec['n']} ev={rec['n_event']} "
            f"C {rec['C_reduced']:.4f}->{rec['C_full']:.4f} dC={rec['delta_C']:.4f} "
            f"LRT p={rec['LRT_p']}",
            flush=True,
        )
        results["pairs"][name] = rec

    boot_targets = [
        "OS_clinical_mrd_plus_ZEB1",
        "EFS_clinical_mrd_plus_ZEB1",
        "OS_mrd_increment_noZEB1",
    ]
    for name in boot_targets:
        dur, ev, red, full, use = pairs[name]
        print(f"bootstrap {name} n={N_BOOT} ...", flush=True)
        boot = bootstrap_delta(use, dur, ev, red, full, N_BOOT, SEED)
        results["pairs"][name]["bootstrap"] = boot
        print(json.dumps(boot), flush=True)

    (OUT / "NESTED_CINDEX.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    rows = []
    for name, rec in results["pairs"].items():
        row = {"pair": name}
        row.update({k: rec[k] for k in rec if k not in ("ZEB1", "bootstrap")})
        if rec.get("ZEB1"):
            row["ZEB1_hr"] = rec["ZEB1"]["hr"]
            row["ZEB1_p_wald"] = rec["ZEB1"]["p_wald"]
        boot = rec.get("bootstrap") or {}
        row["boot_dC_lo"] = boot.get("ci95_lo")
        row["boot_dC_hi"] = boot.get("ci95_hi")
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "nested_cindex.tsv", sep="\t", index=False)
    print("wrote NESTED_CINDEX.json", flush=True)


if __name__ == "__main__":
    main()
