#!/usr/bin/env python3
"""Cox PH for ZEB1 in Pölönen T-ALL. Does not overwrite SU/MC/D."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, proportional_hazard_test

ROOT = Path(r"d:\_bioinformation\ZEB1")
OUT = ROOT / "02_analyses" / "03_q1_bioinfo_paper"
TPM = ROOT / "00_data" / "00_raw" / "TALL_X01_tpm.tsv"
CLIN = ROOT / "00_data" / "00_raw" / "polonen_syn54032669" / "audit_clinical_1309.tsv"
ZEB1 = "ENSG00000148516"


def extract_zeb1() -> pd.Series:
    with TPM.open(encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            gid = line.split("\t", 1)[0].split(".")[0]
            if gid == ZEB1:
                vals = [float(x) for x in line.rstrip("\n").split("\t")[1:]]
                return pd.Series(np.log2(np.array(vals, dtype=np.float64) + 1.0), index=samples, name="ZEB1")
    raise SystemExit("ZEB1 not found")


def group_subtype(s: pd.Series) -> pd.Series:
    counts = s.value_counts()
    keep = set(counts[counts >= 30].index)
    return s.where(s.isin(keep), "Other")


def prepare() -> pd.DataFrame:
    z = extract_zeb1()
    df = pd.read_csv(CLIN, sep="\t")
    df["sample_id"] = df["sample_id"].astype(str)
    df = df.merge(z.reset_index().rename(columns={"index": "sample_id"}), on="sample_id", how="inner")
    df["wbc"] = pd.to_numeric(df["wbc"], errors="coerce")
    df["mrd_d29"] = pd.to_numeric(df["mrd_d29"], errors="coerce")
    df["age_years"] = pd.to_numeric(df["age_years"], errors="coerce")
    df["log_wbc"] = np.log1p(df["wbc"])
    df["log_mrd"] = np.log1p(df["mrd_d29"].clip(lower=0))
    df["sex"] = df["sex"].astype(str)
    df["etp3"] = df["etp"].where(df["etp"].isin(["ETP", "Near-ETP", "Non-ETP"]), "Unknown")
    df["subtype_g"] = group_subtype(df["subtype"].astype(str))
    df["protocol"] = df["protocol"].astype(str)
    df["zeb1_sd"] = (df["ZEB1"] - df["ZEB1"].mean()) / df["ZEB1"].std(ddof=0)
    df["zeb1_tertile"] = pd.qcut(df["ZEB1"], 3, labels=["T1_low", "T2_mid", "T3_high"])
    for col in ["os_time", "efs_time"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df.loc[df[col] <= 0, col] = 0.5
    return df


def fit_cox(df: pd.DataFrame, duration: str, event: str, formula: str, tag: str) -> dict:
    sub = df.dropna(subset=[duration, event]).copy()
    used = [p.strip() for p in formula.replace("~", "").split("+")]
    # drop rows missing formula covariates that are numeric columns
    numeric_needed = [c for c in ["ZEB1", "zeb1_sd", "age_years", "log_wbc", "log_mrd"] if c in formula]
    sub = sub.dropna(subset=numeric_needed)
    cph = CoxPHFitter()
    cph.fit(sub, duration_col=duration, event_col=event, formula=formula, show_progress=False)
    summ = cph.summary.reset_index().rename(columns={"covariate": "term"})
    summ.to_csv(OUT / f"cox_{tag}.tsv", sep="\t", index=False)
    row = None
    for key in ("ZEB1", "zeb1_sd"):
        hit = summ[summ["covariate"] == key] if "covariate" in summ.columns else summ[summ.iloc[:, 0] == key]
        if "covariate" not in summ.columns:
            first = summ.columns[0]
            hit = summ[summ[first] == key]
        if len(hit):
            r = hit.iloc[0]
            row = {
                "term": key,
                "hr": float(r["exp(coef)"]),
                "hr_lo": float(r["exp(coef) lower 95%"]),
                "hr_hi": float(r["exp(coef) upper 95%"]),
                "p": float(r["p"]),
                "coef": float(r["coef"]),
            }
            break
    ph = None
    try:
        ph_res = proportional_hazard_test(cph, sub, time_transform="rank")
        ph = {"p_global": float(ph_res.summary["p"].min()) if len(ph_res.summary) else None}
        zeb = ph_res.summary.loc[ph_res.summary.index.astype(str).str.contains("ZEB1|zeb1")]
        if len(zeb):
            ph["p_ZEB1"] = float(zeb["p"].iloc[0])
    except Exception as exc:
        ph = {"error": str(exc)}
    return {
        "tag": tag,
        "n": int(len(sub)),
        "n_event": int(sub[event].sum()),
        "formula": formula,
        "concordance": float(cph.concordance_index_),
        "ZEB1": row,
        "ph": ph,
    }


def km_tertiles(df: pd.DataFrame, duration: str, event: str, tag: str) -> dict:
    sub = df.dropna(subset=[duration, event, "zeb1_tertile"]).copy()
    km = KaplanMeierFitter()
    rows = []
    for lab, g in sub.groupby("zeb1_tertile", observed=True):
        km.fit(g[duration], g[event], label=str(lab))
        sf = km.survival_function_.reset_index()
        sf.columns = ["time", "surv"]
        sf["tertile"] = str(lab)
        rows.append(sf)
        rows[-1]  # keep
    pd.concat(rows).to_csv(OUT / f"km_{tag}_curves.tsv", sep="\t", index=False)
    counts = sub.groupby("zeb1_tertile", observed=True)[event].agg(n="count", events="sum")
    counts.to_csv(OUT / f"km_{tag}_counts.tsv", sep="\t")
    res = multivariate_logrank_test(sub[duration], sub["zeb1_tertile"].astype(str), sub[event])
    return {
        "n": int(len(sub)),
        "logrank_p": float(res.p_value),
        "test_statistic": float(res.test_statistic),
        "counts": counts.reset_index().to_dict("records"),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = prepare()
    df.to_csv(OUT / "cox_analysis_table.tsv", sep="\t", index=False)

    formulas = {
        "uni": "ZEB1",
        "uni_sd": "zeb1_sd",
        "demo": "ZEB1 + age_years + C(sex) + log_wbc + C(protocol)",
        "etp": "ZEB1 + age_years + C(sex) + log_wbc + C(protocol) + C(etp3)",
        "subtype": "ZEB1 + age_years + C(sex) + log_wbc + C(protocol) + C(etp3) + C(subtype_g)",
        "mrd": "ZEB1 + age_years + C(sex) + log_wbc + C(protocol) + C(etp3) + C(subtype_g) + log_mrd",
        "nonetp_uni": "ZEB1",
        "nonetp_mrd": "ZEB1 + age_years + C(sex) + log_wbc + C(protocol) + C(subtype_g) + log_mrd",
    }

    results = {"n_total": int(len(df)), "zeb1_median": float(df["ZEB1"].median()), "models": {}}
    for outcome, dur, ev in [("OS", "os_time", "os_status"), ("EFS", "efs_time", "efs_status")]:
        results["models"][outcome] = {}
        for name, formula in formulas.items():
            use = df
            if name.startswith("nonetp"):
                use = df[df["etp"] == "Non-ETP"].copy()
            try:
                results["models"][outcome][name] = fit_cox(use, dur, ev, formula, f"{outcome}_{name}")
                print(f"{outcome} {name} n={results['models'][outcome][name]['n']} ZEB1={results['models'][outcome][name]['ZEB1']}", flush=True)
            except Exception as exc:
                results["models"][outcome][name] = {"error": str(exc), "formula": formula}
                print(f"{outcome} {name} FAIL {exc}", flush=True)
        results["models"][outcome]["km_tertile"] = km_tertiles(df, dur, ev, outcome)
        print(f"{outcome} KM p={results['models'][outcome]['km_tertile']['logrank_p']}", flush=True)

    (OUT / "COX_RECEIPT.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    # compact headline
    def grab(outcome, model):
        m = results["models"][outcome].get(model, {})
        z = m.get("ZEB1") or {}
        return {"n": m.get("n"), "n_event" : m.get("n_event"), "hr": z.get("hr"), "p": z.get("p"), "lo": z.get("hr_lo"), "hi": z.get("hr_hi")}

    headline = {
        "OS_uni": grab("OS", "uni"),
        "OS_etp": grab("OS", "etp"),
        "OS_subtype": grab("OS", "subtype"),
        "OS_mrd": grab("OS", "mrd"),
        "EFS_uni": grab("EFS", "uni"),
        "EFS_etp": grab("EFS", "etp"),
        "EFS_subtype": grab("EFS", "subtype"),
        "EFS_mrd": grab("EFS", "mrd"),
        "OS_nonetp_uni": grab("OS", "nonetp_uni"),
        "EFS_nonetp_uni": grab("EFS", "nonetp_uni"),
        "OS_km_p": results["models"]["OS"]["km_tertile"]["logrank_p"],
        "EFS_km_p": results["models"]["EFS"]["km_tertile"]["logrank_p"],
    }
    print(json.dumps(headline, indent=2))
    (OUT / "COX_HEADLINE.json").write_text(json.dumps(headline, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
