#!/usr/bin/env python3
"""Patient-level ZEB-family state in Lim 2025 refractory pediatric T-ALL.

Frozen metrics: ZEB1, ZEB2, LMO2, balance = z(logCPM ZEB1) - z(logCPM ZEB2).
Unit of analysis is the patient, after malignant-blast pseudobulk.
The overlapping P058 H5AD is not read. L086 technical replicates are summed.
"""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.formula.api import mixedlm, ols

META = Path("/data-b/liangfuhua/projects/Human_ZEB1_IF_round1/data/cell_metadata.csv.gz")
SD1 = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/Lim_2025_Refractory_TALL/derived/SD1_Patient_cohort.tsv")
H5 = {
    "discovery": Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/Lim_2025_Refractory_TALL/cellxgene/fed94c35-c0ff-425a-abc8-4ecc30cc445d.h5ad"),
    "extension": Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/Lim_2025_Refractory_TALL/cellxgene/644d2ec8-d99c-4571-a94a-0f89cb33b251.h5ad"),
}
OUT = Path(os.environ.get("GATEA_OUT", "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/gateA_lim2025"))
GENES = {
    "ZEB1": "ENSG00000148516",
    "ZEB2": "ENSG00000169554",
    "LMO2": "ENSG00000135363",
    "ZBTB16": "ENSG00000109906",
}
MIN_BLASTS = 50
MIN_STATE = 20
SEED = 20260925
N_BOOT = 5000


def level1(subtype: str) -> str:
    s = str(subtype)
    if s in {"TLX1", "TLX3"}:
        return "TLX"
    if s == "TAL1":
        return "TAL"
    if s == "HOXA":
        return "HOXA"
    if s.startswith("LMO2"):
        return "LMO"
    return "other"


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan
    spooled = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if spooled == 0 or not np.isfinite(spooled):
        return np.nan
    d = (a.mean() - b.mean()) / spooled
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    return float(j * d)


def gene_columns(path: Path) -> dict[str, np.ndarray]:
    with h5py.File(path, "r") as f:
        ids = np.asarray(f["raw/var/_index"][:]).astype(str)
        g = f["raw/X"]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(int(x) for x in g.attrs["shape"])
        obs_names = np.asarray(f["obs/_index"][:]).astype(str)
    X = sp.csr_matrix((data, indices, indptr), shape=shape)
    out = {"obs_name": obs_names}
    for sym, ensg in GENES.items():
        hit = np.where(ids == ensg)[0]
        if hit.size != 1:
            raise RuntimeError(f"{path.name} {sym} {ensg} hits={hit.size}")
        out[sym] = np.asarray(X[:, int(hit[0])].todense()).ravel().astype(np.float64)
    return out


def load_cells() -> pd.DataFrame:
    usecols = [
        "unique_cell_key", "patient_id", "biological_sample_id", "cohort_cxg",
        "timepoint_std", "induction_outcome_std", "etp_status", "genomic_subtype",
        "age_bin", "sex_std", "is_malignant_main", "n_umi", "ZEB1_umi",
        "author_state_label", "disease_std",
    ]
    meta = pd.read_csv(META, usecols=usecols)
    frames = []
    for cohort, path in H5.items():
        print(f"[genes] {cohort}", flush=True)
        cols = gene_columns(path)
        part = pd.DataFrame(cols)
        part["cohort_file"] = cohort
        frames.append(part)
    genes = pd.concat(frames, ignore_index=True)
    # ZEB1 from the same raw matrix; metadata copy is only a cross-check
    merged = meta.merge(genes, left_on="unique_cell_key", right_on="obs_name", how="inner", validate="one_to_one")
    if len(merged) != len(meta):
        raise RuntimeError(f"cell join {len(merged)} != metadata {len(meta)}")
    delta = np.nanmax(np.abs(merged["ZEB1_umi"].to_numpy() - merged["ZEB1"]))
    if delta > 1e-3:
        raise RuntimeError(f"ZEB1 raw mismatch max abs {delta}")
    return merged


def pseudobulk(cells: pd.DataFrame, mask: np.ndarray, key: str) -> pd.DataFrame:
    sub = cells.loc[mask].copy()
    rows = []
    for sid, g in sub.groupby(key, sort=False):
        lib = float(g["n_umi"].sum())
        rec = {
            "biological_sample_id": sid,
            "n_cells": int(len(g)),
            "lib_size": lib,
            "patient_id": g["patient_id"].iloc[0],
            "cohort": g["cohort_cxg"].iloc[0],
            "timepoint": g["timepoint_std"].iloc[0],
            "induction": g["induction_outcome_std"].iloc[0],
            "etp_status": g["etp_status"].iloc[0],
            "genomic_subtype": g["genomic_subtype"].iloc[0],
            "age_bin": g["age_bin"].iloc[0],
            "sex": g["sex_std"].iloc[0],
            "disease": g["disease_std"].iloc[0],
        }
        rec["subtype_l1"] = level1(rec["genomic_subtype"])
        for sym in GENES:
            rec[f"{sym}_umi"] = float(g[sym].sum())
            rec[f"{sym}_logcpm"] = float(np.log2(rec[f"{sym}_umi"] / lib * 1e6 + 1)) if lib > 0 else np.nan
        rec["zbtb16_pos_cells"] = int((g["ZBTB16"] > 0).sum())
        rows.append(rec)
    return pd.DataFrame(rows)


def add_balance(df: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    mu = {g: float(ref[f"{g}_logcpm"].mean()) for g in ("ZEB1", "ZEB2")}
    sd = {g: float(ref[f"{g}_logcpm"].std(ddof=1)) for g in ("ZEB1", "ZEB2")}
    out["z_ZEB1"] = (out["ZEB1_logcpm"] - mu["ZEB1"]) / sd["ZEB1"]
    out["z_ZEB2"] = (out["ZEB2_logcpm"] - mu["ZEB2"]) / sd["ZEB2"]
    out["balance"] = out["z_ZEB1"] - out["z_ZEB2"]
    out.attrs["z_mean"] = mu
    out.attrs["z_sd"] = sd
    return out


def contrast(a: np.ndarray, b: np.ndarray, label_a: str, label_b: str) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    diff = float(a.mean() - b.mean()) if len(a) and len(b) else np.nan
    welch = stats.ttest_ind(a, b, equal_var=False, alternative="two-sided") if len(a) > 1 and len(b) > 1 else None
    mann = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else None
    rng = np.random.default_rng(SEED)
    boots = np.empty(N_BOOT)
    for i in range(N_BOOT):
        aa = a[rng.integers(0, len(a), len(a))]
        bb = b[rng.integers(0, len(b), len(b))]
        boots[i] = aa.mean() - bb.mean()
    # leave-one-out on the pooled patient list is done by caller for paired ids;
    # here LOO within each arm, recomputing the mean difference.
    signs = []
    for i in range(len(a)):
        signs.append(np.sign(np.delete(a, i).mean() - b.mean()))
    for i in range(len(b)):
        signs.append(np.sign(a.mean() - np.delete(b, i).mean()))
    signs = np.asarray(signs)
    base = np.sign(diff) if diff != 0 else 0
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "label_a": label_a,
        "label_b": label_b,
        "mean_a": float(a.mean()) if len(a) else np.nan,
        "mean_b": float(b.mean()) if len(b) else np.nan,
        "diff_a_minus_b": diff,
        "hedges_g_a_minus_b": hedges_g(a, b),
        "welch_p": float(welch.pvalue) if welch is not None else np.nan,
        "mannwhitney_p": float(mann.pvalue) if mann is not None else np.nan,
        "boot_ci95_lo": float(np.quantile(boots, 0.025)),
        "boot_ci95_hi": float(np.quantile(boots, 0.975)),
        "boot_same_sign_frac": float(np.mean(np.sign(boots) == base)) if base != 0 else np.nan,
        "loo_n": int(len(signs)),
        "loo_sign_flips": int(np.sum(signs != base)) if base != 0 else np.nan,
    }


def one_sample(x: np.ndarray) -> dict:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return {"n": int(len(x)), "mean": float(x.mean()) if len(x) else np.nan}
    wil = stats.wilcoxon(x, alternative="two-sided", zero_method="wilcox")
    t = stats.ttest_1samp(x, 0.0)
    rng = np.random.default_rng(SEED)
    boots = np.array([x[rng.integers(0, len(x), len(x))].mean() for _ in range(N_BOOT)])
    base = np.sign(x.mean())
    flips = 0
    for i in range(len(x)):
        if np.sign(np.delete(x, i).mean()) != base:
            flips += 1
    return {
        "n": int(len(x)),
        "mean": float(x.mean()),
        "median": float(np.median(x)),
        "welch_p": float(t.pvalue),
        "wilcoxon_p": float(wil.pvalue),
        "boot_ci95_lo": float(np.quantile(boots, 0.025)),
        "boot_ci95_hi": float(np.quantile(boots, 0.975)),
        "boot_same_sign_frac": float(np.mean(np.sign(boots) == base)) if base != 0 else np.nan,
        "loo_sign_flips": int(flips),
    }


def ols_sensitivity(df: pd.DataFrame, y: str) -> pd.DataFrame:
    d = df.copy()
    d["refractory"] = (d["induction"] == "IF").astype(int)
    rows = []
    for formula, tag in [
        (f"{y} ~ refractory", "unadjusted"),
        (f"{y} ~ refractory + C(subtype_l1)", "plus_subtype_l1"),
        (f"{y} ~ refractory + C(age_bin)", "plus_age_bin"),
        (f"{y} ~ refractory + C(etp_status)", "plus_etp"),
    ]:
        try:
            fit = ols(formula, data=d).fit()
            ci = fit.conf_int().loc["refractory"]
            rows.append({
                "y": y,
                "model": tag,
                "n": int(fit.nobs),
                "coef_refractory": float(fit.params["refractory"]),
                "ci95_lo": float(ci[0]),
                "ci95_hi": float(ci[1]),
                "p": float(fit.pvalues["refractory"]),
            })
        except Exception as exc:  # noqa: BLE001
            rows.append({"y": y, "model": tag, "n": len(d), "error": str(exc)})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells = load_cells()
    mal = cells["is_malignant_main"].to_numpy()
    print(f"[cells] {len(cells)} malignant {int(mal.sum())}", flush=True)
    pb = pseudobulk(cells, mal, "biological_sample_id")
    pb.to_csv(OUT / "pseudobulk_all_timepoints.tsv", sep="\t", index=False)

    eligible = pb[
        pb["timepoint"].isin(["Day0", "Day28"])
        & pb["induction"].isin(["IF", "responsive"])
        & (pb["n_cells"] >= MIN_BLASTS)
    ].copy()
    # one biological sample per patient x timepoint is already the key
    ref = eligible.copy()
    eligible = add_balance(eligible, ref)
    mu, sd = {g: float(ref[f"{g}_logcpm"].mean()) for g in ("ZEB1", "ZEB2")}, {
        g: float(ref[f"{g}_logcpm"].std(ddof=1)) for g in ("ZEB1", "ZEB2")
    }
    pd.DataFrame([
        {"gene": g, "logcpm_mean": mu[g], "logcpm_sd": sd[g], "n_reference_samples": len(ref)}
        for g in ("ZEB1", "ZEB2")
    ]).to_csv(OUT / "balance_reference.tsv", sep="\t", index=False)
    eligible.to_csv(OUT / "eligible_patient_timepoint.tsv", sep="\t", index=False)

    day0 = eligible[eligible["timepoint"] == "Day0"].copy()
    if day0["patient_id"].duplicated().any():
        raise RuntimeError("duplicate Day0 patients")
    metrics = ["ZEB1_logcpm", "ZEB2_logcpm", "LMO2_logcpm", "balance"]
    q1_rows = []
    for cohort_name, sub in [("combined", day0), ("discovery", day0[day0.cohort == "discovery"]), ("extension", day0[day0.cohort == "extension"])]:
        for metric in metrics:
            a = sub.loc[sub.induction == "IF", metric].to_numpy()
            b = sub.loc[sub.induction == "responsive", metric].to_numpy()
            rec = contrast(a, b, "IF", "Responsive")
            rec.update({"cohort": cohort_name, "metric": metric, "question": "Q1_day0"})
            q1_rows.append(rec)
    q1 = pd.DataFrame(q1_rows)
    q1.to_csv(OUT / "q1_baseline.tsv", sep="\t", index=False)

    sens = pd.concat([ols_sensitivity(day0, m) for m in metrics], ignore_index=True)
    sens.to_csv(OUT / "q1_sensitivity_ols.tsv", sep="\t", index=False)

    # paired Day0-Day28
    d0 = eligible[eligible.timepoint == "Day0"].set_index("patient_id")
    d28 = eligible[eligible.timepoint == "Day28"].set_index("patient_id")
    both = d0.index.intersection(d28.index)
    pair_rows = []
    long_rows = []
    for pid in both:
        r0, r28 = d0.loc[pid], d28.loc[pid]
        pair_rows.append({
            "patient_id": pid,
            "cohort": r0["cohort"],
            "induction": r0["induction"],
            "subtype_l1": r0["subtype_l1"],
            "genomic_subtype": r0["genomic_subtype"],
            "n_day0": int(r0["n_cells"]),
            "n_day28": int(r28["n_cells"]),
            **{f"d_{m}": float(r28[m] - r0[m]) for m in metrics},
            **{f"day0_{m}": float(r0[m]) for m in metrics},
            **{f"day28_{m}": float(r28[m]) for m in metrics},
        })
        for tp, rr in (("Day0", r0), ("Day28", r28)):
            long_rows.append({
                "patient_id": pid,
                "time": 1 if tp == "Day28" else 0,
                "refractory": int(r0["induction"] == "IF"),
                "cohort": r0["cohort"],
                **{m: float(rr[m]) for m in metrics},
            })
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(OUT / "q2_pairs.tsv", sep="\t", index=False)
    q2_rows = []
    for metric in metrics:
        col = f"d_{metric}"
        for name, sub in [
            ("all_paired", pairs),
            ("IF", pairs[pairs.induction == "IF"]),
            ("responsive", pairs[pairs.induction == "responsive"]),
        ]:
            rec = one_sample(sub[col].to_numpy())
            rec.update({"group": name, "metric": metric, "question": "Q2_delta"})
            q2_rows.append(rec)
        if pairs["induction"].nunique() == 2 and len(pairs) >= 6:
            a = pairs.loc[pairs.induction == "IF", col].to_numpy()
            b = pairs.loc[pairs.induction == "responsive", col].to_numpy()
            if len(a) >= 2 and len(b) >= 2:
                rec = contrast(a, b, "IF_delta", "Responsive_delta")
                rec.update({"group": "IF_minus_responsive_delta", "metric": metric, "question": "Q2_interaction_delta"})
                q2_rows.append(rec)
    long = pd.DataFrame(long_rows)
    mix_rows = []
    for metric in metrics:
        try:
            fit = mixedlm(f"{metric} ~ time * refractory", long, groups=long["patient_id"]).fit(reml=False, method="lbfgs")
            for term in ("time", "refractory", "time:refractory"):
                ci = fit.conf_int().loc[term]
                mix_rows.append({
                    "metric": metric,
                    "term": term,
                    "coef": float(fit.params[term]),
                    "ci95_lo": float(ci[0]),
                    "ci95_hi": float(ci[1]),
                    "p": float(fit.pvalues[term]),
                    "n_obs": int(fit.nobs),
                    "converged": bool(fit.converged),
                })
        except Exception as exc:  # noqa: BLE001
            mix_rows.append({"metric": metric, "error": str(exc)})
    pd.DataFrame(q2_rows).to_csv(OUT / "q2_longitudinal.tsv", sep="\t", index=False)
    pd.DataFrame(mix_rows).to_csv(OUT / "q2_mixed_model.tsv", sep="\t", index=False)

    # Q3: ZBTB16+ vs ZBTB16- within the same patient-timepoint, Day0 primary
    pos = pseudobulk(cells, mal & (cells["ZBTB16"].to_numpy() > 0), "biological_sample_id")
    neg = pseudobulk(cells, mal & (cells["ZBTB16"].to_numpy() <= 0), "biological_sample_id")
    pos = pos[pos.n_cells >= MIN_STATE]
    neg = neg[neg.n_cells >= MIN_STATE]
    pos = add_balance(pos, ref)
    neg = add_balance(neg, ref)
    state = pos.merge(neg, on="biological_sample_id", suffixes=("_pos", "_neg"))
    state = state[state["timepoint_pos"].isin(["Day0", "Day28"]) & state["induction_pos"].isin(["IF", "responsive"])]
    for m in metrics:
        state[f"d_{m}"] = state[f"{m}_pos"] - state[f"{m}_neg"]
    state_day0 = state[state.timepoint_pos == "Day0"].copy()
    q3_rows = []
    for metric in metrics:
        rec = one_sample(state_day0[f"d_{metric}"].to_numpy())
        rec.update({"group": "all_day0_pos_minus_neg", "metric": metric})
        q3_rows.append(rec)
        for outcome in ("IF", "responsive"):
            sub = state_day0[state_day0.induction_pos == outcome]
            rec = one_sample(sub[f"d_{metric}"].to_numpy())
            rec.update({"group": outcome, "metric": metric})
            q3_rows.append(rec)
        a = state_day0.loc[state_day0.induction_pos == "IF", f"d_{metric}"].to_numpy()
        b = state_day0.loc[state_day0.induction_pos == "responsive", f"d_{metric}"].to_numpy()
        if len(a) >= 2 and len(b) >= 2:
            rec = contrast(a, b, "IF", "Responsive")
            rec.update({"group": "IF_minus_responsive", "metric": metric})
            q3_rows.append(rec)
    state_day0.to_csv(OUT / "q3_state_day0_pairs.tsv", sep="\t", index=False)
    pd.DataFrame(q3_rows).to_csv(OUT / "q3_state.tsv", sep="\t", index=False)

    # inventory
    inv = {
        "n_cells_joined": int(len(cells)),
        "n_malignant": int(mal.sum()),
        "n_patients_meta": int(cells.patient_id.nunique()),
        "n_biological_samples_malignant": int(pb.biological_sample_id.nunique()),
        "n_eligible_day0": int((eligible.timepoint == "Day0").sum()),
        "n_eligible_day28": int((eligible.timepoint == "Day28").sum()),
        "n_paired": int(len(pairs)),
        "n_paired_IF": int((pairs.induction == "IF").sum()) if len(pairs) else 0,
        "n_paired_responsive": int((pairs.induction == "responsive").sum()) if len(pairs) else 0,
        "min_blasts": MIN_BLASTS,
        "min_state_cells": MIN_STATE,
        "p058_overlap_h5ad": "excluded",
        "l086": "R1+R2 summed before logCPM",
        "matrix": "raw/X integer UMI",
        "logcpm": "log2(UMI/lib_size*1e6+1)",
        "balance_reference": "z-score of logCPM across eligible Day0 and Day28 samples",
    }
    pd.Series(inv).to_csv(OUT / "inventory.tsv", sep="\t", header=False)
    print(pd.Series(inv).to_string())
    print(q1[q1.cohort == "combined"][["metric", "n_a", "n_b", "diff_a_minus_b", "hedges_g_a_minus_b", "boot_ci95_lo", "boot_ci95_hi", "boot_same_sign_frac", "loo_sign_flips", "welch_p"]].to_string(index=False))
    print("DONE", OUT)


if __name__ == "__main__":
    main()
