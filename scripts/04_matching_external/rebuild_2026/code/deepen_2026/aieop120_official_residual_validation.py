"""AIEOP 120: official clinTALL published probabilities + frozen residual.

Assignment rule is the authors' figure rule, unchanged:
    predicted_subtype = argmax(published subtype probabilities)
No post-hoc probability cutoff is created after seeing residuals.
Labels are classifier-assigned, not genomic confirmation.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026")
PROB = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\tools\clinTALL\figures\figures\res\predictions_probs_BFM_italy.csv")
RESID = ROOT / "data/deepen_2026/a1_external_residual/aieop_patient_residual.tsv"
DISC = ROOT / "data/deepen_2026/a1_external_residual/polonen_discovery_subtype_medians.tsv"
S9 = ROOT / "data/deepen_2026/clintall/aieop_s9_labels.tsv"
OUT = ROOT / "data/deepen_2026/aieop120_official"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY = ["BCL11B", "TLX3", "ETP-like"]
MIN_N = 5
N_BOOT = 5000
SEED = 20260927


def norm_sub(x):
    if pd.isna(x):
        return ""
    s = str(x).replace("prob_", "").strip()
    compact = "".join(ch for ch in s.upper() if ch.isalnum())
    if compact.startswith("TAL1") and "DP" in compact:
        return "TAL1 DP-like"
    if compact.startswith("TAL1"):
        return "TAL1 αβ-like"
    if compact.startswith("LMO2"):
        return "LMO2 γδ-like"
    if "STAG2" in compact:
        return "STAG2&LMO2"
    if compact.replace("/", "&") == "STAG2LMO2":
        return "STAG2&LMO2"
    return s


def boot_median(x, rng):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan, np.nan)
    if x.size == 1:
        return (float(x[0]), np.nan, np.nan)
    draws = np.array([np.median(rng.choice(x, size=x.size, replace=True)) for _ in range(N_BOOT)])
    return float(np.median(x)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def boot_spearman(a, b, rng):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    rho = pd.Series(a).corr(pd.Series(b), method="spearman")
    draws = []
    n = len(a)
    for _ in range(N_BOOT):
        i = rng.choice(n, size=n, replace=True)
        draws.append(pd.Series(a[i]).corr(pd.Series(b[i]), method="spearman"))
    draws = np.array(draws, dtype=float)
    draws = draws[np.isfinite(draws)]
    if draws.size == 0:
        return float(rho) if pd.notna(rho) else np.nan, np.nan, np.nan
    return float(rho), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def main():
    probs = pd.read_csv(PROB, index_col=0)
    probs.index = probs.index.astype(str).str.strip()
    if probs.shape[0] != 120:
        raise ValueError(f"expected 120 AIEOP rows, got {probs.shape[0]}")
    pred = probs.idxmax(axis=1).map(norm_sub)
    maxp = probs.max(axis=1)
    labels = pd.DataFrame({
        "sample_id": pred.index,
        "predicted_subtype": pred.values,
        "max_probability": maxp.values,
        "n_subtype_classes": probs.shape[1],
        "assignment_rule": "official_published_argmax_no_extra_cutoff",
        "label_source": "clinTALL_figures_predictions_probs_BFM_italy",
        "label_status": "classifier-assigned",
    })
    labels.to_csv(OUT / "aieop120_official_labels.tsv", sep="\t", index=False)

    s9 = pd.read_csv(S9, sep="\t")
    s9["sample_id"] = s9["sample_id"].astype(str).str.strip()
    s9["predicted_subtype"] = s9["predicted_subtype"].map(norm_sub)
    chk = labels.merge(s9[["sample_id", "predicted_subtype", "fusion"]], on="sample_id", how="inner", suffixes=("_official120", "_s9"))
    chk["exact_match"] = chk["predicted_subtype_official120"] == chk["predicted_subtype_s9"]
    chk.to_csv(OUT / "aieop120_vs_s9.tsv", sep="\t", index=False)

    resid = pd.read_csv(RESID, sep="\t")
    resid["sample_id"] = resid["sample_id"].astype(str).str.strip()
    keep = ["sample_id", "balance", "dev", "expected", "residual"]
    dat = resid[keep].merge(labels, on="sample_id", how="inner")
    if len(dat) != 120:
        raise ValueError(f"residual/label merge lost samples: {len(dat)}")
    dat.to_csv(OUT / "aieop120_residual_with_official_labels.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    rows = []
    for sub, g in dat.groupby("predicted_subtype"):
        med, lo, hi = boot_median(g["residual"].to_numpy(), rng)
        rows.append({
            "predicted_subtype": sub,
            "n": int(len(g)),
            "residual_median": med,
            "residual_median_lo": lo,
            "residual_median_hi": hi,
            "mean_max_probability": float(g["max_probability"].mean()),
            "min_max_probability": float(g["max_probability"].min()),
        })
    aie_med = pd.DataFrame(rows).sort_values("residual_median")
    disc = pd.read_csv(DISC, sep="\t")
    disc["subtype"] = disc["subtype"].map(norm_sub)
    merged = aie_med.merge(disc, left_on="predicted_subtype", right_on="subtype", how="left", suffixes=("_aieop", "_polonen"))
    if "residual_median" in merged.columns and "residual_median_polonen" not in merged.columns:
        merged = merged.rename(columns={"residual_median": "residual_median_aieop"})
    if "n" in merged.columns and "n_aieop" not in merged.columns:
        merged = merged.rename(columns={"n": "n_aieop"})
    aie_n = "n_aieop" if "n_aieop" in merged.columns else "n"
    aie_med_col = "residual_median_aieop" if "residual_median_aieop" in merged.columns else "residual_median"
    pol_med_col = "residual_median_polonen" if "residual_median_polonen" in merged.columns else "residual_median"
    merged["same_sign"] = np.sign(merged[aie_med_col]) == np.sign(merged[pol_med_col])
    merged.to_csv(OUT / "aieop120_subtype_residual_medians.tsv", sep="\t", index=False)

    usable = merged.loc[merged[aie_n] >= MIN_N].copy()
    rho, rho_lo, rho_hi = boot_spearman(
        usable[pol_med_col].to_numpy(),
        usable[aie_med_col].to_numpy(),
        rng,
    )
    sign_conc = float(usable["same_sign"].mean()) if len(usable) else np.nan

    anchors = {}
    for name in PRIMARY:
        g = dat.loc[dat["predicted_subtype"] == name, "residual"]
        med, lo, hi = boot_median(g.to_numpy(), rng)
        anchors[name] = {
            "n": int(g.size),
            "median": med,
            "lo": lo,
            "hi": hi,
            "direction": "negative" if pd.notna(med) and med < 0 else ("positive" if pd.notna(med) and med > 0 else "missing"),
        }

    bcl_ok = anchors["BCL11B"]["n"] >= 1 and anchors["BCL11B"]["median"] < 0
    tlx_ok = anchors["TLX3"]["n"] >= 1 and anchors["TLX3"]["median"] > 0
    poles_ok = bcl_ok and tlx_ok
    poles_fail = (
        (anchors["BCL11B"]["n"] >= 1 and anchors["BCL11B"]["median"] > 0)
        or (anchors["TLX3"]["n"] >= 1 and anchors["TLX3"]["median"] < 0)
    )
    rank_strong = pd.notna(rho) and rho > 0 and pd.notna(rho_lo) and rho_lo > 0
    rank_near0 = (not pd.notna(rho)) or abs(rho) < 0.1

    if poles_ok and rank_strong:
        grade = "Strong"
    elif poles_ok:
        grade = "Moderate"
    elif poles_fail and rank_near0:
        grade = "Fail"
    else:
        grade = "Indeterminate"

    summary = pd.DataFrame([
        {"item": "n_aieop", "value": len(dat)},
        {"item": "n_official_classes", "value": probs.shape[1]},
        {"item": "assignment_rule", "value": "official_published_argmax_no_extra_cutoff"},
        {"item": "s9_n", "value": len(chk)},
        {"item": "s9_exact_match", "value": int(chk["exact_match"].sum())},
        {"item": "s9_exact_match_frac", "value": float(chk["exact_match"].mean()) if len(chk) else None},
        {"item": "BCL11B_n", "value": anchors["BCL11B"]["n"]},
        {"item": "BCL11B_median_residual", "value": anchors["BCL11B"]["median"]},
        {"item": "BCL11B_median_lo", "value": anchors["BCL11B"]["lo"]},
        {"item": "BCL11B_median_hi", "value": anchors["BCL11B"]["hi"]},
        {"item": "TLX3_n", "value": anchors["TLX3"]["n"]},
        {"item": "TLX3_median_residual", "value": anchors["TLX3"]["median"]},
        {"item": "TLX3_median_lo", "value": anchors["TLX3"]["lo"]},
        {"item": "TLX3_median_hi", "value": anchors["TLX3"]["hi"]},
        {"item": "ETP-like_n", "value": anchors["ETP-like"]["n"]},
        {"item": "ETP-like_median_residual", "value": anchors["ETP-like"]["median"]},
        {"item": "ETP-like_median_lo", "value": anchors["ETP-like"]["lo"]},
        {"item": "ETP-like_median_hi", "value": anchors["ETP-like"]["hi"]},
        {"item": "n_subtypes_ge5", "value": int(len(usable))},
        {"item": "subtypes_ge5", "value": ",".join(usable["predicted_subtype"].astype(str))},
        {"item": "spearman_median_residual", "value": rho},
        {"item": "spearman_lo", "value": rho_lo},
        {"item": "spearman_hi", "value": rho_hi},
        {"item": "sign_concordance_ge5", "value": sign_conc},
        {"item": "poles_ok", "value": poles_ok},
        {"item": "rank_CI_excludes_0", "value": rank_strong},
        {"item": "prelocked_grade", "value": grade},
        {"item": "grade_rule", "value": "Strong=BCL11B<0 and TLX3>0 and Spearman CI>0; Moderate=both poles repeat but rank CI wide; Fail=pole reverse and rank~0"},
        {"item": "label_status", "value": "classifier-assigned; not genomic confirmation"},
        {"item": "confidence_rule", "value": "authors publish probabilities and use argmax; no reject cutoff applied"},
    ])
    summary.to_csv(OUT / "aieop120_validation_summary.tsv", sep="\t", index=False)

    print(labels["predicted_subtype"].value_counts().to_string())
    print(summary.to_string(index=False))
    print("S9 mismatches:")
    print(chk.loc[~chk["exact_match"], ["sample_id", "predicted_subtype_official120", "predicted_subtype_s9", "fusion"]].to_string(index=False))
    print("usable n>=5:")
    print(usable[["predicted_subtype", aie_n, aie_med_col, pol_med_col, "same_sign"]].to_string(index=False))
    print("grade", grade)


if __name__ == "__main__":
    main()
