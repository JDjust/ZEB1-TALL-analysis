#!/usr/bin/env python3
"""Formal four-cohort ZEB1 / ZEB2 / LMO2 / balance subtype-vs-rest effects.

Exploratory relative to the original ZEB1-only Module 3 plan, but the
estimator matches 02b_analyze_independent.R: within-cohort subtype versus
the other labeled subtypes, Hedges g, DerSimonian-Laird, plus the REML
modified Hartung-Knapp sensitivity already used for TLX.

ratio_zeb is the within-cohort log2 balance. Cohort values in
independent_keygenes.tsv are already on a compressed scale (array intensity
or log expression), so balance = ZEB1 - ZEB2. That is log2(intensity ratio)
when both genes are log2, and log2((ZEB1+1)/(ZEB2+1)) when both are
log2(x+1). It is not recomputed from linear FPKM.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import mannwhitneyu, norm, t

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "modules" / "module3" / "tables" / "independent_keygenes.tsv"
OUT = ROOT / "data" / "validation" / "zeb_family_state"
SEED = 20260925
N_BOOT = 2000
SUBTYPES = ("HOXA", "TAL", "TLX", "Immature")
GENES = ("ZEB1", "ZEB2", "LMO2", "ratio_zeb")


def hedges_g(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return np.nan, np.nan
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if not np.isfinite(sp) or sp == 0:
        return 0.0, 0.0
    d = (a.mean() - b.mean()) / sp
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    g = j * d
    se = np.sqrt((n1 + n2) / (n1 * n2) + d**2 / (2 * (n1 + n2))) * j
    return float(g), float(se)


def fit_meta(y: np.ndarray, variance: np.ndarray, method: str) -> dict:
    y = np.asarray(y, float)
    variance = np.asarray(variance, float)
    k = len(y)
    if method == "DL_normal":
        w0 = 1 / variance
        mu0 = np.sum(w0 * y) / w0.sum()
        q0 = np.sum(w0 * (y - mu0) ** 2)
        c = w0.sum() - np.sum(w0**2) / w0.sum()
        tau2 = max(0.0, (q0 - k + 1) / c)
        q_fixed = float(q0)
    else:
        def objective(tau):
            w = 1 / (variance + tau)
            mu = np.sum(w * y) / w.sum()
            return 0.5 * (np.log(variance + tau).sum() + np.log(w.sum()) + np.sum(w * (y - mu) ** 2))

        upper = max(1.0, float(np.var(y) * 10))
        result = minimize_scalar(objective, bounds=(0, upper), method="bounded", options={"xatol": 1e-12})
        if not result.success:
            raise RuntimeError(result.message)
        tau2 = 0.0 if objective(0) <= result.fun else float(result.x)
        w0 = 1 / variance
        mu0 = np.sum(w0 * y) / w0.sum()
        q_fixed = float(np.sum(w0 * (y - mu0) ** 2))
    w = 1 / (variance + tau2)
    mu = float(np.sum(w * y) / w.sum())
    q = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    se = float(np.sqrt((max(1.0, q) if method == "REML_modified_HK" else 1.0) / w.sum()))
    cutoff = float(t.ppf(0.975, k - 1)) if method == "REML_modified_HK" else float(norm.ppf(0.975))
    p = float(2 * t.sf(abs(mu / se), k - 1)) if method == "REML_modified_HK" else float(2 * norm.sf(abs(mu / se)))
    i2 = max(0.0, (q_fixed - (k - 1)) / q_fixed) * 100 if q_fixed > 0 else 0.0
    return {
        "method": method,
        "k": k,
        "hedges_g": mu,
        "se": se,
        "tau2": tau2,
        "I2": i2,
        "Q": q_fixed,
        "ci95_lo": mu - cutoff * se,
        "ci95_hi": mu + cutoff * se,
        "p": p,
        "n_positive": int(np.sum(y > 0)),
        "n_negative": int(np.sum(y < 0)),
        "sign_agreement": f"{int(max(np.sum(y > 0), np.sum(y < 0)))}/{k}",
    }


def wilcoxon(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(SOURCE, sep="\t")
    tall = raw[raw["disease"].fillna("").str.contains("T-ALL", case=False)].copy()
    tall = tall[tall["level1"].isin(SUBTYPES)].copy()
    for gene in ("ZEB1", "ZEB2", "LMO2"):
        tall[gene] = pd.to_numeric(tall[gene], errors="coerce")
    tall = tall.dropna(subset=["ZEB1", "ZEB2", "LMO2"])
    med = tall.groupby("cohort")[["ZEB1", "ZEB2"]].median()
    if (med.max(axis=1) > 30).any():
        raise RuntimeError(f"unexpected linear-scale cohort:\n{med}")
    tall["ratio_zeb"] = tall["ZEB1"] - tall["ZEB2"]

    rows = []
    cohorts = sorted(tall["cohort"].unique())
    for cohort in cohorts:
        part = tall[tall["cohort"].eq(cohort)]
        for gene in GENES:
            for subtype in SUBTYPES:
                a = part.loc[part["level1"].eq(subtype), gene].to_numpy()
                b = part.loc[~part["level1"].eq(subtype), gene].to_numpy()
                g, se = hedges_g(a, b)
                rows.append({
                    "cohort": cohort,
                    "platform": part["platform"].iloc[0],
                    "gene": gene,
                    "subtype": subtype,
                    "n_subtype": int(np.isfinite(a).sum()),
                    "n_rest": int(np.isfinite(b).sum()),
                    "mean_subtype": float(np.nanmean(a)),
                    "mean_rest": float(np.nanmean(b)),
                    "hedges_g": g,
                    "se": se,
                    "ci95_lo": g - 1.96 * se,
                    "ci95_hi": g + 1.96 * se,
                    "wilcox_p": wilcoxon(a, b),
                })
    effects = pd.DataFrame(rows)
    effects.to_csv(OUT / "cohort_effects.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    boot_rows = []
    meta_rows = []
    loo_rows = []
    for gene in GENES:
        for subtype in SUBTYPES:
            block = effects[(effects.gene.eq(gene)) & (effects.subtype.eq(subtype))].copy()
            block = block[np.isfinite(block.hedges_g) & (block.se > 0)]
            y = block.hedges_g.to_numpy()
            var = block.se.to_numpy() ** 2
            for method in ("DL_normal", "REML_modified_HK"):
                rec = fit_meta(y, var, method)
                rec.update(gene=gene, subtype=subtype)
                meta_rows.append(rec)
            for cohort in block.cohort:
                keep = block[block.cohort.ne(cohort)]
                rec = fit_meta(keep.hedges_g.to_numpy(), keep.se.to_numpy() ** 2, "DL_normal")
                rec.update(gene=gene, subtype=subtype, omitted_cohort=cohort)
                loo_rows.append(rec)

            # Patient bootstrap: resample within each cohort, recompute g, then DL-pool.
            parts = {
                cohort: part.reset_index(drop=True)
                for cohort, part in tall.groupby("cohort")
            }
            pooled = np.empty(N_BOOT)
            signs = np.empty(N_BOOT)
            kept = 0
            for _ in range(N_BOOT):
                gs, vs = [], []
                ok = True
                for cohort, part in parts.items():
                    draw = part.iloc[rng.integers(0, len(part), len(part))]
                    aa = draw.loc[draw.level1.eq(subtype), gene].to_numpy()
                    bb = draw.loc[~draw.level1.eq(subtype), gene].to_numpy()
                    g, se = hedges_g(aa, bb)
                    if not np.isfinite(g) or not np.isfinite(se) or se <= 0:
                        ok = False
                        break
                    gs.append(g)
                    vs.append(se ** 2)
                if not ok:
                    continue
                pooled[kept] = fit_meta(np.array(gs), np.array(vs), "DL_normal")["hedges_g"]
                signs[kept] = np.sign(np.sum(np.sign(gs)))
                kept += 1
            point = fit_meta(y, var, "DL_normal")
            sample = pooled[:kept]
            boot_rows.append({
                "gene": gene,
                "subtype": subtype,
                "dl_g": point["hedges_g"],
                "boot_n": kept,
                "boot_median": float(np.median(sample)),
                "boot_ci95_lo": float(np.quantile(sample, 0.025)),
                "boot_ci95_hi": float(np.quantile(sample, 0.975)),
                "boot_frac_same_sign_as_point": float(np.mean(np.sign(sample) == np.sign(point["hedges_g"]))),
                "cohort_sign_agreement": point["sign_agreement"],
            })

    meta = pd.DataFrame(meta_rows)
    meta.to_csv(OUT / "random_effects_meta.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows).to_csv(OUT / "leave_one_cohort_out.tsv", sep="\t", index=False)
    boot = pd.DataFrame(boot_rows)
    boot.to_csv(OUT / "patient_bootstrap_meta.tsv", sep="\t", index=False)

    calls = []
    dl = meta[meta.method.eq("DL_normal")].set_index(["subtype", "gene"])
    hk = meta[meta.method.eq("REML_modified_HK")].set_index(["subtype", "gene"])
    for subtype in SUBTYPES:
        for gene in GENES:
            row = dl.loc[(subtype, gene)]
            hkrow = hk.loc[(subtype, gene)]
            brow = boot[(boot.subtype.eq(subtype)) & (boot.gene.eq(gene))].iloc[0]
            same = int(row.sign_agreement.split("/")[0]) == int(row.k)
            dl_excludes = row.ci95_lo > 0 or row.ci95_hi < 0
            boot_excludes = brow.boot_ci95_lo > 0 or brow.boot_ci95_hi < 0
            hk_excludes = hkrow.ci95_lo > 0 or hkrow.ci95_hi < 0
            if same and dl_excludes and boot_excludes:
                level = "direction_stable"
            elif same and dl_excludes:
                level = "analytic_ci_excludes_0_bootstrap_wider"
            elif same:
                level = "same_direction_ci_includes_0"
            else:
                level = "direction_not_shared"
            calls.append({
                "subtype": subtype,
                "gene": gene,
                "dl_g": row.hedges_g,
                "dl_ci95_lo": row.ci95_lo,
                "dl_ci95_hi": row.ci95_hi,
                "dl_p": row.p,
                "I2": row.I2,
                "hk_g": hkrow.hedges_g,
                "hk_ci95_lo": hkrow.ci95_lo,
                "hk_ci95_hi": hkrow.ci95_hi,
                "hk_p": hkrow.p,
                "sign_agreement": row.sign_agreement,
                "boot_ci95_lo": brow.boot_ci95_lo,
                "boot_ci95_hi": brow.boot_ci95_hi,
                "call": level,
                "hk_ci_excludes_0": bool(hk_excludes),
            })
    decision = pd.DataFrame(calls)
    decision.to_csv(OUT / "state_calls.tsv", sep="\t", index=False)
    (OUT / "methods.json").write_text(json.dumps({
        "source": str(SOURCE),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "cohorts": cohorts,
        "n_by_cohort_subtype": {
            f"{cohort}|{subtype}": int(n)
            for (cohort, subtype), n in tall.groupby(["cohort", "level1"]).size().items()
        },
        "balance": "ratio_zeb = ZEB1 - ZEB2 on the cohort-native compressed scale",
        "contrast": "labeled subtype versus the other Immature/TLX/TAL/HOXA samples in the same cohort; NKX and Other excluded from both sides",
        "hedges": "J-corrected Cohen d; SE uses the large-sample Hedges formula from module3/code/02b_analyze_independent.R",
        "meta": ["DL_normal", "REML_modified_HK"],
        "bootstrap": {"seed": SEED, "replicates_requested": N_BOOT, "unit": "patient within cohort, then DL pool"},
        "call_rule": "direction_stable requires 4/4 sign agreement, DL 95% CI excluding 0, and patient-bootstrap 95% CI excluding 0",
        "not_claimed": "causation, TLX1 directly activating ZEB1, or a drug mechanism",
    }, indent=2, default=str), encoding="utf-8")
    show = decision[decision.call.eq("direction_stable")][
        ["subtype", "gene", "dl_g", "dl_ci95_lo", "dl_ci95_hi", "sign_agreement", "boot_ci95_lo", "boot_ci95_hi", "I2", "call"]
    ]
    print(decision.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("\nSTABLE")
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
