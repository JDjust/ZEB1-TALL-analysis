#!/usr/bin/env python3
"""Who accounts for daunorubicin LC50: ZEB1, ZEB2, or their balance?

Primary cohort is the St. Jude pharmacotype T-ALL table already used for
Figure 7. Spearman is the rank association and is invariant to a monotone
transform of LC50. Linear models use log2(LC50) so a few extreme LC50 values
do not dominate ordinary least squares. Huber regression is the robust check
on the same log2 scale. Age and ETP are the requested covariates. Molecular
subtype is a sensitivity, because a ratio effect that is only subtype would
not be an independent pharmacologic association.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
from statsmodels.regression.linear_model import OLS
from statsmodels.robust.robust_linear_model import RLM
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools import add_constant

SOURCE = Path(__file__).resolve().parents[3] / "TALL_resistance_atlas" / "pharmacotype_tall_drugs.tsv"
OUT = Path(__file__).resolve().parents[1] / "data" / "validation" / "zeb_family_state"
SEED = 20260925
N_BOOT = 2000
DRUG = "Daunorubicin (µM)"


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


def spearman_boot(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> dict:
    rho, p = spearmanr(x, y)
    n = len(x)
    draws = np.empty(N_BOOT)
    for i in range(N_BOOT):
        ix = rng.integers(0, n, n)
        draws[i] = spearmanr(x[ix], y[ix]).statistic
    return {
        "n": n,
        "rho": float(rho),
        "p": float(p),
        "boot_ci95_lo": float(np.quantile(draws, 0.025)),
        "boot_ci95_hi": float(np.quantile(draws, 0.975)),
    }


def residual_spearman(frame: pd.DataFrame, outcome: str, predictor: str, covariates: str) -> dict:
    y = smf.ols(f"{outcome} ~ {covariates}", frame).fit().resid
    x = smf.ols(f"{predictor} ~ {covariates}", frame).fit().resid
    rho, p = spearmanr(x, y)
    return {"rho": float(rho), "p": float(p), "n": int(len(frame))}


def fit_ols(frame: pd.DataFrame, formula: str):
    return smf.ols(formula, frame).fit(cov_type="HC3")


def huber_beta(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    model = RLM(y, add_constant(x), M=None)
    # default HuberT
    return np.asarray(model.fit().params, float)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ph = pd.read_csv(SOURCE, sep="\t")
    ph = ph[ph["disease"].eq("T-ALL")].copy()
    ph["lc50"] = pd.to_numeric(ph[DRUG], errors="coerce")
    ph["age"] = pd.to_numeric(ph["Age at diagnosis (years)"], errors="coerce")
    ph["ZEB1_log"] = np.log2(ph["ZEB1"] + 1)
    ph["ZEB2_log"] = np.log2(ph["ZEB2"] + 1)
    ph["ratio_zeb"] = np.log2((ph["ZEB1"] + 1e-6) / (ph["ZEB2"] + 1e-6))
    ph["etp_flag"] = ph["etp"].map({"ETP": 1, "notETP": 0})
    ph["log2_lc50"] = np.log2(ph["lc50"])
    use = ph.dropna(subset=["lc50", "age", "ZEB1", "ZEB2", "etp_flag", "log2_lc50"]).copy()
    use = use[np.isfinite(use["log2_lc50"]) & (use["lc50"] > 0)].copy()
    for col in ("ZEB1_log", "ZEB2_log", "ratio_zeb", "age"):
        use[col + "_z"] = zscore(use[col])

    rng = np.random.default_rng(SEED)
    spearman_rows = []
    for name, col in (("ZEB1_log", "ZEB1_log"), ("ZEB2_log", "ZEB2_log"), ("ratio_zeb", "ratio_zeb")):
        rec = spearman_boot(use[col].to_numpy(), use["lc50"].to_numpy(), rng)
        rec["predictor"] = name
        rec["adjustment"] = "none"
        # leave-one-out rho
        x = use[col].to_numpy()
        y = use["lc50"].to_numpy()
        loos = []
        for i in range(len(use)):
            m = np.ones(len(use), dtype=bool)
            m[i] = False
            loos.append(spearmanr(x[m], y[m]).statistic)
        rec["loo_min"] = float(np.min(loos))
        rec["loo_max"] = float(np.max(loos))
        rec["loo_sign_flip"] = bool(np.min(loos) * np.max(loos) < 0)
        spearman_rows.append(rec)
        adj = residual_spearman(use, "log2_lc50", col, "age + etp_flag")
        spearman_rows.append({
            "predictor": name,
            "adjustment": "rank_residual_age_ETP",
            **adj,
            "boot_ci95_lo": np.nan,
            "boot_ci95_hi": np.nan,
            "loo_min": np.nan,
            "loo_max": np.nan,
            "loo_sign_flip": np.nan,
        })
    # subtype sensitivity on the same complete cases that have a subtype label
    sub = use[use["subtype"].fillna("").ne("")].copy()
    if sub["subtype"].nunique() > 1:
        for name in ("ZEB1_log", "ZEB2_log", "ratio_zeb"):
            adj = residual_spearman(sub, "log2_lc50", name, "age + etp_flag + C(subtype)")
            spearman_rows.append({
                "predictor": name,
                "adjustment": "rank_residual_age_ETP_subtype",
                **adj,
            })
    pd.DataFrame(spearman_rows).to_csv(OUT / "daunorubicin_spearman.tsv", sep="\t", index=False)

    formulas = {
        "ZEB1": "log2_lc50 ~ ZEB1_log_z",
        "ZEB2": "log2_lc50 ~ ZEB2_log_z",
        "ratio": "log2_lc50 ~ ratio_zeb_z",
        "ZEB1_ZEB2": "log2_lc50 ~ ZEB1_log_z + ZEB2_log_z",
        "ZEB1_age_ETP": "log2_lc50 ~ ZEB1_log_z + age_z + etp_flag",
        "ZEB2_age_ETP": "log2_lc50 ~ ZEB2_log_z + age_z + etp_flag",
        "ratio_age_ETP": "log2_lc50 ~ ratio_zeb_z + age_z + etp_flag",
        "ZEB1_ZEB2_age_ETP": "log2_lc50 ~ ZEB1_log_z + ZEB2_log_z + age_z + etp_flag",
    }
    if use["subtype"].nunique() > 1:
        formulas["ratio_age_ETP_subtype"] = "log2_lc50 ~ ratio_zeb_z + age_z + etp_flag + C(subtype)"
        formulas["both_age_ETP_subtype"] = "log2_lc50 ~ ZEB1_log_z + ZEB2_log_z + age_z + etp_flag + C(subtype)"

    coef_rows = []
    fit_rows = []
    y = use["log2_lc50"].to_numpy()
    n = len(use)
    index = np.arange(n)
    for name, formula in formulas.items():
        fit = fit_ols(use, formula)
        plain = smf.ols(formula, use).fit()
        fit_rows.append({
            "model": name,
            "formula": formula,
            "n": int(fit.nobs),
            "r2": float(plain.rsquared),
            "adj_r2": float(plain.rsquared_adj),
            "aic": float(plain.aic),
            "llf": float(plain.llf),
        })
        terms = [t for t in fit.params.index if t != "Intercept" and not str(t).startswith("C(")]
        # design matrix column positions for Huber / bootstrap / LOO
        exog = plain.model.exog
        names = list(plain.model.exog_names)
        huber = huber_beta(y, exog[:, 1:] if exog.shape[1] > 1 else exog)
        # huber_beta adds constant, so params align with exog columns if exog already has constant
        # RLM(add_constant(exog without intercept)) duplicates if exog has intercept.
        # Use exog as-is when it already contains the intercept.
        huber = np.asarray(RLM(y, exog).fit().params, float)
        boot = np.empty((N_BOOT, exog.shape[1]))
        for i in range(N_BOOT):
            ix = rng.integers(0, n, n)
            boot[i] = np.asarray(RLM(y[ix], exog[ix]).fit().params, float)
        loo = np.empty((n, exog.shape[1]))
        for i, drop in enumerate(index):
            keep = np.ones(n, dtype=bool)
            keep[drop] = False
            loo[i] = OLS(y[keep], exog[keep]).fit().params
        for term in names:
            j = names.index(term)
            if term == "Intercept":
                continue
            coef_rows.append({
                "model": name,
                "term": term,
                "beta_z": float(fit.params[term]),
                "se_hc3": float(fit.bse[term]),
                "p_hc3": float(fit.pvalues[term]),
                "ci95_lo": float(fit.conf_int().loc[term, 0]),
                "ci95_hi": float(fit.conf_int().loc[term, 1]),
                "huber_beta": float(huber[j]),
                "huber_boot_lo": float(np.quantile(boot[:, j], 0.025)),
                "huber_boot_hi": float(np.quantile(boot[:, j], 0.975)),
                "loo_min": float(loo[:, j].min()),
                "loo_max": float(loo[:, j].max()),
                "loo_sign_flip": bool(loo[:, j].min() * loo[:, j].max() < 0),
            })
    pd.DataFrame(fit_rows).to_csv(OUT / "daunorubicin_model_fit.tsv", sep="\t", index=False)
    pd.DataFrame(coef_rows).to_csv(OUT / "daunorubicin_coefficients.tsv", sep="\t", index=False)

    both = add_constant(use[["ZEB1_log_z", "ZEB2_log_z"]])
    vif = {
        col: float(variance_inflation_factor(both.to_numpy(), i))
        for i, col in enumerate(both.columns) if col != "const"
    }
    corr = float(use[["ZEB1_log", "ZEB2_log"]].corr(method="spearman").iloc[0, 1])
    (OUT / "daunorubicin_methods.json").write_text(json.dumps({
        "source": str(SOURCE),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "drug": DRUG,
        "n_complete_age_etp": int(len(use)),
        "n_source_tall": int(ph["disease"].eq("T-ALL").sum()) if "disease" in ph else int(len(ph)),
        "outcome_spearman": "raw LC50; rank association",
        "outcome_regression": "log2(LC50), predictors z-scored on the analysis set",
        "ratio": "log2((ZEB1+1e-6)/(ZEB2+1e-6)), the Module 7 definition",
        "covariates": ["age", "ETP vs notETP"],
        "sensitivity": "molecular subtype as an additional covariate",
        "robust": "Huber RLM on the same z-scored design; percentile bootstrap of Huber coefficients",
        "loo": "ordinary least-squares coefficient after dropping one patient",
        "spearman_ZEB1_ZEB2": corr,
        "vif_ZEB1_ZEB2": vif,
        "seed": SEED,
        "bootstrap": N_BOOT,
        "interpretation_limit": "One pharmacotype cohort. Association after age/ETP is not a mechanism and is not independent replication.",
    }, indent=2), encoding="utf-8")
    print(pd.DataFrame(spearman_rows).to_string(index=False))
    print(pd.DataFrame(fit_rows).to_string(index=False))
    print(pd.DataFrame(coef_rows).to_string(index=False))
    print("VIF", vif, "spearman ZEB1-ZEB2", corr)


if __name__ == "__main__":
    main()
