"""Four-cohort TLX sensitivity: DL-normal, REML modified HK, prediction interval."""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm, t

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / "modules/module3/tables/M3_r2_subtype_vs_rest.tsv"
OUT = ROOT / "data/validation"


def fit(y, variance, method):
    k = len(y)
    if method == "DL_normal":
        w0 = 1 / variance
        mu0 = np.sum(w0 * y) / w0.sum()
        q0 = np.sum(w0 * (y - mu0) ** 2)
        tau2 = max(0., (q0 - k + 1) / (w0.sum() - np.sum(w0 ** 2) / w0.sum()))
    else:
        def objective(tau):
            w = 1 / (variance + tau)
            mu = np.sum(w * y) / w.sum()
            return .5 * (np.log(variance + tau).sum() + np.log(w.sum()) + np.sum(w * (y - mu) ** 2))
        upper = max(1., float(np.var(y) * 10))
        result = minimize_scalar(objective, bounds=(0, upper), method="bounded", options={"xatol": 1e-12})
        assert result.success
        tau2 = 0. if objective(0) <= result.fun else float(result.x)
    w = 1 / (variance + tau2)
    mu = float(np.sum(w * y) / w.sum())
    q = float(np.sum(w * (y - mu) ** 2) / (k - 1))
    se = float(np.sqrt((max(1., q) if method == "REML_modified_HK" else 1.) / w.sum()))
    cutoff = float(t.ppf(.975, k - 1)) if method == "REML_modified_HK" else float(norm.ppf(.975))
    p = float(2 * t.sf(abs(mu / se), k - 1)) if method == "REML_modified_HK" else float(2 * norm.sf(abs(mu / se)))
    pred = float(t.ppf(.975, k - 2) * np.sqrt(tau2 + se ** 2))
    return dict(method=method, k=k, estimate=mu, se=se, tau2=tau2,
                ci95_lo=mu-cutoff*se, ci95_hi=mu+cutoff*se, p=p,
                prediction95_lo=mu-pred, prediction95_hi=mu+pred,
                hk_scale=q, hk_scale_used=max(1., q) if method == "REML_modified_HK" else 1.)


def main():
    data = pd.read_csv(SOURCE, sep="\t")
    data = data[data.gene.eq("ZEB1") & data.subtype.eq("TLX")].copy()
    assert len(data) == 4 and data.cohort.is_unique and (data.se > 0).all()
    data.to_csv(OUT / "tlx_meta_input.tsv", sep="\t", index=False)
    records = [fit(data.hedges_g.to_numpy(), data.se.to_numpy() ** 2, method)
               for method in ["DL_normal", "REML_modified_HK"]]
    original = pd.read_csv(ROOT / "data/SourceData_Fig3E_meta.tsv", sep="\t")
    reference = original[original.contrast.eq("TLX vs rest")].iloc[0]
    assert abs(records[0]["estimate"] - reference.hedges_g) < 1e-6
    assert abs(records[0]["ci95_lo"] - reference.ci95_lo) < 1e-4
    pd.DataFrame(records).to_csv(OUT / "tlx_meta_sensitivity.tsv", sep="\t", index=False)
    omissions = []
    for cohort in data.cohort:
        remaining = data[data.cohort.ne(cohort)]
        result = fit(remaining.hedges_g.to_numpy(), remaining.se.to_numpy() ** 2, "REML_modified_HK")
        result["omitted_cohort"] = cohort
        omissions.append(result)
    pd.DataFrame(omissions).to_csv(OUT / "tlx_meta_leave_one_out.tsv", sep="\t", index=False)
    (OUT / "tlx_meta_methods.json").write_text(json.dumps({
        "input_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "reml": "Restricted normal likelihood, bounded scalar tau-squared optimization",
        "modified_HK": "Variance inflation max(1,Q_RE/(k-1)); t interval df=k-1",
        "prediction_interval": "Exploratory t(k-2)*sqrt(tau2+SE_mean^2); imprecise for k=4",
        "interpretation": "Same direction across four cohorts; label harmonization and small-k uncertainty retained",
        "retrospective_sensitivity": True}, indent=2), encoding="utf-8")
    print(pd.DataFrame(records).to_string(index=False))


if __name__ == "__main__":
    main()
