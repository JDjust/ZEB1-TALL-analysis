#!/usr/bin/env python3
"""Patient-independent cell-line audit of the two GDSC methotrexate screens.

Monte Carlo permutation P values are exploratory, with a fixed seed. The two
screens share most models and are never counted as independent biological cohorts.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/F7_mtx_points.tsv"
OUTDIR = ROOT / "data/validation"
SEED = 20260924
N_PERM = 200_000
N_BOOT = 20_000


def rho(a, b):
    return float(spearmanr(a, b).statistic)


def perm_p(x, y, rng):
    xr = rankdata(x); yr = rankdata(y)
    xr = (xr - xr.mean()) / np.linalg.norm(xr - xr.mean())
    yr = (yr - yr.mean()) / np.linalg.norm(yr - yr.mean())
    observed = float(xr @ yr)
    extreme = 0
    for _ in range(N_PERM):
        extreme += abs(float(xr @ yr[rng.permutation(len(yr))])) >= abs(observed) - 1e-12
    return (extreme + 1) / (N_PERM + 1)


def boot_ci(x, y, rng):
    vals = []
    for _ in range(N_BOOT):
        ix = rng.integers(0, len(x), size=len(x))
        if len(np.unique(x[ix])) < 2 or len(np.unique(y[ix])) < 2:
            continue
        vals.append(rho(x[ix], y[ix]))
    return np.quantile(vals, [.025, .975]).tolist(), len(vals)


def main():
    d = pd.read_csv(SOURCE, sep="\t")
    d = d[d.lineage.eq("T-ALL")].dropna(subset=["ZEB1", "value", "model_id"]).copy()
    if d.duplicated(["source", "model_id"]).any():
        raise ValueError("duplicate ModelID within screen")
    screens = {s: v.sort_values("model_id").reset_index(drop=True)
               for s, v in d.groupby("source")}
    if set(screens) != {"GDSC1", "GDSC2"}:
        raise ValueError("expected GDSC1 and GDSC2")
    rng = np.random.default_rng(SEED)
    rows = []
    loo = []
    for name, frame in screens.items():
        x = frame.ZEB1.to_numpy(); y = frame.value.to_numpy()
        ci, nvalid = boot_ci(x, y, rng)
        rows.append(dict(screen=name, n_models=len(frame), rho=rho(x, y),
                         permutation_p=perm_p(x, y, rng), bootstrap_ci_lo=ci[0],
                         bootstrap_ci_hi=ci[1], bootstrap_valid=nvalid,
                         unit="cell line", screen_overlap="not independent"))
        for i, model_id in enumerate(frame.model_id):
            keep = np.arange(len(frame)) != i
            loo.append(dict(screen=name, omitted_model_id=model_id,
                            rho_omit=rho(x[keep], y[keep])))
    wide = d.pivot(index="model_id", columns="source", values="value")
    meta = d.drop_duplicates("model_id").set_index("model_id")[["cell_line", "ZEB1"]]
    paired = meta.join(wide).sort_index()
    paired.to_csv(OUTDIR / "mtx_model_overlap.tsv", sep="\t")
    shared = paired.dropna(subset=["GDSC1", "GDSC2"])
    summary = {"n_shared": len(shared),
               "rho_screen_response_shared": rho(shared.GDSC1, shared.GDSC2),
               "rho_zeb1_gdsc1_shared": rho(shared.ZEB1, shared.GDSC1),
               "rho_zeb1_gdsc2_shared": rho(shared.ZEB1, shared.GDSC2),
               "nonshared_gdsc1": paired[paired.GDSC1.notna() & paired.GDSC2.isna()].index.tolist(),
               "nonshared_gdsc2": paired[paired.GDSC2.notna() & paired.GDSC1.isna()].index.tolist(),
               "n_permutations": N_PERM, "n_bootstrap": N_BOOT, "seed": SEED,
               "caution": "shared cell lines; cross-screen consistency, not independent validation"}
    pd.DataFrame(rows).to_csv(OUTDIR / "mtx_screen_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(loo).to_csv(OUTDIR / "mtx_leave_one_model_out.tsv", sep="\t", index=False)
    (OUTDIR / "mtx_overlap_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(pd.DataFrame(rows).to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
