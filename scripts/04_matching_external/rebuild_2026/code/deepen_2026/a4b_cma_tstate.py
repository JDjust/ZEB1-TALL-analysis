"""A4-B spatial half: author T-state CMA localization. No raw Visium ZEB1-ZEB2."""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
adata = ad.read_h5ad(P / "merged_thymus_visium_pediatric.h5ad", backed="r")
obs = adata.obs.copy()
adata.file.close()


def rename_t(c: str) -> str:
    if not c.startswith("T_"):
        return c
    if "CD8" in c and "entry" in c:
        return "T_CD8ab(entry)"
    if "CD8" in c and "(II)" in c:
        return "T_CD8ab(II)"
    if "CD8" in c and "(I)" in c:
        return "T_CD8ab(I)"
    if "entry" in c and "CD8" not in c and "DN" not in c and "DP" not in c:
        return "T_gdT(entry)"
    known = {
        "T_CD4", "T_CD8", "T_CD8_memory", "T_DN(P)", "T_DN(Q)", "T_DN(early)",
        "T_DP(P)", "T_DP(Q)", "T_DP(Q)-CD99", "T_DP(Q)-early", "T_DP(Q)-late_vdj",
        "T_ETP", "T_NK", "T_Treg(agonist)", "T_Treg-diff_1", "T_Treg-diff_2",
        "T_Treg_CD8", "T_Treg_mature", "T_Treg_recirc",
    }
    if c in known:
        return c
    if "Treg" in c or "DN" in c or "DP" in c or "ETP" in c or "NK" in c:
        return c
    return "T_gdT"


obs = obs.rename(columns={c: rename_t(c) for c in obs.columns})
tcols = [c for c in obs.columns if c.startswith("T_")]
print("T_states", tcols, flush=True)

bin_cma = (
    obs.groupby("manual_bin_cma_v2", observed=True)["cma_v2"]
    .agg(n="size", median="median", mean="mean")
    .sort_values("median")
)
bin_cma.to_csv(P / "a4b_cma_by_author_bin.tsv", sep="\t")
print(bin_cma, flush=True)

cma = obs["cma_v2"].to_numpy(float)
donors = obs["donor_id"].astype(str)


def wmedian(values, weights):
    values = np.asarray(values, float)
    weights = np.asarray(weights, float)
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[ok], weights[ok]
    if weights.size == 0:
        return np.nan
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights)
    return float(values[np.searchsorted(cdf, 0.5 * cdf[-1])])


rows = []
for st in tcols:
    w = pd.to_numeric(obs[st], errors="coerce").to_numpy(float)
    w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
    if w.sum() <= 0:
        continue
    rows.append({
        "state": st,
        "n_spots_pos": int((w > 0).sum()),
        "abundance_sum": float(w.sum()),
        "cma_weighted_mean": float(np.average(cma, weights=w)),
        "cma_weighted_median": wmedian(cma, w),
    })
loc = pd.DataFrame(rows).sort_values("cma_weighted_median")
loc.to_csv(P / "a4b_tstate_cma_localization.tsv", sep="\t", index=False)
print(loc.to_string(index=False), flush=True)

drows = []
for don in sorted(donors.unique()):
    m = donors.to_numpy() == don
    for st in tcols:
        w = pd.to_numeric(obs.loc[m, st], errors="coerce").to_numpy(float)
        med = wmedian(cma[m], w)
        if np.isnan(med):
            continue
        drows.append({
            "donor_id": don,
            "state": st,
            "cma_weighted_median": med,
            "abundance_sum": float(np.nansum(w)),
        })
dloc = pd.DataFrame(drows)
dloc.to_csv(P / "a4b_tstate_cma_by_donor.tsv", sep="\t", index=False)

early = ["T_ETP", "T_DN(early)"]
cortical = ["T_DP(P)", "T_DP(Q)", "T_DP(Q)-early", "T_DP(Q)-CD99"]
late = ["T_CD4", "T_CD8", "T_CD8ab(I)", "T_CD8ab(II)"]
idx = loc.set_index("state")["cma_weighted_median"]


def med_of(names):
    hit = [n for n in names if n in idx.index]
    return float(np.median([idx[n] for n in hit])) if hit else np.nan


early_cma = med_of(early)
cort_cma = med_of(cortical)
late_cma = med_of(late)
cortex_bins = [b for b in bin_cma.index.astype(str) if ("Cortical" in b or "Capsular" in b) and "CMJ" not in b]
medulla_bins = [b for b in bin_cma.index.astype(str) if "Medullar" in b and "CMJ" not in b]
cortex_sign = float(np.sign(
    bin_cma.loc[cortex_bins, "median"].median() - bin_cma.loc[medulla_bins, "median"].median()
))
dp_more_cortical = ((cort_cma - early_cma) * cortex_sign) > 0
sp_return = ((late_cma - cort_cma) * cortex_sign) < 0
rep = []
for don, sub in dloc.groupby("donor_id"):
    s = sub.set_index("state")["cma_weighted_median"]
    if "T_ETP" in s.index and "T_DP(Q)" in s.index:
        delta = (s["T_DP(Q)"] - s["T_ETP"]) * cortex_sign
        rep.append((don, float(delta), bool(delta > 0)))
n_rep = sum(1 for _, _, ok in rep if ok)
spatial = "Pass" if dp_more_cortical and n_rep >= 4 else ("Partial" if dp_more_cortical else "Fail")
lock = pd.DataFrame([
    {"item": "n_donors", "value": int(obs["donor_id"].nunique())},
    {"item": "n_sections", "value": int(obs["library_id"].nunique())},
    {"item": "n_spots", "value": int(len(obs))},
    {"item": "cma_column", "value": "cma_v2"},
    {"item": "deconvolution", "value": "author_cell2location_T_states"},
    {"item": "cortex_minus_medulla_sign", "value": cortex_sign},
    {"item": "early_cma_median", "value": early_cma},
    {"item": "cortical_DP_cma_median", "value": cort_cma},
    {"item": "late_SP_cma_median", "value": late_cma},
    {"item": "DP_more_cortical_than_early", "value": dp_more_cortical},
    {"item": "SP_returns_vs_DP", "value": sp_return},
    {"item": "donors_DPQ_vs_ETP_cortical", "value": f"{n_rep}/{len(rep)}"},
    {"item": "gate_B_spatial_half", "value": spatial},
    {"item": "gate_B_full", "value": "blocked_until_A4A_ZEB_balance_on_same_states"},
    {"item": "raw_ZEB_on_spots", "value": "not_used"},
])
lock.to_csv(P / "a4b_spatial_half_lock.tsv", sep="\t", index=False)
print(lock.to_string(index=False), flush=True)
print("donor_DPQ_vs_ETP", rep, flush=True)
