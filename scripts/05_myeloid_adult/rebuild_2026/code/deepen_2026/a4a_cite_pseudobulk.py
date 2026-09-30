"""A4-A: HTSA_Ghent paediatric CITE-seq, donor x author-stage pseudobulk.

Frozen:
  B = z(ZEB1) - z(ZEB2)
  D = z(CD1A) - [z(CD34) + z(LYL1)] / 2
z is within this CITE-seq pseudobulk table, not T-ALL.
Unit is donor x author stage, never cells.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import mannwhitneyu

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
adata = ad.read_h5ad(P / "thymus_scrna_tcell_subset.h5ad", backed="r")
obs = adata.obs
ghent = obs["study"].astype(str) == "HTSA_Ghent"
print("HTSA_Ghent", int(ghent.sum()), "donors", obs.loc[ghent, "donor_id"].nunique(), flush=True)
print(obs.loc[ghent, "donor_id"].astype(str).value_counts().to_string(), flush=True)
print("age", obs.loc[ghent, "age_group"].astype(str).value_counts().to_string(), flush=True)

genes = ["ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"]
sym = adata.var["feature_name"].astype(str)
idx = {s: i for i, s in enumerate(sym)}
print("genes", {g: g in idx for g in genes}, flush=True)
# backed CSR cannot take fancy row index; pull 5 gene columns in consecutive chunks
js = [idx[g] for g in genes]
n = adata.n_obs
block = np.zeros((n, len(genes)), dtype=np.float32)
X = adata.X
step = 20000
for start in range(0, n, step):
    stop = min(start + step, n)
    sl = X[start:stop]
    if hasattr(sl, "tocsr"):
        sl = sl.tocsr()
    block[start:stop] = np.asarray(sl[:, js].todense() if hasattr(sl[:, js], "todense") else sl[:, js])
    if start % 100000 == 0:
        print("chunk", start, flush=True)
adata.file.close()
keep = ghent.to_numpy()
mat = {g: block[keep, i] for i, g in enumerate(genes)}
del block

meta = obs.loc[ghent, ["donor_id", "cell_type_level_4_explore", "author_cell_type", "age_group"]].copy()
meta = meta.reset_index(drop=True)
for g in genes:
    meta[g] = mat[g]

# author high-res stage
meta["stage"] = meta["cell_type_level_4_explore"].astype(str)
# drop see_lv4 leftovers if any
meta = meta[meta["stage"] != "see_lv4_explore"]
meta["donor_id"] = meta["donor_id"].astype(str)

# pseudobulk: mean of existing normalized X per donor x stage
# then z across the pseudobulk table
rows = []
for (don, st), sub in meta.groupby(["donor_id", "stage"], observed=True):
    if len(sub) < 10:
        continue
    rec = {"donor_id": don, "stage": st, "n_cells": int(len(sub))}
    for g in genes:
        rec[g] = float(sub[g].mean())
    rows.append(rec)
pb = pd.DataFrame(rows)
for g in genes:
    x = pb[g].to_numpy(float)
    s = x.std(ddof=0)
    pb["z_" + g] = (x - x.mean()) / s if s > 0 else np.nan
pb["balance"] = pb["z_ZEB1"] - pb["z_ZEB2"]
pb["dev"] = pb["z_CD1A"] - (pb["z_CD34"] + pb["z_LYL1"]) / 2

EARLY = {"T_ETP", "T_DN(early)", "T_DN(P)", "T_DN(Q)", "T_DN(Q)-intermediate", "T_DN(Q)-stress_1", "T_DN(Q)-stress_2"}
DP = {"T_DP(P)", "T_DP(Q)", "T_DP(Q)-early", "T_DP(Q)-CD99", "T_DP(Q)-late_vdj", "T_DP(Q)-HSPH1"}
SP = {"T_CD4", "T_CD8", "T_CD8-Prolif", "T_CD8_memory"}


def pole(st):
    if st in EARLY:
        return "early_DN"
    if st in DP:
        return "cortical_DP"
    if st in SP:
        return "SP"
    return "other"


pb["pole"] = pb["stage"].map(pole)
pb.to_csv(P / "a4a_donor_stage_pseudobulk.tsv", sep="\t", index=False)

stage = (
    pb.groupby("stage", observed=True)
    .agg(
        n_donors=("donor_id", "nunique"),
        n_cells=("n_cells", "sum"),
        balance_median=("balance", "median"),
        dev_median=("dev", "median"),
        z_ZEB1=("z_ZEB1", "median"),
        z_ZEB2=("z_ZEB2", "median"),
        pole=("pole", "first"),
    )
    .sort_values("balance_median")
)
stage.to_csv(P / "a4a_stage_balance.tsv", sep="\t")
print(stage.to_string(), flush=True)

a = pb.loc[pb.pole == "early_DN", "balance"].to_numpy(float)
b = pb.loc[pb.pole == "cortical_DP", "balance"].to_numpy(float)
c = pb.loc[pb.pole == "SP", "balance"].to_numpy(float)
p_dp_vs_early = float(mannwhitneyu(b, a, alternative="greater").pvalue) if len(a) and len(b) else np.nan
delta = float(np.median(b) - np.median(a)) if len(a) and len(b) else np.nan

# donor repeat: each donor, median DP balance > median early balance
rep = []
for don, sub in pb.groupby("donor_id"):
    ea = sub.loc[sub.pole == "early_DN", "balance"]
    dp = sub.loc[sub.pole == "cortical_DP", "balance"]
    if len(ea) == 0 or len(dp) == 0:
        continue
    ok = float(dp.median()) > float(ea.median())
    rep.append((don, float(ea.median()), float(dp.median()), ok))
n_ok = sum(1 for *_, ok in rep if ok)
shape_ok = delta > 0
grade = "Pass" if shape_ok and n_ok >= 4 else ("Partial" if shape_ok else "Fail")
lock = pd.DataFrame([
    {"item": "object", "value": "CELLxGENE thymus T-cell subset; study=HTSA_Ghent"},
    {"item": "n_cells_ghent", "value": int(ghent.sum())},
    {"item": "n_donors", "value": int(pb.donor_id.nunique())},
    {"item": "n_pseudobulks", "value": int(len(pb))},
    {"item": "min_cells_per_pseudobulk", "value": 10},
    {"item": "stage_field", "value": "cell_type_level_4_explore"},
    {"item": "early_median_balance", "value": float(np.median(a)) if len(a) else np.nan},
    {"item": "DP_median_balance", "value": float(np.median(b)) if len(b) else np.nan},
    {"item": "SP_median_balance", "value": float(np.median(c)) if len(c) else np.nan},
    {"item": "DP_minus_early", "value": delta},
    {"item": "MWU_DP_gt_early_p", "value": p_dp_vs_early},
    {"item": "donors_DP_gt_early", "value": f"{n_ok}/{len(rep)}"},
    {"item": "gate_A", "value": grade},
    {"item": "SP_monotonic_required", "value": False},
])
lock.to_csv(P / "a4a_gate_lock.tsv", sep="\t", index=False)
print("donor_repeat", rep, flush=True)
print(lock.to_string(index=False), flush=True)
print("grade", grade, flush=True)
