"""A4-A Gate A from streamed HTSA_Ghent five-gene table."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a4_yayon")
meta = pd.read_csv(P / "a4a_ghent_five_genes.tsv", sep="\t")
meta["stage"] = meta["cell_type_level_4_explore"].astype(str)
meta = meta[meta["stage"] != "see_lv4_explore"]
genes = ["ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"]

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
p = float(mannwhitneyu(b, a, alternative="greater").pvalue) if len(a) and len(b) else np.nan
delta = float(np.median(b) - np.median(a)) if len(a) and len(b) else np.nan
rep = []
for don, sub in pb.groupby("donor_id"):
    ea = sub.loc[sub.pole == "early_DN", "balance"]
    dp = sub.loc[sub.pole == "cortical_DP", "balance"]
    if len(ea) == 0 or len(dp) == 0:
        continue
    ok = float(dp.median()) > float(ea.median())
    rep.append((don, float(ea.median()), float(dp.median()), ok))
n_ok = sum(1 for *_, ok in rep if ok)
grade = "Pass" if delta > 0 and n_ok >= 4 else ("Partial" if delta > 0 else "Fail")
lock = pd.DataFrame([
    {"item": "object", "value": "CELLxGENE T-cell subset; study=HTSA_Ghent"},
    {"item": "n_cells_ghent", "value": int(len(meta))},
    {"item": "n_donors", "value": int(pb.donor_id.nunique())},
    {"item": "n_pseudobulks", "value": int(len(pb))},
    {"item": "min_cells_per_pseudobulk", "value": 10},
    {"item": "stage_field", "value": "cell_type_level_4_explore"},
    {"item": "early_median_balance", "value": float(np.median(a)) if len(a) else np.nan},
    {"item": "DP_median_balance", "value": float(np.median(b)) if len(b) else np.nan},
    {"item": "SP_median_balance", "value": float(np.median(c)) if len(c) else np.nan},
    {"item": "DP_minus_early", "value": delta},
    {"item": "MWU_DP_gt_early_p", "value": p},
    {"item": "donors_DP_gt_early", "value": f"{n_ok}/{len(rep)}"},
    {"item": "gate_A", "value": grade},
    {"item": "SP_monotonic_required", "value": False},
])
lock.to_csv(P / "a4a_gate_lock.tsv", sep="\t", index=False)
print("donor_repeat", rep, flush=True)
print(lock.to_string(index=False), flush=True)
print("grade", grade, flush=True)
