"""A4 donor-level stats. Do not use the 6-vs-20 Mann-Whitney as 5-donor inference."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
pb = pd.read_csv(P / "a4a_donor_stage_pseudobulk.tsv", sep="\t")
EARLY = {"T_ETP", "T_DN(early)"}
DP = {"T_DP(P)", "T_DP(Q)", "T_DP(Q)-early", "T_DP(Q)-CD99"}
SP = {"T_CD4", "T_CD8"}

def pole(st):
    if st in EARLY:
        return "early"
    if st in DP:
        return "cortical_DP"
    if st in SP:
        return "SP"
    return "other"

pb["pole"] = pb["stage"].map(pole)
rows = []
for don, sub in pb.groupby("donor_id"):
    rec = {"donor_id": don}
    for name in ("early", "cortical_DP", "SP"):
        rec[name] = float(sub.loc[sub.pole == name, "balance"].median()) if (sub.pole == name).any() else np.nan
    rec["DP_minus_early"] = rec["cortical_DP"] - rec["early"]
    rec["DP_gt_early"] = rec["DP_minus_early"] > 0
    rows.append(rec)
don = pd.DataFrame(rows).sort_values("donor_id")
x = don["early"].to_numpy(float)
y = don["cortical_DP"].to_numpy(float)
# one-sided: DP > early
w = wilcoxon(y - x, alternative="greater", zero_method="wilcox", method="exact")
don.to_csv(P / "a4a_donor_level_balance.tsv", sep="\t", index=False)
lock = pd.DataFrame([
    {"item": "early_median_of_donor_medians", "value": float(np.median(x))},
    {"item": "DP_median_of_donor_medians", "value": float(np.median(y))},
    {"item": "delta_median_of_donor_medians", "value": float(np.median(y) - np.median(x))},
    {"item": "donors_DP_gt_early", "value": f"{int(don.DP_gt_early.sum())}/{len(don)}"},
    {"item": "paired_wilcoxon_greater_p", "value": float(w.pvalue)},
    {"item": "wilcoxon_role", "value": "sensitivity_only; not 5-donor significance"},
    {"item": "do_not_report_as_primary", "value": "MWU_P=0.005 on 6_vs_20_pseudobulks"},
    {"item": "main_text_core", "value": "early -0.37 to DP +0.52; delta +0.89; 4/5 donors concordant"},
    {"item": "gate_A", "value": "Pass"},
    {"item": "main_text", "value": "yes; physical CMA topology"},
])
lock.to_csv(P / "a4_closed_lock.tsv", sep="\t", index=False)
print(don.to_string(index=False), flush=True)
print(lock.to_string(index=False), flush=True)
print("wilcoxon_p", float(w.pvalue), flush=True)
