"""A4-A using the same migration poles already locked for A4-B spatial half."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
pb = pd.read_csv(P / "a4a_donor_stage_pseudobulk.tsv", sep="\t")

# identical to a4b_cma_tstate.py
EARLY = {"T_ETP", "T_DN(early)"}
DP = {"T_DP(P)", "T_DP(Q)", "T_DP(Q)-early", "T_DP(Q)-CD99"}
SP = {"T_CD4", "T_CD8", "T_CD8ab(I)", "T_CD8ab(II)"}
# CITE-seq uses unicode alpha names
SP |= {s for s in pb.stage.unique() if s.startswith("T_CD8") and s != "T_CD8_memory"}
SP |= {"T_CD4", "T_CD8"}


def pole(st):
    if st in EARLY:
        return "early"
    if st in DP:
        return "cortical_DP"
    if st in {"T_CD4", "T_CD8"}:
        return "SP"
    return "other"


pb["pole_locked"] = pb["stage"].map(pole)
a = pb.loc[pb.pole_locked == "early", "balance"].to_numpy(float)
b = pb.loc[pb.pole_locked == "cortical_DP", "balance"].to_numpy(float)
c = pb.loc[pb.pole_locked == "SP", "balance"].to_numpy(float)
p = float(mannwhitneyu(b, a, alternative="greater").pvalue) if len(a) and len(b) else np.nan
delta = float(np.median(b) - np.median(a)) if len(a) and len(b) else np.nan
rep = []
for don, sub in pb.groupby("donor_id"):
    ea = sub.loc[sub.pole_locked == "early", "balance"]
    dp = sub.loc[sub.pole_locked == "cortical_DP", "balance"]
    if len(ea) == 0 or len(dp) == 0:
        continue
    ok = float(dp.median()) > float(ea.median())
    rep.append((don, int(len(ea)), float(ea.median()), float(dp.median()), ok))
n_ok = sum(1 for *_, ok in rep if ok)
grade = "Pass" if delta > 0 and n_ok >= 4 else ("Partial" if delta > 0 else "Fail")
# join CMA from A4-B
cma = pd.read_csv(P / "a4b_tstate_cma_localization.tsv", sep="\t")
st = pb.groupby("stage", observed=True).agg(
    n_donors=("donor_id", "nunique"),
    n_cells=("n_cells", "sum"),
    balance_median=("balance", "median"),
    pole=("pole_locked", "first"),
).reset_index()
# name alignment for CMA join
alias = {
    "T_αβT(entry)": "T_gdT(entry)",
}
# try join on exact then common prefixes
cma_map = dict(zip(cma["state"], cma["cma_weighted_median"]))
st["cma_weighted_median"] = st["stage"].map(cma_map)
st.to_csv(P / "a4ab_stage_balance_and_cma.tsv", sep="\t", index=False)
lock = pd.DataFrame([
    {"item": "poles", "value": "same as A4-B: ETP+DN(early) vs DP(P/Q/early/CD99) vs CD4/CD8"},
    {"item": "early_n_pseudobulks", "value": int(len(a))},
    {"item": "DP_n_pseudobulks", "value": int(len(b))},
    {"item": "early_median_balance", "value": float(np.median(a)) if len(a) else np.nan},
    {"item": "DP_median_balance", "value": float(np.median(b)) if len(b) else np.nan},
    {"item": "SP_median_balance", "value": float(np.median(c)) if len(c) else np.nan},
    {"item": "DP_minus_early", "value": delta},
    {"item": "MWU_DP_gt_early_p", "value": p},
    {"item": "donors_DP_gt_early", "value": f"{n_ok}/{len(rep)}"},
    {"item": "gate_A", "value": grade},
    {"item": "all_DN_as_early", "value": "not_used; DN(P)/DN(Q) are cortical on CMA"},
])
lock.to_csv(P / "a4a_gate_lock.tsv", sep="\t", index=False)
print(st.sort_values("balance_median").to_string(index=False), flush=True)
print("donor_repeat", rep, flush=True)
print(lock.to_string(index=False), flush=True)
print("grade", grade, flush=True)
