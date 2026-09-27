"""A6 formal: inheritance vs reconfiguration. Donor-level stats only."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026")
out = P / "a6_yayon"
frozen = pd.read_csv(P / "gate2_program/frozen_ZEB_side50.tsv", sep="\t")
pb = pd.read_csv(out / "a6_donor_stage_pseudobulk.tsv", sep="\t")
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
z1 = [s for s in frozen.loc[frozen.side == "ZEB1-side50", "symbol"].astype(str) if s in pb.columns]
z2 = [s for s in frozen.loc[frozen.side == "ZEB2-side50", "symbol"].astype(str) if s in pb.columns]
rec = z1 + z2
print("recovered ZEB1-side", len(z1), "/50", "ZEB2-side", len(z2), "/50", flush=True)

# gene-wise z across donor x stage table
for g in rec:
    x = pb[g].to_numpy(float)
    s = x.std(ddof=0)
    pb["z_" + g] = (x - x.mean()) / s if s > 0 else np.nan
pb["ZEB1_side50"] = pb[["z_" + g for g in z1]].mean(axis=1)
pb["ZEB2_side50"] = pb[["z_" + g for g in z2]].mean(axis=1)
pb.to_csv(out / "a6_donor_stage_scores.tsv", sep="\t", index=False)

# A6-1 donor-level
drows = []
for don, sub in pb.groupby("donor_id"):
    recd = {"donor_id": don}
    for name, col in [("early", "early"), ("DP", "cortical_DP"), ("SP", "SP")]:
        sl = sub.loc[sub.pole == col]
        recd[f"{name}_ZEB1_side50"] = float(sl["ZEB1_side50"].median()) if len(sl) else np.nan
        recd[f"{name}_ZEB2_side50"] = float(sl["ZEB2_side50"].median()) if len(sl) else np.nan
    recd["d_ZEB1_DP_minus_early"] = recd["DP_ZEB1_side50"] - recd["early_ZEB1_side50"]
    recd["d_ZEB2_DP_minus_early"] = recd["DP_ZEB2_side50"] - recd["early_ZEB2_side50"]
    recd["ZEB1_DP_gt_early"] = recd["d_ZEB1_DP_minus_early"] > 0
    recd["ZEB2_early_gt_DP"] = recd["d_ZEB2_DP_minus_early"] < 0
    drows.append(recd)
don = pd.DataFrame(drows).sort_values("donor_id")
don.to_csv(out / "a6_1_donor_program.tsv", sep="\t", index=False)

def ntrue(s):
    return f"{int(s.sum())}/{int(s.notna().sum())}"

a61 = pd.DataFrame([
    {"item": "ZEB1_side50_recovered", "value": f"{len(z1)}/50"},
    {"item": "ZEB2_side50_recovered", "value": f"{len(z2)}/50"},
    {"item": "early_ZEB1_median", "value": float(don.early_ZEB1_side50.median())},
    {"item": "DP_ZEB1_median", "value": float(don.DP_ZEB1_side50.median())},
    {"item": "donors_ZEB1_DP_gt_early", "value": ntrue(don.ZEB1_DP_gt_early)},
    {"item": "early_ZEB2_median", "value": float(don.early_ZEB2_side50.median())},
    {"item": "DP_ZEB2_median", "value": float(don.DP_ZEB2_side50.median())},
    {"item": "donors_ZEB2_early_gt_DP", "value": ntrue(don.ZEB2_early_gt_DP)},
    {"item": "disease_implied_ZEB1", "value": "DP_higher"},
    {"item": "disease_implied_ZEB2", "value": "early_higher"},
])
# classification
z1_track = int(don.ZEB1_DP_gt_early.sum())
z2_track = int(don.ZEB2_early_gt_DP.sum())
if z1_track >= 4 and z2_track >= 4:
    klass = "A_inheritance"
elif z1_track <= 1 and z2_track <= 1:
    klass = "B_reconfiguration"
else:
    klass = "C_coexist_or_decoupled"
a61 = pd.concat([a61, pd.DataFrame([{"item": "A6_1_class", "value": klass}])], ignore_index=True)
a61.to_csv(out / "a6_1_lock.tsv", sep="\t", index=False)
print(don.to_string(index=False), flush=True)
print(a61.to_string(index=False), flush=True)

# A6-2 gene-level Δ_normal vs β_R
beta = frozen.set_index("symbol")["beta_R"].to_dict()
side = frozen.set_index("symbol")["side"].to_dict()
grows = []
for g in rec:
    dvals = []
    npos = 0
    for _, sub in pb.groupby("donor_id"):
        ea = sub.loc[sub.pole == "early", g]
        dp = sub.loc[sub.pole == "cortical_DP", g]
        if len(ea) == 0 or len(dp) == 0:
            continue
        d = float(dp.mean()) - float(ea.mean())
        dvals.append(d)
        npos += int(d > 0)
    if not dvals:
        continue
    delta = float(np.mean(dvals))
    grows.append({
        "symbol": g,
        "side": side[g],
        "beta_R": float(beta[g]),
        "delta_normal": delta,
        "n_donors": len(dvals),
        "n_donors_DP_gt_early": npos,
    })
ge = pd.DataFrame(grows)
thr = float(ge["delta_normal"].abs().median())
ge["normal_weak"] = ge["delta_normal"].abs() < thr
same = np.sign(ge["delta_normal"]) == np.sign(ge["beta_R"])
ge["class"] = np.where(ge["normal_weak"], "neutral", np.where(same, "concordant", "discordant"))
rho, rp = spearmanr(ge["delta_normal"], ge["beta_R"])
ge.to_csv(out / "a6_2_gene_effects.tsv", sep="\t", index=False)
a62 = pd.DataFrame([
    {"item": "n_genes", "value": len(ge)},
    {"item": "neutral_threshold_abs_delta", "value": thr},
    {"item": "n_concordant", "value": int((ge["class"] == "concordant").sum())},
    {"item": "n_discordant", "value": int((ge["class"] == "discordant").sum())},
    {"item": "n_neutral", "value": int((ge["class"] == "neutral").sum())},
    {"item": "spearman_delta_vs_betaR", "value": float(rho)},
    {"item": "spearman_descriptive_p", "value": float(rp)},
    {"item": "ZEB1_side_discordant_or_neutral", "value": int(((ge.side == "ZEB1-side50") & (ge["class"] != "concordant")).sum())},
    {"item": "ZEB2_side_discordant_or_neutral", "value": int(((ge.side == "ZEB2-side50") & (ge["class"] != "concordant")).sum())},
])
a62.to_csv(out / "a6_2_lock.tsv", sep="\t", index=False)
print(ge["class"].value_counts().to_string(), flush=True)
print(a62.to_string(index=False), flush=True)

# A6-3 split lists for cameraPR
z2g = ge.loc[ge.side == "ZEB2-side50"]
early_conc = z2g.loc[(z2g["class"] == "concordant") & (z2g["delta_normal"] < 0), "symbol"]
# concordant + Δ<0 is the normal-early class; also include concordant with Δ<0 only
# ZEB2-side β is negative, so concordant means Δ<0 automatically
other = z2g.loc[~z2g.symbol.isin(set(early_conc)), "symbol"]
pd.DataFrame({"symbol": early_conc}).to_csv(out / "a6_3_zeb2_early_concordant.tsv", sep="\t", index=False)
pd.DataFrame({"symbol": other}).to_csv(out / "a6_3_zeb2_discordant_or_neutral.tsv", sep="\t", index=False)
print("A6-3 early_concordant", len(early_conc), "discordant_or_neutral", len(other), flush=True)
