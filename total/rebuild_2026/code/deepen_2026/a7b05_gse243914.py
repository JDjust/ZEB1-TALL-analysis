"""A7-B0.5: GSE243914 blast-enriched primary T-ALL. Cell lines out."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026")
out = P / "a7b_malignant"
fr = pd.read_csv(P / "gate2_program/frozen_ZEB_side50.tsv", sep="\t")
tpm = pd.read_csv(out / "GSE243914_counts_pool_tpm.csv.gz")
tpm["symbol_u"] = tpm["gene"].astype(str).str.upper()
mat = tpm.groupby("symbol_u").mean(numeric_only=True)

# author titles from series matrix; pool columns ALL-pXX / IA-*
etp = ["ALL-p04", "ALL-p15", "ALL-p16", "ALL-p17", "ALL-p18", "IA-06_ribozero"]
non = [
    "ALL-p01", "ALL-p02", "ALL-p03", "ALL-p05", "ALL-p06", "ALL-p07", "ALL-p08",
    "ALL-p09", "ALL-p10", "ALL-p11", "ALL-p12", "ALL-p13", "ALL-p14",
    "ALL-p19", "ALL-p20", "ALL-p21", "ALL-p22", "IA-15_ribozero",
]
use = etp + non
assert all(c in mat.columns for c in use)

# log2(TPM+1)
X = np.log2(mat[use].clip(lower=0) + 1)

def zrow(s):
    v = s.to_numpy(float)
    sd = v.std(ddof=0)
    return (v - v.mean()) / sd if sd > 0 else np.full_like(v, np.nan)

need = ["ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"]
have_axis = {g: g in X.index for g in need}
axis = pd.DataFrame({g: zrow(X.loc[g]) if have_axis[g] else np.nan for g in need}, index=use)
axis["B"] = axis["ZEB1"] - axis["ZEB2"]
axis["D"] = axis["CD1A"] - (axis["CD34"] + axis["LYL1"]) / 2
axis["group"] = ["ETP" if s in etp else "nonETP" for s in use]

fr["symbol_u"] = fr["symbol"].astype(str).str.upper()
rec = fr[fr.symbol_u.isin(X.index)].copy()
z1 = rec.loc[rec.side == "ZEB1-side50", "symbol_u"]
z2 = rec.loc[rec.side == "ZEB2-side50", "symbol_u"]
Z = pd.DataFrame({g: zrow(X.loc[g]) for g in rec.symbol_u}, index=use)
axis["ZEB1_side50"] = Z[z1].mean(axis=1)
axis["ZEB2_side50"] = Z[z2].mean(axis=1)
axis.to_csv(out / "a7b05_sample_scores.tsv", sep="\t")

def mwu(col, alt):
    a = axis.loc[etp, col].to_numpy(float)
    b = axis.loc[non, col].to_numpy(float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    p = mannwhitneyu(a, b, alternative=alt).pvalue if len(a) and len(b) else np.nan
    return float(np.median(a)), float(np.median(b)), float(p)

rows = []
for col, alt in [
    ("ZEB2_side50", "greater"),
    ("ZEB1_side50", "less"),
    ("B", "less"),
]:
    med_e, med_n, p = mwu(col, alt)
    rows.append({
        "item": col, "ETP_median": med_e, "nonETP_median": med_n,
        "delta_ETP_minus_non": med_e - med_n, "mwu_onesided_p": p,
        "n_ETP": 6, "n_nonETP": 18,
    })
# gene-level Δ vs -β_R
ge = rec.copy()
ge["delta"] = [float(X.loc[g, etp].mean() - X.loc[g, non].mean()) for g in ge.symbol_u]
ge["neg_beta"] = -ge.beta_R
ge["expected"] = np.sign(ge["delta"]) == np.sign(ge["neg_beta"])
ge.loc[ge.delta == 0, "expected"] = False
rho, rp = spearmanr(ge["delta"], ge["neg_beta"])
ge.to_csv(out / "a7b05_genes.tsv", sep="\t", index=False)
sm = pd.DataFrame(rows)
sm.loc[len(sm)] = {
    "item": "recovered", "ETP_median": np.nan, "nonETP_median": np.nan,
    "delta_ETP_minus_non": np.nan, "mwu_onesided_p": np.nan,
    "n_ETP": int((rec.side == "ZEB1-side50").sum()),
    "n_nonETP": int((rec.side == "ZEB2-side50").sum()),
}
extra = pd.DataFrame([
    {"item": "recovered_Z1_Z2", "ETP_median": (rec.side == "ZEB1-side50").sum(),
     "nonETP_median": (rec.side == "ZEB2-side50").sum()},
    {"item": "gene_expected_frac", "ETP_median": float(ge.expected.mean()),
     "nonETP_median": float(ge.loc[ge.side == "ZEB1-side50", "expected"].mean()),
     "delta_ETP_minus_non": float(ge.loc[ge.side == "ZEB2-side50", "expected"].mean()),
     "mwu_onesided_p": float(rho)},
    {"item": "rho_delta_vs_negbeta_p", "mwu_onesided_p": float(rp)},
    {"item": "axis_genes_present", "ETP_median": float(all(have_axis.values()))},
])
pd.concat([sm, extra], ignore_index=True).to_csv(out / "a7b05_summary.tsv", sep="\t", index=False)
print(axis[["group", "B", "ZEB1_side50", "ZEB2_side50"]].to_string())
print("recovered", len(rec), "Z1", (rec.side == "ZEB1-side50").sum(), "Z2", (rec.side == "ZEB2-side50").sum())
print(sm.to_string(index=False))
print("gene expected", float(ge.expected.mean()), "rho", rho, "p", rp)
