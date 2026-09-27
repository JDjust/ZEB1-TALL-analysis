"""A7-C: freeze ZEB2 split from A7-A, then apply to existing leukemia."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

ROOT = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026")
OUT = ROOT / "a7c_decompose"
OUT.mkdir(exist_ok=True)

a7a = pd.read_csv(ROOT / "a7_lineage/scores/a7a_genes_min20.tsv", sep="\t")
z2 = a7a[a7a.side == "ZEB2-side50"].copy()
z2["subset"] = np.where(z2["delta_Myeloid_minus_T"] > 0, "myeloid_concordant", "non_myeloid_concordant")
z2["subset"] = np.where(z2["recovered"], z2["subset"], "unrecovered_A7A")
split = z2[["symbol", "gene_id", "beta_R", "recovered", "delta_Myeloid_minus_T", "subset"]].copy()
split.to_csv(OUT / "a7c_zeb2_split.tsv", sep="\t", index=False)
n_my = int(((split.subset == "myeloid_concordant")).sum())
n_nm = int((split.subset == "non_myeloid_concordant").sum())
print("FROZEN ZEB2 subsets myeloid", n_my, "non_myeloid", n_nm, flush=True)

my_sym = set(split.loc[split.subset == "myeloid_concordant", "symbol"].str.upper())
nm_sym = set(split.loc[split.subset == "non_myeloid_concordant", "symbol"].str.upper())


def zrow(v):
    v = np.asarray(v, float)
    sd = v.std(ddof=0)
    return (v - v.mean()) / sd if sd > 0 else np.full_like(v, np.nan)


def score_from_log(X, genes):
    """X: genes x samples, already expression; return mean z across samples."""
    have = [g for g in genes if g in X.index]
    if not have:
        return pd.Series(np.nan, index=X.columns), 0
    Z = pd.DataFrame({g: zrow(X.loc[g].to_numpy()) for g in have}, index=X.columns)
    return Z.mean(axis=1), len(have)


rows = []

# --- GSE248287 ---
lcpm = pd.read_csv(ROOT / "a7b_malignant/a7b_patient_logcpm.tsv", sep="\t")
lcpm["symbol_u"] = lcpm["symbol"].astype(str).str.upper()
X248 = lcpm.drop(columns=["symbol"]).groupby("symbol_u").mean(numeric_only=True)
sc248 = pd.read_csv(ROOT / "a7b_malignant/a7b_patient_scores.tsv", sep="\t").set_index("patient")
s_my, n1 = score_from_log(X248, my_sym)
s_nm, n2 = score_from_log(X248, nm_sym)
idx = sc248.index.intersection(s_my.index)
B = sc248.loc[idx, "B"]
r_my, p_my = spearmanr(s_my.loc[idx], B)
r_nm, p_nm = spearmanr(s_nm.loc[idx], B)
pd.DataFrame({
    "patient": idx,
    "B": B.to_numpy(),
    "ZEB2_myeloid_concordant": s_my.loc[idx].to_numpy(),
    "ZEB2_non_myeloid": s_nm.loc[idx].to_numpy(),
}).to_csv(OUT / "a7c_GSE248287_scores.tsv", sep="\t", index=False)
rows += [
    {"cohort": "GSE248287", "subset": "myeloid_concordant", "n_genes": n1,
     "stat": "spearman_vs_B", "estimate": float(r_my), "p": float(p_my),
     "direction_ok": bool(r_my < 0)},
    {"cohort": "GSE248287", "subset": "non_myeloid_concordant", "n_genes": n2,
     "stat": "spearman_vs_B", "estimate": float(r_nm), "p": float(p_nm),
     "direction_ok": bool(r_nm < 0)},
]

# --- GSE243914 ---
fr = pd.read_csv(ROOT / "gate2_program/frozen_ZEB_side50.tsv", sep="\t")
tpm = pd.read_csv(ROOT / "a7b_malignant/GSE243914_counts_pool_tpm.csv.gz")
tpm["symbol_u"] = tpm["gene"].astype(str).str.upper()
mat = tpm.groupby("symbol_u").mean(numeric_only=True)
etp = ["ALL-p04", "ALL-p15", "ALL-p16", "ALL-p17", "ALL-p18", "IA-06_ribozero"]
non = [
    "ALL-p01", "ALL-p02", "ALL-p03", "ALL-p05", "ALL-p06", "ALL-p07", "ALL-p08",
    "ALL-p09", "ALL-p10", "ALL-p11", "ALL-p12", "ALL-p13", "ALL-p14",
    "ALL-p19", "ALL-p20", "ALL-p21", "ALL-p22", "IA-15_ribozero",
]
use = etp + non
X914 = np.log2(mat[use].clip(lower=0) + 1)
s_my, n1 = score_from_log(X914, my_sym)
s_nm, n2 = score_from_log(X914, nm_sym)
pd.DataFrame({
    "sample": use,
    "group": ["ETP" if s in etp else "nonETP" for s in use],
    "ZEB2_myeloid_concordant": s_my.loc[use].to_numpy(),
    "ZEB2_non_myeloid": s_nm.loc[use].to_numpy(),
}).to_csv(OUT / "a7c_GSE243914_scores.tsv", sep="\t", index=False)
for name, ser, ng in [
    ("myeloid_concordant", s_my, n1),
    ("non_myeloid_concordant", s_nm, n2),
]:
    a = ser.loc[etp].to_numpy(float)
    b = ser.loc[non].to_numpy(float)
    p = mannwhitneyu(a, b, alternative="greater").pvalue
    rows.append({
        "cohort": "GSE243914", "subset": name, "n_genes": ng,
        "stat": "MWU_ETP_greater", "estimate": float(np.median(a) - np.median(b)),
        "p": float(p), "direction_ok": bool(np.median(a) > np.median(b)),
    })

# --- optional gene-level external ---
for fn, cohort in [
    ("a65_external/a65_GSE146901_genes.tsv", "GSE146901"),
    ("a65_external/a65_GSE234608_genes.tsv", "GSE234608"),
]:
    ge = pd.read_csv(ROOT / fn, sep="\t")
    ge["symbol_u"] = ge["symbol"].astype(str).str.upper()
    ge = ge.merge(split.assign(symbol_u=split.symbol.str.upper())[["symbol_u", "subset"]],
                  on="symbol_u", how="inner")
    ge = ge[ge.subset.isin(["myeloid_concordant", "non_myeloid_concordant"])]
    ge["neg_beta"] = -ge["beta_R"]
    ge["expected"] = np.sign(ge["delta_external"]) == np.sign(ge["neg_beta"])
    ge.loc[ge.delta_external == 0, "expected"] = False
    ge.to_csv(OUT / f"a7c_{cohort}_genes.tsv", sep="\t", index=False)
    for name, sub in ge.groupby("subset"):
        rho, rp = spearmanr(sub["delta_external"], sub["neg_beta"]) if len(sub) > 2 else (np.nan, np.nan)
        rows.append({
            "cohort": cohort, "subset": name, "n_genes": int(len(sub)),
            "stat": "gene_expected_frac", "estimate": float(sub.expected.mean()),
            "p": float(rp) if pd.notna(rp) else np.nan,
            "direction_ok": bool(sub.expected.mean() > 0.5) if len(sub) else False,
        })
        rows.append({
            "cohort": cohort, "subset": name, "n_genes": int(len(sub)),
            "stat": "rho_delta_vs_negbeta", "estimate": float(rho) if pd.notna(rho) else np.nan,
            "p": float(rp) if pd.notna(rp) else np.nan,
            "direction_ok": bool(rho > 0) if pd.notna(rho) else False,
        })

sm = pd.DataFrame(rows)
sm.to_csv(OUT / "a7c_summary.tsv", sep="\t", index=False)
print(sm.to_string(index=False), flush=True)
