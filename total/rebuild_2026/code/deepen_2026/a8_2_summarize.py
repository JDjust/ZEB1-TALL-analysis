"""A8-2: merge ALLCatchR2 calls with A8-1 scores. Subtype context only.

Primary = high-confidence. Sensitivity = high-confidence + candidate.
TLX3 DP-like and TLX3 immature stay separate in primary reports.
TLX3 family merge is a pre-specified sensitivity only.
Rank: mapped labels with n>=5 and a unique Pölönen counterpart.
If two ALLCatchR2 labels map to the same Pölönen subtype they enter
rank only after the pre-specified family merge.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026")
out = P / "a8_adult"
pred = pd.read_csv(out / "allcatchr2_predictions.tsv", sep="\t")
sc = pd.read_csv(out / "a8_patient_scores.tsv", sep="\t")

sid = pred.columns[0]
pred = pred.rename(columns={sid: "sample"})
if pred["sample"].astype(str).str.startswith("GSM").any():
    sc["key"] = sc["geo_accession"]
else:
    sc["key"] = sc["title"]
m = sc.merge(pred, left_on="key", right_on="sample", how="left")
m.to_csv(out / "a8_2_merged.tsv", sep="\t", index=False)

HC_MAIN = "T-ALL main-cluster high-confidence"
HC_SUB = "T-ALL sub-cluster high-confidence"
HC_IMM = "T-ALL immature high-confidence"
CD_MAIN = "T-ALL main-cluster candidate"
CD_SUB = "T-ALL sub-cluster candidate"
CD_IMM = "T-ALL immature candidate"

# Pre-specified 1:1 map (A82_LOCK). Match the official inner label only.
# Compound clusters (C1, C9, C12, C13) are not forced onto one Pölönen class.
UNIQUE_INNER = {
    "BCL11B": "BCL11B",
    "TLX1": "TLX1",
    "SPI1": "SPI1",
    "KMT2A": "KMT2A",
    "TAL1 DP-like": "TAL1 DP-like",
    "TAL1 αβ-like": "TAL1 αβ-like",
    "LMO2 γδ-like": "LMO2 γδ-like",
    "STAG2/LMO2": "STAG2&LMO2",
    "HOXA9/10 TCR": "HOXA9 TCR",
    "HOXA9 TCR": "HOXA9 TCR",
    "ETP-like": "ETP-like",
    "NKX2-1 other": "NKX2-1",
    "NKX2-1 TCR": "NKX2-1",
    "TLX3 DP-like": "TLX3",
    "TLX3 immature": "TLX3",
}
SHARED_POLO = {"TLX3", "NKX2-1"}


def nonempty(s: pd.Series) -> pd.Series:
    s = s.fillna("").astype(str).str.strip().str.strip(";")
    return s[(s != "") & (s != "NA") & (s != "nan")]


def has_label(series: pd.Series, token: str) -> pd.Series:
    """Substring match for pre-specified unique tokens (BCL11B, TLX3 DP-like, ...)."""
    return series.fillna("").astype(str).str.contains(token, regex=False)


def has_exact_call(series: pd.Series, label: str) -> pd.Series:
    """True if the official ALLCatchR2 label is one of the semicolon parts."""

    def _ok(raw: object) -> bool:
        parts = [p.strip() for p in str(raw).split(";")]
        return label in parts

    return series.fillna("").map(_ok)


def explode_labels(series: pd.Series) -> list[str]:
    out_l = []
    for raw in series.fillna(""):
        for part in str(raw).split(";"):
            part = part.strip().strip(";")
            if part and part not in {"NA", "nan"}:
                out_l.append(part)
    return out_l


def pole(mask: pd.Series, name: str) -> dict:
    n = int(mask.sum())
    if n == 0:
        row = {"item": name, "n": 0, "R_median": np.nan, "Z2_median": np.nan, "Z1_median": np.nan}
        print(f"{name} n=0")
        return row
    sub = m.loc[mask]
    row = {
        "item": name,
        "n": n,
        "R_median": float(np.nanmedian(sub["R"])),
        "Z2_median": float(np.nanmedian(sub["ZEB2_side50"])),
        "Z1_median": float(np.nanmedian(sub["ZEB1_side50"])),
        "ages": ",".join(str(x) for x in sub["age"].tolist()) if "age" in sub.columns else "",
    }
    print(
        f"{name} n={n} R_median={row['R_median']:.3f} "
        f"Z1={row['Z1_median']:.3f} Z2={row['Z2_median']:.3f}"
    )
    return row


rows = []
print("n_merged", len(m), "n_pred_matched", int(m["sample"].notna().sum()) if "sample" in m.columns else 0)

for col, title in [
    (HC_MAIN, "HC_main"),
    (HC_SUB, "HC_sub"),
    (HC_IMM, "HC_immature"),
    (CD_MAIN, "CD_main"),
    (CD_SUB, "CD_sub"),
    (CD_IMM, "CD_immature"),
]:
    if col not in m.columns:
        print(title, "missing")
        continue
    vc = Counter(explode_labels(m[col]))
    print(title)
    for k, v in vc.most_common():
        print(f"  {v}\t{k}")
        rows.append({"item": f"count_{title}", "label": k, "n": v})

unclassified = pd.Series(True, index=m.index)
if HC_MAIN in m.columns:
    unclassified &= nonempty(m[HC_MAIN]).reindex(m.index).isna() | (
        m[HC_MAIN].fillna("").astype(str).str.strip().isin(["", "NA", "nan"])
    )
n_unclass = int(
    (
        m[HC_MAIN].fillna("").astype(str).str.strip().isin(["", "NA", "nan"])
        if HC_MAIN in m.columns
        else pd.Series(True, index=m.index)
    ).sum()
)
print("HC_main_empty", n_unclass)
rows.append({"item": "HC_main_empty", "n": n_unclass})

if HC_MAIN in m.columns:
    rows.append(pole(has_label(m[HC_MAIN], "BCL11B"), "HC_BCL11B"))
    rows.append(pole(has_label(m[HC_MAIN], "TLX3 DP-like"), "HC_TLX3_DPlike"))
    rows.append(pole(has_label(m[HC_MAIN], "TLX3 immature"), "HC_TLX3_immature"))
    rows.append(pole(has_label(m[HC_MAIN], "TLX3"), "HC_TLX3_family_sensitivity"))
    rows.append(pole(has_label(m[HC_MAIN], "TLX1"), "HC_TLX1"))
if HC_IMM in m.columns:
    rows.append(pole(has_label(m[HC_IMM], "ETP-like"), "HC_ETPlike"))
if HC_SUB in m.columns:
    rows.append(pole(has_label(m[HC_SUB], "CH-related"), "HC_CHrelated"))

def union_call(hc_col: str, cd_col: str, token: str) -> pd.Series:
    mask = pd.Series(False, index=m.index)
    if hc_col in m.columns:
        mask = mask | has_label(m[hc_col], token)
    if cd_col in m.columns:
        mask = mask | has_label(m[cd_col], token)
    return mask


rows.append(pole(union_call(HC_MAIN, CD_MAIN, "BCL11B"), "HCplusCand_BCL11B"))
rows.append(pole(union_call(HC_MAIN, CD_MAIN, "TLX3 DP-like"), "HCplusCand_TLX3_DPlike"))
rows.append(pole(union_call(HC_MAIN, CD_MAIN, "TLX3 immature"), "HCplusCand_TLX3_immature"))
rows.append(pole(union_call(HC_MAIN, CD_MAIN, "TLX3"), "HCplusCand_TLX3_family"))
rows.append(pole(union_call(HC_IMM, CD_IMM, "ETP-like"), "HCplusCand_ETPlike"))

disc = pd.read_csv(P / "a1_external_residual/polonen_discovery_subtype_medians.tsv", sep="\t")
disc = disc.set_index("subtype")["residual_median"]


def inner_label(lab: str) -> str:
    lab = str(lab).strip()
    if "(" in lab and lab.endswith(")"):
        return lab[lab.find("(") + 1 : -1]
    return lab


def mapped_target(label: str) -> str | None:
    inner = inner_label(label)
    if inner in UNIQUE_INNER:
        return UNIQUE_INNER[inner]
    if inner == "Immature T-ALL (ETP-like)" or label == "Immature T-ALL (ETP-like)":
        return "ETP-like"
    return None


def collect_unique_hc_medians() -> pd.DataFrame:
    """One row per official HC label that maps uniquely to Pölönen."""
    labels = explode_labels(m[HC_MAIN]) if HC_MAIN in m.columns else []
    if HC_IMM in m.columns:
        labels += explode_labels(m[HC_IMM])
    recs = []
    for lab in sorted(set(labels)):
        target = mapped_target(lab)
        if target is None:
            recs.append({"allcatchr": lab, "polonen": "", "shared_map": False, "n": 0})
            continue
        src = m[HC_IMM] if target == "ETP-like" and HC_IMM in m.columns else m[HC_MAIN]
        mask = has_exact_call(src, lab)
        recs.append(
            {
                "allcatchr": lab,
                "polonen": target,
                "shared_map": target in SHARED_POLO,
                "n": int(mask.sum()),
                "adult_R_median": float(np.nanmedian(m.loc[mask, "R"])) if mask.any() else np.nan,
            }
        )
    return pd.DataFrame(recs)


map_df = collect_unique_hc_medians() if HC_MAIN in m.columns else pd.DataFrame()
if len(map_df):
    map_df.to_csv(out / "a8_2_subtype_map.tsv", sep="\t", index=False)

rank_rows = map_df[(map_df.get("shared_map") == False) & (map_df.get("n", 0) >= 5)].copy() if len(map_df) else pd.DataFrame()
if len(rank_rows):
    rank_rows["disc_R_median"] = rank_rows["polonen"].map(disc)
    rank_rows = rank_rows.dropna(subset=["adult_R_median", "disc_R_median"])
    if len(rank_rows) >= 3:
        rho, p = spearmanr(rank_rows["adult_R_median"], rank_rows["disc_R_median"])
        print(f"PRIMARY_RANK n_subtypes={len(rank_rows)} rho={rho:.3f} p={p:.3g}")
        rows.append({"item": "primary_rank_n", "n": len(rank_rows), "rho": float(rho), "p": float(p)})
        for _, r in rank_rows.iterrows():
            rows.append(
                {
                    "item": "primary_rank_point",
                    "label": r["allcatchr"],
                    "polonen": r["polonen"],
                    "n": int(r["n"]),
                    "adult_R_median": float(r["adult_R_median"]),
                    "disc_R_median": float(r["disc_R_median"]),
                }
            )
    else:
        print(f"PRIMARY_RANK skipped: only {len(rank_rows)} mapped subtypes with n>=5")
        rows.append({"item": "primary_rank_n", "n": len(rank_rows), "rho": np.nan})
else:
    print("PRIMARY_RANK skipped: no unique mapped subtype with n>=5")
    rows.append({"item": "primary_rank_n", "n": 0, "rho": np.nan})

# Sensitivity: TLX3 family as one mapped point + unique maps
if HC_MAIN in m.columns:
    fam = has_label(m[HC_MAIN], "TLX3")
    sens = []
    if int(fam.sum()) >= 5 and "TLX3" in disc.index:
        sens.append(
            {
                "polonen": "TLX3",
                "n": int(fam.sum()),
                "adult_R_median": float(np.nanmedian(m.loc[fam, "R"])),
                "disc_R_median": float(disc["TLX3"]),
            }
        )
    if len(rank_rows):
        for _, r in rank_rows.iterrows():
            if r["polonen"] != "TLX3":
                sens.append(
                    {
                        "polonen": r["polonen"],
                        "n": int(r["n"]),
                        "adult_R_median": float(r["adult_R_median"]),
                        "disc_R_median": float(r["disc_R_median"]),
                    }
                )
    sens_df = pd.DataFrame(sens)
    if len(sens_df) >= 3:
        rho, p = spearmanr(sens_df["adult_R_median"], sens_df["disc_R_median"])
        print(f"SENS_RANK_Tlx3family n_subtypes={len(sens_df)} rho={rho:.3f} p={p:.3g}")
        rows.append({"item": "sens_rank_tlx3family_n", "n": len(sens_df), "rho": float(rho), "p": float(p)})
    else:
        print(f"SENS_RANK_Tlx3family skipped: {len(sens_df)} points")
        rows.append({"item": "sens_rank_tlx3family_n", "n": len(sens_df), "rho": np.nan})

pd.DataFrame(rows).to_csv(out / "a8_2_summary.tsv", sep="\t", index=False)
print("wrote", out / "a8_2_summary.tsv")
