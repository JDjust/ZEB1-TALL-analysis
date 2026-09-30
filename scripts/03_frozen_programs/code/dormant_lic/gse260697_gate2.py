#!/usr/bin/env python3
"""Gate 2: ZEB1 activity in human CD7+CD1a- L-IC versus CD7+CD1a+ blasts.

REL13 and REL14 are two PDX models. Libraries inside a model are not patients.
Primary endpoint is the CollecTRI ZEB1 weighted-target score. RNA is auxiliary.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE260697/suppl")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/dormant_lic/gate2")
NET = Path("/data-b/liangfuhua/projects/bioinfo_direction_screening_20260903/config/collectri_human.tsv")
CYCLE = {
    "MKI67", "TOP2A", "PCNA", "MCM5", "MCM2", "MCM4", "MCM6", "TYMS", "RRM2",
    "CDK1", "CCNB2", "UBE2C", "BIRC5", "NUSAP1", "CENPF", "AURKA", "AURKB",
}


def log(msg: str) -> None:
    print(msg, flush=True)


def compartment(column: str) -> str:
    name = column.upper()
    if "CD1AN" in name or "CD1A_MINUS" in name or "CD1N" in name or "MINUS" in name:
        return "LIC"
    if "CD1AP" in name or "CD1A_PLUS" in name or "CD1P" in name or "PLUS" in name:
        return "DP"
    raise ValueError(column)


def load_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df[df["Type"].eq("protein_coding")].drop_duplicates("Name")
    mat = df.set_index("Name").filter(regex="COUNT$").astype(float)
    mat.columns = [c.replace(":COUNT", "") for c in mat.columns]
    return mat


def log_cpm(mat: pd.DataFrame) -> pd.DataFrame:
    lib = mat.sum(axis=0)
    return np.log2(mat.div(lib, axis=1) * 1e6 + 1)


def zeb1_targets(symbols: set[str]) -> pd.DataFrame:
    net = pd.read_csv(NET, sep="\t")
    stim = net["consensus_stimulation"].astype(str).str.lower().isin({"true", "1", "1.0"})
    inhib = net["consensus_inhibition"].astype(str).str.lower().isin({"true", "1", "1.0"})
    net = net.copy()
    net["weight"] = np.where(stim & ~inhib, 1.0, np.where(inhib & ~stim, -1.0, np.nan))
    zeb = net[net["source_genesymbol"].astype(str).str.upper().eq("ZEB1") & net["weight"].notna()]
    zeb = zeb[["target_genesymbol", "weight"]].drop_duplicates("target_genesymbol")
    zeb = zeb[zeb["target_genesymbol"].isin(symbols)]
    return zeb


def score(expr: pd.DataFrame, reg: pd.DataFrame) -> pd.Series:
    x = expr.loc[reg["target_genesymbol"]].T
    z = (x - x.mean()) / (x.std(ddof=0) + 1e-8)
    w = reg.set_index("target_genesymbol").loc[z.columns, "weight"].astype(float)
    return z.mul(w, axis=1).sum(axis=1) / np.sum(np.abs(w))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for model, fname in [("REL13", "GSE260697_REL13_COUNT.txt.gz"), ("REL14", "GSE260697_REL14_COUNT.txt.gz")]:
        mat = load_counts(DATA / fname)
        expr = log_cpm(mat)
        labels = pd.Series({c: compartment(c) for c in expr.columns}, name="compartment")
        reg = zeb1_targets(set(expr.index))
        reg.to_csv(OUT / f"{model}_zeb1_targets.tsv", sep="\t", index=False)
        act = score(expr, reg)
        reg_nc = reg[~reg["target_genesymbol"].isin(CYCLE)]
        act_nc = score(expr, reg_nc)
        rna = expr.loc["ZEB1"]
        for sample in expr.columns:
            rows.append({
                "model": model,
                "sample": sample,
                "compartment": labels[sample],
                "ZEB1_log2cpm": float(rna[sample]),
                "ZEB1_activity": float(act[sample]),
                "ZEB1_activity_nocycle": float(act_nc[sample]),
                "n_targets": int(len(reg)),
                "n_targets_nocycle": int(len(reg_nc)),
            })
        log(f"{model} libraries={list(expr.columns)} targets={len(reg)}")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "library_scores.tsv", sep="\t", index=False)
    summary = []
    for model, sub in tab.groupby("model"):
        for feature in ["ZEB1_log2cpm", "ZEB1_activity", "ZEB1_activity_nocycle"]:
            a = sub.loc[sub.compartment.eq("LIC"), feature]
            b = sub.loc[sub.compartment.eq("DP"), feature]
            summary.append({
                "model": model,
                "feature": feature,
                "n_LIC": int(a.size),
                "n_DP": int(b.size),
                "mean_LIC": float(a.mean()),
                "mean_DP": float(b.mean()),
                "delta_LIC_minus_DP": float(a.mean() - b.mean()),
            })
    smry = pd.DataFrame(summary)
    smry.to_csv(OUT / "model_contrasts.tsv", sep="\t", index=False)
    print(tab.to_string(index=False), flush=True)
    print(smry.to_string(index=False), flush=True)
    act = smry[smry.feature.eq("ZEB1_activity")]
    signs = np.sign(act["delta_LIC_minus_DP"])
    decision = "PASS" if set(signs) == {-1.0} else "FAIL"
    (OUT / "gate2_decision.txt").write_text(decision + "\n")
    log(f"GATE2 {decision}")


if __name__ == "__main__":
    main()
