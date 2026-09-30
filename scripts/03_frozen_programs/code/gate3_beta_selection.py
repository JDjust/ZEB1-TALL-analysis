#!/usr/bin/env python3
"""Gate 3. Signatures are defined on GSE151079 only, before any DepMap score.

beta-activated: log2 mean of CD34+CD4+ and ISP CD28+ exceeds CD34+CD1-/CD34+CD1+
by >= 1, and every beta sample exceeds every early sample.
DP-transition: same rule for DP CD3-/CD3+ versus the beta-selection samples.
ZEB1, the frozen NOTCH lists, and Hallmark G2M/E2F/MYC targets are removed
before the call.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

VAL = Path(r"d:\_bioinformation\ZEB1\total\data\validation")
COUNTS = VAL / "GSE151079_norm_counts.txt.gz"
GSE29959 = VAL / "GSE29959_series_matrix.txt.gz"
ANNOT = None  # filled on hulu; local GPL annotation is not required for this half

DIRECT9 = ["CA3", "CD300C", "DTX1", "FGR", "GPR157", "HES4", "LOC100129447", "RLN2", "TMEM158"]
NOTCH10 = ["HES1", "HES4", "HEY1", "DTX1", "NRARP", "MYC", "PTCRA", "NOTCH3", "IL7R", "SHQ1"]
CANON = ["HES1", "HES4", "HEY1", "DTX1", "NRARP"]
STAGES = ["CD34+CD1-", "CD34+CD1+", "CD34+CD4+", "ISP CD28+", "DP CD3-", "DP CD3+", "SP CD4+", "SP CD8+"]
EARLY = ["CD34+CD1-", "CD34+CD1+"]
BETA = ["CD34+CD4+", "ISP CD28+"]
DP = ["DP CD3-", "DP CD3+"]


def col(stage, rep):
    # Header uses "CD34+ CD4+" with a space; other CD34 labels do not.
    label = "CD34+ CD4+" if stage == "CD34+CD4+" else stage
    return f"{label} R{rep}"


def load_counts() -> pd.DataFrame:
    df = pd.read_csv(COUNTS, sep="\t", decimal=",")
    df["symbol"] = df["symbol"].astype(str).str.split(".").str[0]
    cols = [col(s, r) for s in STAGES for r in (1, 2)]
    use = df[["symbol"] + cols].dropna(subset=["symbol"])
    use = use[use["symbol"].ne("") & use["symbol"].ne("nan")]
    num = use.groupby("symbol")[cols].mean()
    return np.log2(num.clip(lower=0) + 1)


def exclusion() -> set[str]:
    banned = set(DIRECT9) | set(NOTCH10) | {"ZEB1"}
    for name in ["HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]:
        text = (VAL / f"{name}.txt").read_text(encoding="utf-8")
        banned.update(ln.strip() for ln in text.splitlines() if ln.strip() and not ln.startswith("HALLMARK") and not ln.startswith(">"))
    return banned


def separated(mat, high_stages, low_stages, banned):
    high = [col(s, r) for s in high_stages for r in (1, 2)]
    low = [col(s, r) for s in low_stages for r in (1, 2)]
    hi = mat[high]
    lo = mat[low]
    diff = hi.mean(axis=1) - lo.mean(axis=1)
    perfect = hi.min(axis=1) > lo.max(axis=1)
    keep = perfect & (diff >= 1) & ~mat.index.isin(banned)
    return mat.index[keep].tolist(), diff


def zmean(mat, genes):
    present = [g for g in genes if g in mat.index]
    if len(present) < 2:
        return pd.Series(np.nan, index=mat.columns), present
    sub = mat.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    return z.mean(axis=0), present


def trajectory(logc, genesets):
    rows = []
    for stage in STAGES:
        for rep in (1, 2):
            c = col(stage, rep)
            rec = {"stage": stage, "replicate": rep}
            for name, genes in genesets.items():
                present = [g for g in genes if g in logc.index]
                # z across the 16 sorted samples, then the value at this sample
                rec[name] = np.nan
                rec[name + "_n"] = len(present)
            rows.append(rec)
    # compute properly
    out_scores = {}
    for name, genes in genesets.items():
        score, present = zmean(logc, genes)
        out_scores[name] = (score, present)
    table = []
    for stage in STAGES:
        for rep in (1, 2):
            c = col(stage, rep)
            rec = {"stage": stage, "replicate": f"R{rep}"}
            for name, (score, present) in out_scores.items():
                rec[name] = float(score[c])
                rec[name + "_n"] = len(present)
            if "ZEB1" in logc.index:
                rec["ZEB1_log2"] = float(logc.loc["ZEB1", c])
            table.append(rec)
    return pd.DataFrame(table), out_scores


def gse29959_state(beta_genes, dp_genes):
    expr = pd.read_csv(GSE29959, sep="\t", comment="!", quotechar='"')
    expr = expr.rename(columns={expr.columns[0]: "probe"})
    expr = expr[expr["probe"].ne("ID_REF")].set_index("probe")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    # GPL570 symbols were resolved on hulu; repeat the local-free mapping by gene symbol
    # stored from the previous run if present, else skip until hulu.
    sym_path = VAL / "gpl570_gene_symbols.tsv"
    if not sym_path.exists():
        print("no local GPL570 map; GSE29959 scoring deferred")
        return
    sym = pd.read_csv(sym_path, sep="\t")
    # not used in this pass
    print("symbol map", len(sym))


def main():
    logc = load_counts()
    banned = exclusion()
    beta_genes, beta_diff = separated(logc, BETA, EARLY, banned)
    dp_genes, dp_diff = separated(logc, DP, BETA, banned)
    # disjoint
    dp_genes = [g for g in dp_genes if g not in set(beta_genes)]
    print("beta-activated", len(beta_genes), "head", beta_genes[:15])
    print("dp-transition", len(dp_genes), "head", dp_genes[:15])
    for g in ["PTCRA", "RAG1", "DNTT", "CD3E", "TRBC1", "ZEB1", "DTX1", "HES1"]:
        if g not in logc.index:
            print(g, "absent")
            continue
        early = np.mean([logc.loc[g, col(s, r)] for s in EARLY for r in (1, 2)])
        beta = np.mean([logc.loc[g, col(s, r)] for s in BETA for r in (1, 2)])
        dp = np.mean([logc.loc[g, col(s, r)] for s in DP for r in (1, 2)])
        print(f"  {g} early {early:.2f} beta {beta:.2f} dp {dp:.2f} d_beta {beta-early:.2f} in_beta {g in beta_genes} banned {g in banned}")
    sets = {
        "NOTCH_direct_9": DIRECT9,
        "canonical_NOTCH_TF": CANON,
        "beta_activated": beta_genes,
        "dp_transition": dp_genes,
    }
    table, _ = trajectory(logc, sets)
    table.to_csv(VAL / "gate3_gse151079_trajectory.tsv", sep="\t", index=False)
    print(table.round(3).to_string(index=False))
    pd.DataFrame({"signature": "beta_activated", "gene": beta_genes}).to_csv(VAL / "gate3_beta_activated_genes.tsv", sep="\t", index=False)
    pd.DataFrame({"signature": "dp_transition", "gene": dp_genes}).to_csv(VAL / "gate3_dp_transition_genes.tsv", sep="\t", index=False)
    # stage centroids on log2, developmental genes only, exclusions removed
    stage_mean = pd.DataFrame({s: logc[[col(s, 1), col(s, 2)]].mean(axis=1) for s in STAGES})
    vary = stage_mean.max(axis=1) - stage_mean.min(axis=1)
    keep = stage_mean.index[(vary >= 1) & ~stage_mean.index.isin(banned)]
    stage_mean.loc[keep].to_csv(VAL / "gate3_stage_centroids.tsv", sep="\t")
    print("centroid genes", len(keep), "banned", len(banned))


if __name__ == "__main__":
    main()
