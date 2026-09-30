#!/usr/bin/env python3
"""Score frozen GSE151079 signatures on DepMap T-ALL and on GSE29959."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

VAL = Path("/tmp/gate3")
DEPMAP = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
MATRIX = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE29959/GSE29959_series_matrix.txt.gz")
ANNOT = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE13159/GPL570.annot.gz")
OUT = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/_zeb1_forensic")

GROUP = {
    "ACH-000981": ("DND-41", "dependent"),
    "ACH-000197": ("TALL-1", "dependent"),
    "ACH-000995": ("JURKAT", "dependent"),
    "ACH-000953": ("SUP-T1", "dependent"),
    "ACH-000937": ("PF-382", "nondependent"),
    "ACH-000519": ("PEER", "nondependent"),
    "ACH-000101": ("KE-37", "nondependent"),
    "ACH-001737": ("CCRF-HSB-2-DM", "nondependent"),
}


def gene_symbol(col):
    return str(col).split(" (")[0].strip()


def auc(values, dep, non):
    d = values.reindex(dep).to_numpy(float)
    n = values.reindex(non).to_numpy(float)
    wins = ties = 0
    for a in d:
        for b in n:
            if a > b:
                wins += 1
            elif a == b:
                ties += 1
    return (wins + 0.5 * ties) / (len(d) * len(n))


def load_depmap(genes):
    df = pd.read_csv(DEPMAP, low_memory=False)
    df = df[df["ModelID"].astype(str).isin(GROUP)].copy()
    df = df.rename(columns={c: gene_symbol(c) for c in df.columns if c != "ModelID"})
    df["ModelID"] = df["ModelID"].astype(str)
    mat = df.groupby("ModelID").mean(numeric_only=True)
    want = [g for g in genes if g in mat.columns]
    return mat.reindex(list(GROUP))[want]


def score_sets(mat, sets, dep, non):
    rows = []
    for name, genes in sets.items():
        present = [g for g in genes if g in mat.columns]
        sub = mat[present]
        z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
        score = z.mean(axis=1)
        rows.append({
            "signature": name,
            "n_defined": len(genes),
            "n_present": len(present),
            "auc_dep_higher": auc(score, dep, non),
            **{GROUP[i][0]: float(score.loc[i]) for i in GROUP},
        })
    return pd.DataFrame(rows)


def project(centroids):
    genes = [g for g in centroids.index if g != "symbol"]
    mat = load_depmap(list(centroids.index))
    shared = [g for g in centroids.index if g in mat.columns]
    cent = centroids.loc[shared]
    expr = mat[shared]
    rows = []
    for mid, (name, group) in GROUP.items():
        rec = {"CellLineName": name, "group": group}
        for stage in cent.columns:
            r = spearmanr(expr.loc[mid], cent[stage], nan_policy="omit")
            rec[stage] = float(r.statistic)
        rows.append(rec)
    return pd.DataFrame(rows)


def gpl_symbols():
    sym = {}
    with gzip.open(ANNOT, "rt", errors="replace") as fh:
        header = None
        for line in fh:
            if not line.startswith("ID\t") and header is None:
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene = parts[2].split("///")[0].strip()
            if gene and gene not in {"---", "NA"}:
                sym[parts[0]] = gene
    return sym


def gse29959(sets):
    expr = pd.read_csv(MATRIX, sep="\t", comment="!", quotechar='"')
    expr = expr.rename(columns={expr.columns[0]: "probe"})
    expr = expr[expr["probe"].ne("ID_REF")].set_index("probe").apply(pd.to_numeric, errors="coerce")
    sym = pd.Series(gpl_symbols())
    use = expr.loc[expr.index.intersection(sym.index)].copy()
    use["symbol"] = sym.reindex(use.index).to_numpy()
    gene = use.groupby("symbol").mean(numeric_only=True)
    gene = np.log2(np.clip(gene, 1, None))
    lines = ["ALLSIL", "DND41", "HPBALL", "KOPTK1", "TALL-1"]
    roles = ["DMSO", "GSI", "ICN_DMSO", "ICN_GSI", "MYC_DMSO", "MYC_GSI"]
    cols = list(gene.columns)
    rows = []
    for name, genes in sets.items():
        present = [g for g in genes if g in gene.index]
        sub = gene.loc[present]
        # z across the 30 arrays, then mean: a state score, not a refit
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
        score = z.mean(axis=0)
        for i, line in enumerate(lines):
            block = {roles[k]: score[cols[i * 6 + k]] for k in range(6)}
            rows.append({
                "signature": name,
                "n_genes": len(present),
                "line": line,
                **{k: float(v) for k, v in block.items()},
                "GSI_minus_DMSO": float(block["GSI"] - block["DMSO"]),
                "ICN_GSI_minus_GSI": float(block["ICN_GSI"] - block["GSI"]),
                "MYC_GSI_minus_GSI": float(block["MYC_GSI"] - block["GSI"]),
            })
    return pd.DataFrame(rows)


def main():
    beta = pd.read_csv(VAL / "gate3_beta_activated_genes.tsv", sep="\t")["gene"].tolist()
    dp = pd.read_csv(VAL / "gate3_dp_transition_genes.tsv", sep="\t")["gene"].tolist()
    sets = {"beta_activated": beta, "dp_transition": dp}
    dep = [i for i, (_, g) in GROUP.items() if g == "dependent"]
    non = [i for i, (_, g) in GROUP.items() if g == "nondependent"]
    genes = sorted(set(beta) | set(dp))
    mat = load_depmap(genes)
    scored = score_sets(mat, sets, dep, non)
    scored.to_csv(OUT / "gate3_depmap_auc.tsv", sep="\t", index=False)
    print(scored.round(3).to_string(index=False))
    cent = pd.read_csv(VAL / "gate3_stage_centroids.tsv", sep="\t", index_col=0)
    proj = project(cent)
    proj.to_csv(OUT / "gate3_stage_projection.tsv", sep="\t", index=False)
    print(proj.round(3).to_string(index=False))
    pert = gse29959(sets)
    pert.to_csv(OUT / "gate3_gse29959_state.tsv", sep="\t", index=False)
    print(pert.groupby("signature")[["GSI_minus_DMSO", "ICN_GSI_minus_GSI", "MYC_GSI_minus_GSI"]].mean().round(3))
    print(pert.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
