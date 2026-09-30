#!/usr/bin/env python3
"""Gate 2: do frozen NOTCH signatures share a thymocyte window with ZEB1?

Signatures are not refit. A missing gene is omitted, not reweighted.
GSE195812 is one pooled library per FACS stage, so it can place a trajectory
and cannot support a within-stage donor correlation.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

OUT = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/_zeb1_forensic")
COLLECTRI = Path("/data-b/liangfuhua/projects/bioinfo_direction_screening_20260903/config/collectri_human.tsv")
MIN_CELLS = 20

DIRECT9 = ["CA3", "CD300C", "DTX1", "FGR", "GPR157", "HES4", "LOC100129447", "RLN2", "TMEM158"]
DIRECT_INDEP = [g for g in DIRECT9 if g not in {"HES4", "DTX1"}]
CANONICAL = ["HES1", "HES4", "HEY1", "DTX1", "NRARP"]
SETS = {
    "NOTCH_direct_9": DIRECT9,
    "NOTCH_direct_without_HES4_DTX1": DIRECT_INDEP,
    "canonical_NOTCH_TF": CANONICAL,
}

ORDER_195812 = ["DN1", "DN2", "DN3", "ISP", "DP_CD3neg", "DP_CD3pos", "CD4SP", "CD8SP"]
ORDER_142522 = ["CD34plus", "ISP", "EC", "LC", "SP4", "SP8"]
ORDER_PARK = ["DN", "DP", "CD4SP", "CD8SP"]
PARK_MAP = {
    "double negative thymocyte": "DN",
    "double-positive, alpha-beta thymocyte": "DP",
    "CD4-positive, alpha-beta T cell": "CD4SP",
    "CD8-positive, alpha-beta T cell": "CD8SP",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def symbols_of(adata) -> pd.Index:
    for col in ("gene_symbols", "feature_name", "gene_symbol"):
        if col in adata.var.columns:
            return pd.Index(adata.var[col].astype(str))
    return pd.Index(adata.var_names.astype(str))


def mean_or_cpm(adata, mask, gene_idx, lib_col):
    """Return one vector. Counts use sum / total n_counts * 1e4 then log1p. Log data use the mean."""
    sub = adata[mask, gene_idx]
    x = sub.X
    if sparse.issparse(x):
        total = np.asarray(x.sum(axis=0)).ravel()
        xmax = x.max()
    else:
        total = np.asarray(x).sum(axis=0)
        xmax = np.nanmax(x) if x.size else 0
    if xmax > 30 and lib_col in adata.obs.columns:
        lib = float(np.asarray(adata.obs[lib_col])[np.asarray(mask)].sum())
        if lib <= 0:
            return np.full(len(gene_idx), np.nan)
        return np.log1p(total / lib * 1e4)
    n = int(mask.sum())
    return total / n


def pseudobulk(adata, group_cols, genes, lib_candidates):
    sym = symbols_of(adata)
    present = [g for g in genes if g in set(sym)]
    gene_idx = [int(np.where(sym == g)[0][0]) for g in present]
    lib_col = next((c for c in lib_candidates if c in adata.obs.columns), None)
    # detect scale on a small slice
    probe = adata[: min(200, adata.n_obs), gene_idx[:1] if gene_idx else [0]]
    xmax = probe.X.max()
    use_counts = (xmax > 30) and lib_col is not None
    obs = adata.obs.reset_index(drop=True)
    obs["_row"] = np.arange(len(obs))
    rows = []
    for key, part in obs.groupby(group_cols, observed=True):
        if not isinstance(key, tuple):
            key = (key,)
        if len(part) < MIN_CELLS:
            continue
        mask = np.zeros(adata.n_obs, dtype=bool)
        mask[part["_row"].to_numpy()] = True
        if use_counts:
            vec = mean_or_cpm(adata, mask, gene_idx, lib_col)
        else:
            sub = adata[mask, gene_idx].X
            vec = np.asarray(sub.mean(axis=0)).ravel()
        rec = {c: v for c, v in zip(group_cols, key)}
        rec["n_cells"] = int(len(part))
        rec.update({g: float(vec[i]) for i, g in enumerate(present)})
        rows.append(rec)
    out = pd.DataFrame(rows)
    return out, present, "log1p_cpm_pseudobulk" if use_counts else "mean_of_stored_matrix"


def zmean(df, genes):
    present = [g for g in genes if g in df.columns]
    if len(present) < 2:
        return pd.Series(np.nan, index=df.index), present
    sub = df[present].astype(float)
    z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1), present


def stage_table(df, stage_col, order, score_cols):
    rows = []
    for stage in order:
        part = df[df[stage_col] == stage]
        if part.empty:
            continue
        rec = {"stage": stage, "n_profiles": int(len(part))}
        for c in score_cols:
            rec[c] = float(part[c].mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def residual_spearman(df, stage_col, a, b):
    d = df.dropna(subset=[a, b, stage_col]).copy()
    if d[stage_col].nunique() < 2 or len(d) < 8:
        return {"n": int(len(d)), "rho": np.nan, "p": np.nan}
    d["_a"] = d[a] - d.groupby(stage_col)[a].transform("mean")
    d["_b"] = d[b] - d.groupby(stage_col)[b].transform("mean")
    # drop stages with a single profile; their residual is zero and inflates n
    counts = d[stage_col].map(d[stage_col].value_counts())
    d = d[counts >= 3]
    if len(d) < 8:
        return {"n": int(len(d)), "rho": np.nan, "p": np.nan}
    r = spearmanr(d["_a"], d["_b"])
    return {"n": int(len(d)), "rho": float(r.statistic), "p": float(r.pvalue)}


def collectri_weights():
    net = pd.read_csv(COLLECTRI, sep="\t")
    z = net[net["source_genesymbol"].astype(str) == "ZEB1"]
    w = {}
    for _, row in z.iterrows():
        gene = str(row["target_genesymbol"])
        if int(row["consensus_stimulation"]) == 1:
            w[gene] = 1.0
        elif int(row["consensus_inhibition"]) == 1:
            w[gene] = -1.0
    return w


def activity(df, weights):
    present = [g for g in weights if g in df.columns]
    if len(present) < 15:
        return pd.Series(np.nan, index=df.index), present
    sub = df[present].astype(float)
    z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
    signed = z.mul([weights[g] for g in present], axis=1)
    return signed.mean(axis=1), present


def gse29959_zeb1():
    path = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE29959/GSE29959_series_matrix.txt.gz")
    print("GSE29959 sha256", sha256(path), "bytes", path.stat().st_size)
    expr = pd.read_csv(path, sep="\t", comment="!", quotechar='"')
    expr = expr.rename(columns={expr.columns[0]: "probe"})
    expr = expr[expr["probe"].ne("ID_REF")].set_index("probe")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    import gzip
    sym = {}
    annot = "/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE13159/GPL570.annot.gz"
    with gzip.open(annot, "rt", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("#") or line.startswith("!") or not line.startswith("ID\t") and header is None:
                if line.startswith("ID\t"):
                    header = line.rstrip("\n").split("\t")
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            gene = parts[2].split("///")[0].strip()
            if gene and gene not in {"---", "NA"}:
                sym[parts[0]] = gene
    ann = pd.Series(sym)
    use = expr.loc[expr.index.intersection(ann.index)].copy()
    use["symbol"] = ann.reindex(use.index).values
    gene = use.groupby("symbol").mean(numeric_only=True)
    lines = ["ALLSIL", "DND41", "HPBALL", "KOPTK1", "TALL-1"]
    roles = ["DMSO", "GSI", "ICN_DMSO", "ICN_GSI", "MYC_DMSO", "MYC_GSI"]
    cols = list(gene.columns)
    rows = []
    if "ZEB1" not in gene.index:
        print("ZEB1 absent")
        return
    for i, line in enumerate(lines):
        block = {roles[k]: cols[i * 6 + k] for k in range(6)}
        vals = {role: float(np.log2(max(gene.loc["ZEB1", block[role]], 1))) for role in roles}
        vals["line"] = line
        vals["GSI_minus_DMSO"] = vals["GSI"] - vals["DMSO"]
        vals["ICN_GSI_minus_GSI"] = vals["ICN_GSI"] - vals["GSI"]
        rows.append(vals)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "gate3_gse29959_ZEB1_rna.tsv", sep="\t", index=False)
    print(out.round(3).to_string(index=False))


def main():
    weights = collectri_weights()
    print("collectri ZEB1 targets", len(weights))
    frames = []
    coupling = []

    # GSE195812: pooled FACS libraries
    a = ad.read_h5ad("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/processed/GSE195812_qc.h5ad", backed="r")
    genes = sorted(set(sum(SETS.values(), [])) | {"ZEB1"} | set(weights))
    pb, present, how = pseudobulk(a, ["facs_stage"], genes, ["total_counts", "n_counts"])
    print("195812", how, "genes_hit", len(set(present) & set(DIRECT9)), "of", len(DIRECT9))
    for name, gs in SETS.items():
        pb[name], used = zmean(pb, gs)
        print(" ", name, "used", used)
    pb["ZEB1_RNA"] = pb["ZEB1"] if "ZEB1" in pb.columns else np.nan
    pb["ZEB1_activity"], used_act = activity(pb, weights)
    print("  activity targets", len(used_act))
    pb["dataset"] = "GSE195812"
    pb["stage"] = pb["facs_stage"]
    pb["unit"] = "pooled_library"
    frames.append(pb)
    a.file.close()

    # GSE206710: donor x CD3 fraction. Not a DN1-SP FACS axis.
    a = ad.read_h5ad("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/processed/GSE206710_qc.h5ad", backed="r")
    pb, present, how = pseudobulk(a, ["donor", "fraction"], genes, ["total_counts", "n_counts"])
    print("206710", how, "genes_hit", len(set(present) & set(DIRECT9)))
    for name, gs in SETS.items():
        pb[name], used = zmean(pb, gs)
        print(" ", name, "used", len(used))
    pb["ZEB1_RNA"] = pb["ZEB1"] if "ZEB1" in pb.columns else np.nan
    pb["ZEB1_activity"], used_act = activity(pb, weights)
    print("  activity targets", len(used_act))
    pb["dataset"] = "GSE206710"
    pb["stage"] = pb["fraction"]
    pb["unit"] = "donor_fraction"
    frames.append(pb)
    for score in ["NOTCH_direct_9", "ZEB1_RNA", "ZEB1_activity"]:
        if score not in pb.columns:
            continue
        stats = residual_spearman(pb, "fraction", "NOTCH_direct_9", score if score != "NOTCH_direct_9" else "ZEB1_RNA")
        if score == "NOTCH_direct_9":
            continue
        coupling.append({"dataset": "GSE206710", "score": score, **stats, "note": "within CD3 fraction, 3 donors"})
    a.file.close()

    # Park: donor x coarse thymocyte stage
    a = ad.read_h5ad("/data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F/Park_HTA/Human_T_cells.h5ad", backed="r")
    sym = symbols_of(a)
    print("park symbol example", list(sym[:3]), "ZEB1" in set(sym), "DTX1" in set(sym))
    keep = a.obs["cell_type"].isin(PARK_MAP).to_numpy()
    # subset cells before pseudobulk to limit memory
    idx = np.where(keep)[0]
    # read gene slice only
    wanted = [g for g in genes if g in set(sym)]
    gidx = [int(np.where(sym == g)[0][0]) for g in wanted]
    sub = a[idx, gidx].to_memory()
    sub.var_names = wanted
    sub.obs["stage"] = sub.obs["cell_type"].map(PARK_MAP).values
    a.file.close()
    pb, present, how = pseudobulk(sub, ["donor_id", "stage"], wanted, ["n_counts", "total_counts"])
    print("park", how, "profiles", len(pb), "genes_hit", len(set(present) & set(DIRECT9)))
    for name, gs in SETS.items():
        pb[name], used = zmean(pb, gs)
        print(" ", name, "used", used)
    pb["ZEB1_RNA"] = pb["ZEB1"] if "ZEB1" in pb.columns else np.nan
    pb["ZEB1_activity"], used_act = activity(pb, weights)
    print("  activity targets", len(used_act))
    pb["dataset"] = "ParkHTA"
    pb["unit"] = "donor_stage"
    frames.append(pb)
    for score in ["ZEB1_RNA", "ZEB1_activity"]:
        stats = residual_spearman(pb, "stage", "NOTCH_direct_9", score)
        coupling.append({"dataset": "ParkHTA", "score": score, "contrast": "NOTCH_direct_9", **stats, "note": "within DN/DP/CD4SP/CD8SP"})
        stats2 = residual_spearman(pb, "stage", "NOTCH_direct_without_HES4_DTX1", score)
        coupling.append({"dataset": "ParkHTA", "score": score, "contrast": "NOTCH_direct_without_HES4_DTX1", **stats2, "note": "within stage"})
        stats3 = residual_spearman(pb, "stage", "canonical_NOTCH_TF", score)
        coupling.append({"dataset": "ParkHTA", "score": score, "contrast": "canonical_NOTCH_TF", **stats3, "note": "reference only"})

    # GSE142522 bulk sorts, n=1
    raw = pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE142522/GSE142522_Table_S1.tsv.gz", sep="\t", header=None)
    header = [str(x).strip() if pd.notna(x) else "" for x in raw.iloc[6].tolist()]
    data = raw.iloc[7:].copy()
    data.columns = header
    data["gene"] = data["ID"].astype(str).str.strip().str.upper()
    sample_cols = [c for c in ORDER_142522 if c in data.columns]
    expr = data.set_index("gene")[sample_cols].apply(pd.to_numeric, errors="coerce")
    expr = expr[~expr.index.duplicated(keep="first")].T
    expr.index.name = "stage"
    pb = expr.reset_index()
    print("142522 genes_hit", [g for g in DIRECT9 if g in expr.columns])
    for name, gs in SETS.items():
        pb[name], used = zmean(pb, gs)
        print(" ", name, "used", used)
    pb["ZEB1_RNA"] = pb["ZEB1"] if "ZEB1" in pb.columns else np.nan
    pb["ZEB1_activity"], used_act = activity(pb, weights)
    print("  activity targets", len(used_act))
    pb["dataset"] = "GSE142522"
    pb["unit"] = "one_bulk_sort"
    pb["n_cells"] = 1
    frames.append(pb)

    allpb = pd.concat(frames, ignore_index=True)
    keep_cols = ["dataset", "unit", "stage", "donor_id", "donor", "fraction", "n_cells", "ZEB1_RNA", "ZEB1_activity"] + list(SETS)
    keep_cols = [c for c in keep_cols if c in allpb.columns]
    allpb[keep_cols].to_csv(OUT / "gate2_profiles.tsv", sep="\t", index=False)
    pd.DataFrame(coupling).to_csv(OUT / "gate2_stage_adjusted.tsv", sep="\t", index=False)
    print("--- stage means ---")
    for ds, order, stage_name in [
        ("GSE195812", ORDER_195812, "stage"),
        ("ParkHTA", ORDER_PARK, "stage"),
        ("GSE142522", ORDER_142522, "stage"),
        ("GSE206710", ["CD3neg", "CD3pos"], "stage"),
    ]:
        part = allpb[allpb.dataset == ds]
        print(ds)
        print(stage_table(part, stage_name, order, ["NOTCH_direct_9", "NOTCH_direct_without_HES4_DTX1", "canonical_NOTCH_TF", "ZEB1_RNA", "ZEB1_activity"]).round(3).to_string(index=False))
    print("--- coupling ---")
    print(pd.DataFrame(coupling).round(3).to_string(index=False))
    gse29959_zeb1()


if __name__ == "__main__":
    main()
