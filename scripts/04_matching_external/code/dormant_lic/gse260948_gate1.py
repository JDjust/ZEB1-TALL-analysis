#!/usr/bin/env python3
"""Gate 1 for the dormant L-IC hypothesis on GSE260948.

Reconstruct leukemic DN3 states in Tal1/Lmo2 T-ALL#1 and T-ALL#2 before
testing ZEB1. Zeb1 is excluded from the features used to define the state.

This is an independent Scanpy reconstruction. It is not the authors'
unpublished Seurat object. A state is accepted only if, in both leukemia
libraries, the low-cycle cluster also has lower Mki67, lower Myc, and
higher quiescence than the high-cycle cluster.

Primary ZEB1 endpoint: CollecTRI weighted-target activity, same direction
in both libraries. ZEB1 RNA is reported and is not the gate.
"""
from __future__ import annotations

import gzip
import shutil
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

RAW = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE260948/suppl/GSE260948_RAW.tar")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/dormant_lic/gate1")
TENX = OUT / "10x"
DECO = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/total_fig467F/pydeps")
MIN_UMI = 500
MAX_UMI = 50000

LIBS = {
    "GSM8128573": {"library": "WT", "genotype": "WT", "stage": "normal"},
    "GSM8128574": {"library": "PL", "genotype": "Tal1_Lmo2", "stage": "preleukemic"},
    "GSM8128575": {"library": "TALL1", "genotype": "Tal1_Lmo2", "stage": "leukemia"},
    "GSM8128576": {"library": "TALL2", "genotype": "Tal1_Lmo2", "stage": "leukemia"},
}
LEUK = ["TALL1", "TALL2"]
QUIESCENCE = ["Txnip", "Btg1", "Trib2", "Ccr7", "Btg2", "Klf6", "Pnrc1", "Zfp36", "Klf2"]
NOTCH = ["Notch1", "Dtx1", "Hes1", "Il2ra"]
DN3_POS = ["Il2ra", "Ptcra", "Dntt", "Rag1", "Notch1"]
DN3_NEG = ["Cd4", "Cd8a", "Cd44"]
EXCLUDE_FROM_STATE = {"Zeb1", "Zeb2"}


def log(msg: str) -> None:
    print(msg, flush=True)


def mouse_symbol(gene: str) -> str:
    g = str(gene)
    if g.startswith("MT-"):
        return "mt-" + g[3:].lower()
    return g[:1].upper() + g[1:].lower()


def extract_10x() -> None:
    TENX.mkdir(parents=True, exist_ok=True)
    with tarfile.open(RAW, "r") as tf:
        for member in tf.getmembers():
            name = Path(member.name).name
            dest = TENX / name
            if dest.exists() and dest.stat().st_size > 100:
                continue
            src = tf.extractfile(member)
            if src is None:
                continue
            dest.write_bytes(src.read())
            log(f"extracted {name}")


def n_barcodes(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh)


def umi_per_barcode(mtx_gz: Path, n_bc: int) -> np.ndarray:
    colsum = np.zeros(n_bc, dtype=np.float64)
    header_seen = False
    with gzip.open(mtx_gz, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if not header_seen:
                header_seen = True
                continue
            if len(parts) < 3:
                continue
            b = int(parts[1]) - 1
            if 0 <= b < n_bc:
                colsum[b] += float(parts[2])
    return colsum


def write_filtered(mtx_gz: Path, bc_gz: Path, ft_gz: Path, out_dir: Path, keep_idx: np.ndarray) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(bc_gz, "rt") as fh:
        barcodes = [ln.rstrip("\n") for ln in fh]
    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as fo:
        for i in keep_idx:
            fo.write(barcodes[int(i)] + "\n")
    shutil.copyfile(ft_gz, out_dir / "features.tsv.gz")
    mapping = {int(old): new + 1 for new, old in enumerate(keep_idx)}
    entries = []
    n_genes = None
    with gzip.open(mtx_gz, "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            gene, barcode, value = line.split()[:3]
            if n_genes is None:
                n_genes = int(gene)
                continue
            bi = int(barcode) - 1
            if bi in mapping:
                entries.append((int(gene), mapping[bi], value))
    with gzip.open(out_dir / "matrix.mtx.gz", "wt") as fo:
        fo.write("%%MatrixMarket matrix coordinate real general\n")
        fo.write(f"{n_genes} {len(keep_idx)} {len(entries)}\n")
        for gene, barcode, value in entries:
            fo.write(f"{gene} {barcode} {value}\n")
    log(f"  wrote {out_dir.name} cells={len(keep_idx)} nnz={len(entries)}")


def prepare_libraries() -> dict[str, Path]:
    extract_10x()
    paths = {}
    for gsm, meta in LIBS.items():
        mtx = next(TENX.glob(f"{gsm}_*_matrix.mtx.gz"))
        bcs = next(TENX.glob(f"{gsm}_*_barcodes.tsv.gz"))
        fts = next(TENX.glob(f"{gsm}_*_features.tsv.gz"))
        n_bc = n_barcodes(bcs)
        lib = TENX / f"{meta['library']}_filt"
        log(f"{meta['library']} barcodes={n_bc}")
        if n_bc > 30000:
            totals = umi_per_barcode(mtx, n_bc)
            keep = np.flatnonzero((totals >= MIN_UMI) & (totals < MAX_UMI))
            pd.DataFrame({"umi": totals}).describe().to_csv(OUT / f"{meta['library']}_umi_describe.tsv", sep="\t")
            log(f"  UMI {MIN_UMI}-{MAX_UMI}: {n_bc} -> {len(keep)}")
            if not (lib / "matrix.mtx.gz").exists():
                write_filtered(mtx, bcs, fts, lib, keep)
        else:
            lib.mkdir(exist_ok=True)
            for src, name in ((mtx, "matrix.mtx.gz"), (bcs, "barcodes.tsv.gz"), (fts, "features.tsv.gz")):
                tgt = lib / name
                if not tgt.exists():
                    shutil.copyfile(src, tgt)
        paths[meta["library"]] = lib
    return paths


def dense(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def gene_vector(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    return dense(adata[:, gene].X)


def present(adata, genes: list[str]) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def module_score(adata, genes: list[str], name: str) -> None:
    use = present(adata, genes)
    adata.obs[name] = np.nan
    adata.obs[f"{name}_n"] = len(use)
    if len(use) >= 3:
        sc.tl.score_genes(adata, use, score_name=name, random_state=1)


def qc_and_norm(adata):
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    n0 = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= 200)
        & (adata.obs["n_genes_by_counts"] < 8000)
        & (adata.obs["pct_counts_mt"] < 20)
    ].copy()
    log(f"  QC {n0} -> {adata.n_obs}")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


# Tirosh/Satija lists used by Seurat, with retired symbols mapped to current mouse genes.
S_GENES = [
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6",
    "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", "Hells", "Rfc2", "Rpa2", "Nasp",
    "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2",
    "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2",
    "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8",
]
G2M_GENES = [
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2",
    "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Pimreg", "Smc4", "Ccnb2",
    "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b",
    "Hjurp", "Cdca3", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2",
    "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln",
    "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa",
]


def cell_cycle(adata) -> None:
    s = present(adata, S_GENES)
    g2m = present(adata, G2M_GENES)
    log(f"  cell-cycle genes S={len(s)} G2M={len(g2m)}")
    if len(s) < 10 or len(g2m) < 10:
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan
        adata.obs["phase"] = "NA"
        return
    sc.tl.score_genes_cell_cycle(adata, s_genes=s, g2m_genes=g2m)
    adata.uns["cycle_genes"] = sorted(set(s) | set(g2m))


def cluster(adata, key: str, resolution: float) -> None:
    work = adata.copy()
    drop = [g for g in EXCLUDE_FROM_STATE if g in work.var_names]
    if drop:
        work = work[:, ~work.var_names.isin(drop)].copy()
    sc.pp.highly_variable_genes(work, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(work, n_comps=30, mask_var="highly_variable", random_state=1)
    sc.pp.neighbors(work, n_pcs=20, random_state=1)
    try:
        sc.tl.leiden(work, resolution=resolution, key_added=key, random_state=1, flavor="igraph", n_iterations=2)
    except TypeError:
        sc.tl.leiden(work, resolution=resolution, key_added=key, random_state=1)
    adata.obs[key] = work.obs[key].values
    adata.obsm["X_pca"] = work.obsm["X_pca"]


def dn3_cluster_table(adata) -> pd.DataFrame:
    rows = []
    for cl, idx in adata.obs.groupby("leiden_all").groups.items():
        sub = adata[idx]
        def mean_gene(g):
            v = gene_vector(sub, g)
            return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
        pos = np.nanmean([mean_gene(g) for g in present(adata, DN3_POS)])
        neg = np.nanmean([mean_gene(g) for g in present(adata, DN3_NEG)])
        rows.append({
            "cluster": str(cl),
            "n": int(sub.n_obs),
            "dn3_score": float(pos - neg),
            "Il2ra": mean_gene("Il2ra"),
            "Cd4": mean_gene("Cd4"),
            "Cd8a": mean_gene("Cd8a"),
            "Cd44": mean_gene("Cd44"),
            "Mki67": mean_gene("Mki67"),
            "S_score": float(sub.obs["S_score"].mean()),
            "G2M_score": float(sub.obs["G2M_score"].mean()),
        })
    return pd.DataFrame(rows).sort_values("dn3_score", ascending=False)


def label_states(adata) -> str:
    """Label the DN3 subset. Returns the rule used."""
    if adata.n_obs < 80:
        adata.obs["state"] = "too_few"
        return "too_few"
    cluster(adata, "leiden_dn3", resolution=0.4)
    tab = []
    for cl, idx in adata.obs.groupby("leiden_dn3").groups.items():
        sub = adata.obs.loc[idx]
        tab.append((str(cl), len(idx), float(sub["S_score"].mean() + sub["G2M_score"].mean()), float(gene_vector(adata[idx], "Mki67").mean())))
    tab = sorted(tab, key=lambda r: r[2])
    if len(tab) < 2:
        adata.obs["state"] = np.where(adata.obs["phase"].eq("G1"), "dDN3_proxy", "pDN3_proxy")
        return "phase_split_because_one_cluster"
    low, high = tab[0][0], tab[-1][0]
    # Intermediate clusters, if any, stay unlabeled and are excluded from the contrast.
    state = np.full(adata.n_obs, "intermediate", dtype=object)
    state[adata.obs["leiden_dn3"].astype(str).eq(low).to_numpy()] = "dDN3"
    state[adata.obs["leiden_dn3"].astype(str).eq(high).to_numpy()] = "pDN3"
    adata.obs["state"] = state
    return f"leiden_low_cycle={low};high_cycle={high};n_clusters={len(tab)}"


def mean_diff(adata, col: str) -> dict:
    a = adata.obs.loc[adata.obs["state"].eq("dDN3"), col].astype(float)
    b = adata.obs.loc[adata.obs["state"].eq("pDN3"), col].astype(float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    return {
        "feature": col,
        "n_dDN3": int(a.size),
        "n_pDN3": int(b.size),
        "mean_dDN3": float(a.mean()) if len(a) else np.nan,
        "mean_pDN3": float(b.mean()) if len(b) else np.nan,
        "delta_d_minus_p": float(a.mean() - b.mean()) if len(a) and len(b) else np.nan,
    }


def load_zeb1_regulon(var_names) -> pd.DataFrame:
    import sys
    sys.path.insert(0, str(DECO))
    import decoupler as dc
    net = dc.op.collectri(organism="mouse", verbose=True)
    net.to_csv(OUT / "collectri_mouse.tsv.gz", sep="\t", index=False)
    zeb = net[net["source"].astype(str).str.lower().eq("zeb1")].copy()
    zeb["target"] = zeb["target"].map(mouse_symbol)
    zeb = zeb[zeb["target"].isin(var_names)]
    zeb = zeb.drop_duplicates("target")
    log(f"ZEB1 CollecTRI targets in matrix: {len(zeb)}")
    return zeb


def weighted_activity(adata, regulon: pd.DataFrame, name: str, drop_genes: set[str] | None = None) -> None:
    use = regulon.copy()
    if drop_genes:
        use = use[~use["target"].isin(drop_genes)]
    genes = [g for g in use["target"] if g in adata.var_names]
    adata.obs[name] = np.nan
    adata.obs[f"{name}_ntargets"] = len(genes)
    if len(genes) < 8:
        log(f"  {name}: only {len(genes)} targets")
        return
    w = use.set_index("target").loc[genes, "weight"].astype(float).to_numpy()
    x = adata[:, genes].X
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x, dtype=float)
    z = (x - x.mean(axis=0)) / (x.std(axis=0) + 1e-8)
    adata.obs[name] = z @ w / np.sum(np.abs(w))


def phenotype_ok(row_d: dict, row_p_mki: dict, row_myc: dict, row_q: dict) -> bool:
    return (
        row_d["delta_d_minus_p"] > 0
        and row_p_mki["delta_d_minus_p"] < 0
        and row_myc["delta_d_minus_p"] < 0
        and row_q["delta_d_minus_p"] > 0
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.verbosity = 1
    paths = prepare_libraries()
    leuk_objects = {}
    phenotype_rows = []
    rules = []
    for name in LEUK:
        log(f"=== {name} ===")
        adata = sc.read_10x_mtx(paths[name], var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata = qc_and_norm(adata)
        cell_cycle(adata)
        module_score(adata, QUIESCENCE, "quiescence")
        module_score(adata, NOTCH, "notch")
        cluster(adata, "leiden_all", resolution=0.6)
        tab = dn3_cluster_table(adata)
        tab.insert(0, "library", name)
        tab.to_csv(OUT / f"{name}_cluster_dn3_scores.tsv", sep="\t", index=False)
        # Cell-level DN3 gate on log1p(CP10k). The cluster score is diagnostic only:
        # Il2ra minus Cd4/Cd8a stayed positive in DP-like clusters because Il2ra is broadly detected.
        cd4 = gene_vector(adata, "Cd4")
        cd8 = gene_vector(adata, "Cd8a")
        cd44 = gene_vector(adata, "Cd44")
        il2ra = gene_vector(adata, "Il2ra")
        ptcra = gene_vector(adata, "Ptcra")
        dntt = gene_vector(adata, "Dntt")
        gate = (cd4 < 0.5) & (cd8 < 0.5) & (cd44 < 0.1) & ((il2ra > 0) | (ptcra > 0) | (dntt > 0))
        dn3 = adata[gate].copy()
        log(f"  DN3 gate {int(gate.sum())} / {adata.n_obs}")
        rule = label_states(dn3) if dn3.n_obs >= 80 else "too_few"
        if dn3.n_obs < 80:
            dn3.obs["state"] = "too_few"
        rules.append({
            "library": name,
            "rule": rule,
            "n_qc": int(adata.n_obs),
            "n_dn3": int(dn3.n_obs),
            "gate": "Cd4<0.5,Cd8a<0.5,Cd44<0.1,Il2ra|Ptcra|Dntt>0",
        })
        dn3.obs["G1_indicator"] = (dn3.obs["phase"] == "G1").astype(float)
        dn3.obs["Mki67_expr"] = gene_vector(dn3, "Mki67")
        dn3.obs["Myc_expr"] = gene_vector(dn3, "Myc")
        dn3.obs["Zeb1_expr"] = gene_vector(dn3, "Zeb1")
        dn3.obs["Zeb1_detected"] = (dn3.obs["Zeb1_expr"] > 0).astype(float)
        leuk_objects[name] = dn3
        dn3.obs.to_csv(OUT / f"{name}_dn3_obs.tsv.gz", sep="\t")

    # Phenotype gate, computed before ZEB1 is interpreted.
    pheno = []
    passed = {}
    for name, dn3 in leuk_objects.items():
        contrasts = {
            "G1_indicator": mean_diff(dn3, "G1_indicator"),
            "Mki67_expr": mean_diff(dn3, "Mki67_expr"),
            "Myc_expr": mean_diff(dn3, "Myc_expr"),
            "quiescence": mean_diff(dn3, "quiescence"),
            "notch": mean_diff(dn3, "notch"),
            "S_score": mean_diff(dn3, "S_score"),
        }
        for row in contrasts.values():
            row["library"] = name
            row["layer"] = "phenotype"
            pheno.append(row)
        passed[name] = phenotype_ok(contrasts["G1_indicator"], contrasts["Mki67_expr"], contrasts["Myc_expr"], contrasts["quiescence"])
        log(f"PHENOTYPE {name} {'PASS' if passed[name] else 'FAIL'}")

    both_states = all(passed.values())
    # ZEB1 only after the state definition is recorded. Still computed so a
    # failed phenotype reconstruction is visible rather than hidden.
    reg = None
    try:
        common = set.intersection(*(set(a.var_names) for a in leuk_objects.values()))
        reg = load_zeb1_regulon(common)
        reg.to_csv(OUT / "zeb1_collectri_targets_in_matrix.tsv", sep="\t", index=False)
    except Exception as exc:
        log(f"CollecTRI unavailable: {type(exc).__name__}: {exc}")
        (OUT / "collectri_error.txt").write_text(f"{type(exc).__name__}: {exc}\n")

    zeb_rows = []
    for name, dn3 in leuk_objects.items():
        cycle = set(dn3.uns.get("cycle_genes", []))
        if reg is not None:
            weighted_activity(dn3, reg, "ZEB1_activity")
            weighted_activity(dn3, reg, "ZEB1_activity_nocycle", drop_genes=cycle)
        for feature in ["Zeb1_expr", "Zeb1_detected", "ZEB1_activity", "ZEB1_activity_nocycle"]:
            if feature not in dn3.obs:
                continue
            row = mean_diff(dn3, feature)
            row["library"] = name
            row["layer"] = "ZEB1"
            # Within-G1 control. pDN3 G1 cells are a minority by the published definition.
            g1 = dn3[dn3.obs["phase"].eq("G1")].copy()
            if set(g1.obs["state"]).issuperset({"dDN3", "pDN3"}) and (g1.obs["state"] == "pDN3").sum() >= 20:
                g1row = mean_diff(g1, feature)
                g1row["library"] = name
                g1row["layer"] = "ZEB1_within_G1"
                zeb_rows.append(g1row)
            zeb_rows.append(row)
        dn3.obs.to_csv(OUT / f"{name}_dn3_obs.tsv.gz", sep="\t")

    pd.DataFrame(rules).to_csv(OUT / "state_rules.tsv", sep="\t", index=False)
    pd.DataFrame(pheno).to_csv(OUT / "phenotype_contrasts.tsv", sep="\t", index=False)
    pd.DataFrame(zeb_rows).to_csv(OUT / "zeb1_contrasts.tsv", sep="\t", index=False)

    activity = [r for r in zeb_rows if r["feature"] == "ZEB1_activity" and r["layer"] == "ZEB1"]
    nocycle = [r for r in zeb_rows if r["feature"] == "ZEB1_activity_nocycle" and r["layer"] == "ZEB1"]
    dirs = [np.sign(r["delta_d_minus_p"]) for r in activity if np.isfinite(r["delta_d_minus_p"])]
    dirs_nc = [np.sign(r["delta_d_minus_p"]) for r in nocycle if np.isfinite(r["delta_d_minus_p"])]
    gate = {
        "phenotype_both_libraries": bool(both_states),
        "phenotype_TALL1": bool(passed.get("TALL1", False)),
        "phenotype_TALL2": bool(passed.get("TALL2", False)),
        "activity_both_negative": bool(len(dirs) == 2 and set(dirs) == {-1.0}),
        "activity_nocycle_both_negative": bool(len(dirs_nc) == 2 and set(dirs_nc) == {-1.0}),
    }
    gate["gate1"] = "PASS" if gate["phenotype_both_libraries"] and gate["activity_both_negative"] else "FAIL"
    pd.Series(gate).to_csv(OUT / "gate1_decision.tsv", sep="\t", header=False)
    log("GATE1 " + str(gate))
    log("DONE")


if __name__ == "__main__":
    main()
