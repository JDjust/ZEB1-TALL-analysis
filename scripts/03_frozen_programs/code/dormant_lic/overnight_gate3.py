#!/usr/bin/env python3
"""Finish the pre-specified first round after Gate 2 failed.

Gate 3: GSE287751 Lmo2 KO versus control, four author contexts.
Ask whether Zeb1 RNA, ZEB1 CollecTRI activity, a frozen dDN3 signature,
and T-lineage genes move in the predicted directions.

Counterfactual: the same DN3 gate and low-versus-high cycle split in WT
and preleukemic thymi from GSE260948. If the ZEB1 activity drop is as large
in normal DN3, the leukemic result is a developmental correlate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CODE = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/dormant_lic/code")
sys.path.insert(0, str(CODE))
from gse260948_gate1 import G2M_GENES, S_GENES, gene_vector, mouse_symbol, present  # noqa: E402

GATE1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/dormant_lic/gate1")
OUT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/dormant_lic/overnight")
NET = Path("/data-b/liangfuhua/projects/bioinfo_direction_screening_20260903/config/collectri_human.tsv")
SCRNA = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/processed/gse287751_scrna")
META = SCRNA / "author_cell_metadata.tsv"
CYCLE = set(S_GENES) | set(G2M_GENES)
TLIN = ["Tcf7", "Gata3", "Runx1", "Bcl11b", "Il7r"]
QUIESCENCE = ["Txnip", "Btg1", "Trib2", "Ccr7", "Btg2", "Klf6", "Pnrc1", "Zfp36", "Klf2"]


def log(msg: str) -> None:
    print(msg, flush=True)


def zeb1_regulon(symbols: set[str], species: str) -> pd.DataFrame:
    net = pd.read_csv(NET, sep="\t")
    stim = net["consensus_stimulation"].astype(str).str.lower().isin({"true", "1", "1.0"})
    inhib = net["consensus_inhibition"].astype(str).str.lower().isin({"true", "1", "1.0"})
    net = net.copy()
    net["weight"] = np.where(stim & ~inhib, 1.0, np.where(inhib & ~stim, -1.0, np.nan))
    zeb = net[net["source_genesymbol"].astype(str).str.upper().eq("ZEB1") & net["weight"].notna()]
    zeb = zeb[["target_genesymbol", "weight"]].drop_duplicates("target_genesymbol")
    if species == "mouse":
        zeb["target"] = zeb["target_genesymbol"].map(mouse_symbol)
    else:
        zeb["target"] = zeb["target_genesymbol"]
    zeb = zeb[zeb["target"].isin(symbols)].drop_duplicates("target")
    return zeb


def weighted(expr_by_sample: pd.DataFrame, reg: pd.DataFrame) -> pd.Series:
    """expr_by_sample is genes x samples. Activity is the weighted mean of log values."""
    use = reg[reg["target"].isin(expr_by_sample.index)]
    x = expr_by_sample.loc[use["target"]]
    w = use.set_index("target").loc[x.index, "weight"].astype(float).to_numpy()
    return pd.Series(np.asarray(x.T) @ w / np.sum(np.abs(w)), index=expr_by_sample.columns)


def load_ko() -> tuple[sc.AnnData, pd.DataFrame]:
    link = OUT / "gse287751_10x"
    link.mkdir(parents=True, exist_ok=True)
    for kind, srcname in [
        ("matrix.mtx.gz", "GSM9286835_cDNA_matrix.mtx.gz"),
        ("barcodes.tsv.gz", "GSM9286835_cDNA_barcodes.tsv.gz"),
        ("features.tsv.gz", "GSM9286835_cDNA_features.tsv.gz"),
    ]:
        dest = link / kind
        if not dest.exists():
            dest.symlink_to(SCRNA / srcname)
    adata = sc.read_10x_mtx(link, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    meta = pd.read_csv(META, sep="\t")
    meta = meta.set_index("index")
    common = adata.obs_names.intersection(meta.index)
    log(f"GSE287751 matrix={adata.n_obs} author={len(meta)} matched={len(common)}")
    adata = adata[common].copy()
    adata.obs = meta.loc[common].copy()
    parts = adata.obs["hto"].str.extract(r"(?P<time>preNotch|d3)_(?P<priming>flt3l|unprimed)_(?P<geno>KO|EV)")
    adata.obs = pd.concat([adata.obs, parts], axis=1)
    return adata, meta


def pseudobulk(adata: sc.AnnData) -> pd.DataFrame:
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if sparse.issparse(counts):
        counts = counts.tocsr()
    frames = []
    for hto, idx in adata.obs.groupby("hto").groups.items():
        pos = adata.obs_names.get_indexer(idx)
        block = counts[pos]
        total = np.asarray(block.sum(axis=0)).ravel()
        frames.append(pd.Series(total, index=adata.var_names, name=str(hto)))
    mat = pd.concat(frames, axis=1)
    cpm = np.log2(mat.div(mat.sum(axis=0), axis=1) * 1e6 + 1)
    return cpm


def signature_from_gate1() -> pd.DataFrame:
    rows = []
    for name in ["TALL1", "TALL2"]:
        obs = pd.read_csv(GATE1 / f"{name}_dn3_obs.tsv.gz", sep="\t", index_col=0)
        adata = sc.read_10x_mtx(GATE1 / "10x" / f"{name}_filt", var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        adata = adata[
            (adata.obs.n_genes_by_counts >= 200)
            & (adata.obs.n_genes_by_counts < 8000)
            & (adata.obs.pct_counts_mt < 20)
        ].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        common = obs.index.intersection(adata.obs_names)
        adata = adata[common]
        state = obs.loc[common, "state"]
        d = state.eq("dDN3").to_numpy()
        p = state.eq("pDN3").to_numpy()
        x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        delta = x[d].mean(axis=0) - x[p].mean(axis=0)
        rows.append(pd.Series(delta, index=adata.var_names, name=name))
        log(f"signature source {name} d={int(d.sum())} p={int(p.sum())}")
    both = pd.concat(rows, axis=1).dropna()
    both = both.loc[~both.index.isin(CYCLE | {"Zeb1", "Zeb2"})]
    same = np.sign(both["TALL1"]) == np.sign(both["TALL2"])
    both = both.loc[same & (both.min(axis=1).abs() >= 0.25)]
    both["mean_delta"] = both.mean(axis=1)
    both["arm"] = np.where(both["mean_delta"] > 0, "dDN3_high", "pDN3_high")
    both = both.reindex(both["mean_delta"].abs().sort_values(ascending=False).index)
    kept = []
    for arm in ["dDN3_high", "pDN3_high"]:
        kept.append(both.loc[both.arm.eq(arm)].head(50))
    out = pd.concat(kept)
    out.to_csv(OUT / "frozen_ddn3_signature.tsv", sep="\t")
    log(f"frozen signature up={(out.arm=='dDN3_high').sum()} down={(out.arm=='pDN3_high').sum()}")
    return out


def score_signature(expr: pd.DataFrame, sig: pd.DataFrame) -> pd.Series:
    up = [g for g in sig.index[sig.arm.eq("dDN3_high")] if g in expr.index]
    down = [g for g in sig.index[sig.arm.eq("pDN3_high")] if g in expr.index]
    if len(up) < 10 or len(down) < 10:
        return pd.Series(np.nan, index=expr.columns)
    return expr.loc[up].mean(axis=0) - expr.loc[down].mean(axis=0)


def gate3(sig: pd.DataFrame) -> pd.DataFrame:
    adata, _ = load_ko()
    # Store raw counts before normalization for pseudobulk.
    adata.layers["counts"] = adata.X.copy()
    cpm = pseudobulk(adata)
    cpm.to_csv(OUT / "gse287751_pseudobulk_log2cpm.tsv.gz", sep="\t")
    reg = zeb1_regulon(set(cpm.index), "mouse")
    reg.to_csv(OUT / "gse287751_zeb1_targets.tsv", sep="\t", index=False)
    act = weighted(cpm, reg)
    dorm = score_signature(cpm, sig)
    labels = adata.obs[["hto", "time", "priming", "geno"]].drop_duplicates().set_index("hto")
    tab = labels.copy()
    tab["Zeb1"] = cpm.loc["Zeb1"] if "Zeb1" in cpm.index else np.nan
    tab["Lmo2"] = cpm.loc["Lmo2"] if "Lmo2" in cpm.index else np.nan
    tab["ZEB1_activity"] = act
    tab["dormant_signature"] = dorm
    for gene in TLIN + QUIESCENCE:
        tab[gene] = cpm.loc[gene] if gene in cpm.index else np.nan
    tab.to_csv(OUT / "gse287751_condition_scores.tsv", sep="\t")
    contrasts = []
    for (time, priming), sub in tab.groupby(["time", "priming"]):
        if not {"KO", "EV"}.issubset(set(sub["geno"])):
            continue
        ko = sub.loc[sub.geno.eq("KO")].iloc[0]
        ev = sub.loc[sub.geno.eq("EV")].iloc[0]
        for feature in ["Lmo2", "Zeb1", "ZEB1_activity", "dormant_signature", *TLIN]:
            contrasts.append({
                "time": time,
                "priming": priming,
                "feature": feature,
                "KO": float(ko[feature]),
                "EV": float(ev[feature]),
                "delta_KO_minus_EV": float(ko[feature] - ev[feature]),
            })
    con = pd.DataFrame(contrasts)
    con.to_csv(OUT / "gse287751_contrasts.tsv", sep="\t", index=False)
    return con


def counterfactual(sig: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name in ["WT", "PL", "TALL1", "TALL2"]:
        path = GATE1 / "10x" / f"{name}_filt"
        if name in {"WT", "PL"} and not (path / "matrix.mtx.gz").exists():
            # WT and PL were copied during Gate 1 preparation.
            log(f"missing {name}")
            continue
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        adata = adata[
            (adata.obs.n_genes_by_counts >= 200)
            & (adata.obs.n_genes_by_counts < 8000)
            & (adata.obs.pct_counts_mt < 20)
        ].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        cd4, cd8, cd44 = gene_vector(adata, "Cd4"), gene_vector(adata, "Cd8a"), gene_vector(adata, "Cd44")
        il2ra, ptcra, dntt = gene_vector(adata, "Il2ra"), gene_vector(adata, "Ptcra"), gene_vector(adata, "Dntt")
        gate = (cd4 < 0.5) & (cd8 < 0.5) & (cd44 < 0.1) & ((il2ra > 0) | (ptcra > 0) | (dntt > 0))
        dn3 = adata[gate].copy()
        from gse260948_gate1 import cell_cycle
        cell_cycle(dn3)
        low = dn3.obs["phase"].eq("G1")
        high = dn3.obs["phase"].isin(["S", "G2M"])
        reg = zeb1_regulon(set(dn3.var_names), "mouse")
        # Per-cell activity inside this library.
        genes = reg["target"].tolist()
        w = reg.set_index("target").loc[genes, "weight"].astype(float).to_numpy()
        x = dn3[:, genes].X
        if sparse.issparse(x):
            x = x.toarray()
        z = (x - x.mean(0)) / (x.std(0) + 1e-8)
        act = z @ w / np.sum(np.abs(w))
        zeb = gene_vector(dn3, "Zeb1")
        for label, mask in [("G1", low.to_numpy()), ("SG2M", high.to_numpy())]:
            rows.append({
                "library": name,
                "state": label,
                "n": int(mask.sum()),
                "ZEB1_activity": float(np.nanmean(act[mask])) if mask.any() else np.nan,
                "Zeb1_expr": float(np.nanmean(zeb[mask])) if mask.any() else np.nan,
                "Zeb1_detected": float(np.nanmean(zeb[mask] > 0)) if mask.any() else np.nan,
            })
        if low.any() and high.any():
            rows.append({
                "library": name,
                "state": "G1_minus_SG2M",
                "n": int(low.sum() + high.sum()),
                "ZEB1_activity": float(np.nanmean(act[low]) - np.nanmean(act[high])),
                "Zeb1_expr": float(np.nanmean(zeb[low]) - np.nanmean(zeb[high])),
                "Zeb1_detected": float(np.nanmean(zeb[low] > 0) - np.nanmean(zeb[high] > 0)),
            })
        log(f"counterfactual {name} DN3={dn3.n_obs} G1={int(low.sum())} SG2M={int(high.sum())}")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "developmental_counterfactual.tsv", sep="\t", index=False)
    return tab


def decide(con: pd.DataFrame, cf: pd.DataFrame) -> None:
    def deltas(feature: str) -> list[float]:
        return con.loc[con.feature.eq(feature), "delta_KO_minus_EV"].tolist()

    zeb = deltas("Zeb1")
    act = deltas("ZEB1_activity")
    dorm = deltas("dormant_signature")
    lmo = deltas("Lmo2")
    lines = [
        "Gate 3 unit is one pseudobulk per author condition, not cells.",
        f"Lmo2 KO-EV deltas: {lmo}",
        f"Zeb1 KO-EV deltas: {zeb}",
        f"ZEB1 activity KO-EV deltas: {act}",
        f"dormant signature KO-EV deltas: {dorm}",
    ]
    zeb_up = sum(v > 0 for v in zeb)
    act_up = sum(v > 0 for v in act)
    dorm_down = sum(v < 0 for v in dorm)
    lines.append(f"Zeb1 up in {zeb_up}/{len(zeb)} contexts")
    lines.append(f"ZEB1 activity up in {act_up}/{len(act)} contexts")
    lines.append(f"dormant signature down in {dorm_down}/{len(dorm)} contexts")
    passed = zeb_up == len(zeb) and act_up == len(act) and dorm_down == len(dorm) and len(zeb) == 4
    lines.append("GATE3 " + ("PASS" if passed else "FAIL"))
    lines.append("Developmental counterfactual is G1 minus S/G2M inside the same DN3 marker gate.")
    lines.append(cf.to_string(index=False))
    text = "\n".join(lines) + "\n"
    (OUT / "overnight_decision.txt").write_text(text)
    log(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sig = signature_from_gate1()
    con = gate3(sig)
    cf = counterfactual(sig)
    decide(con, cf)


if __name__ == "__main__":
    main()
