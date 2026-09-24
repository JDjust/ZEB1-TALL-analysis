#!/usr/bin/env python3
"""Module 5 round-2: map to M2 thymus reference, ZEB1-core scores,
patient-level composition, CITE-seq ADT (no muon/WNN).
Biological n is patient. ZEB1 stays continuous.
"""
from __future__ import annotations

import gzip
import math
import shutil
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from scipy.stats import chi2_contingency, spearmanr

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
AN = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
ROOT = AN / "module5"
M2 = AN / "module2"
M4 = AN / "module4"
PROC, TAB = ROOT / "processed", ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)


def log(m):
    print(m, flush=True)


def fdr_bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    mask = np.isfinite(p)
    if not mask.any():
        return out
    pv = p[mask]
    n = len(pv)
    order = np.argsort(pv)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    q = pv * n / ranked
    q_sorted = np.minimum.accumulate(q[order[::-1]])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(q_sorted, 0, 1)
    out[mask] = adj
    return out


def write_cores():
    st = pd.read_csv(M4 / "tables" / "M4_r2_subtype_adjusted_ZEB1.tsv", sep="\t")
    pos = st[(st["fdr_adj"] < 0.05) & (st["t_adj"] > 0)].sort_values("t_adj", ascending=False)
    neg = st[(st["fdr_adj"] < 0.05) & (st["t_adj"] < 0)].sort_values("t_adj")
    pos_g = pos["gene"].astype(str).head(80).tolist()
    neg_g = neg["gene"].astype(str).head(80).tolist()
    (PROC / "m4_adj_pos_core.txt").write_text("\n".join(pos_g), encoding="utf-8")
    (PROC / "m4_adj_neg_core.txt").write_text("\n".join(neg_g), encoding="utf-8")
    log(f"cores pos={len(pos_g)} neg={len(neg_g)}")
    return pos_g, neg_g


def score_h5ad(pos_g, neg_g):
    h5 = PROC / "GSE227122_annotated.h5ad"
    if not h5.exists():
        log("no GSE227122 h5ad")
        return pd.DataFrame()
    log(f"read {h5}")
    adata = sc.read_h5ad(h5)
    log(f"h5ad {adata.n_obs} x {adata.n_vars}")
    for name, genes in [("core_pos", pos_g), ("core_neg", neg_g)]:
        present = [g for g in genes if g in adata.var_names]
        log(f"  {name} present {len(present)}/{len(genes)}")
        if len(present) >= 8:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=False)
    keep = [c for c in ["patient", "timepoint", "sample", "lineage", "malignant"] if c in adata.obs]
    extra = [c for c in ["core_pos", "core_neg", "g_ZEB1"] if c in adata.obs.columns]
    if "g_ZEB1" not in adata.obs.columns and "ZEB1" in adata.var_names:
        vec = adata[:, "ZEB1"].X
        adata.obs["g_ZEB1"] = vec.toarray().ravel() if sparse.issparse(vec) else np.asarray(vec).ravel()
        extra.append("g_ZEB1")
    obs = adata.obs[keep + extra].copy()
    mal = obs.copy()
    if "malignant" in mal.columns:
        mal = mal[mal["malignant"].astype(str).str.lower().isin(["true", "1", "malignant", "yes"])]
    elif "lineage" in mal.columns:
        mal = mal[mal["lineage"].astype(str).eq("T")]
    dx = mal[mal["timepoint"].astype(str).eq("Dx")] if "timepoint" in mal.columns else mal
    gcols = [c for c in ["g_ZEB1", "core_pos", "core_neg"] if c in dx.columns]
    pb = dx.groupby("patient", observed=True)[gcols].mean().reset_index()
    pb["n_cells"] = dx.groupby("patient", observed=True).size().values
    pb.to_csv(TAB / "M5_r2_core_patient.tsv", sep="\t", index=False)
    log(f"core patient n={len(pb)}")
    return pb


def merge_mapping(pb: pd.DataFrame):
    mapp = M2 / "tables" / "M2_r2_map_GSE227122_patient.tsv"
    comp = M2 / "tables" / "M2_r2_map_GSE227122_stage_comp.tsv"
    old = TAB / "M5_5.36_patient_pseudobulk.tsv"
    if mapp.exists():
        shutil.copy2(mapp, TAB / "M5_r2_map_patient.tsv")
        mp = pd.read_csv(mapp, sep="\t")
    else:
        log("no M2 mapping table")
        mp = pd.DataFrame()
    if comp.exists():
        shutil.copy2(comp, TAB / "M5_r2_map_stage_comp.tsv")
        cp = pd.read_csv(comp, sep="\t")
    else:
        cp = pd.DataFrame()
    if old.exists():
        pb0 = pd.read_csv(old, sep="\t")
    else:
        pb0 = pd.DataFrame()
    d = mp.copy() if len(mp) else pd.DataFrame()
    if len(pb0):
        d = pb0.merge(d, on="patient", how="outer", suffixes=("", "_map")) if len(d) else pb0
    if len(pb):
        d = d.merge(pb, on="patient", how="outer", suffixes=("", "_core")) if len(d) else pb
    if not len(d):
        return d
    d.to_csv(TAB / "M5_r2_patient_merged.tsv", sep="\t", index=False)
    rows = []
    if "g_ZEB1" in d.columns:
        feats = [c for c in d.columns if c not in ("patient", "g_ZEB1") and pd.api.types.is_numeric_dtype(d[c])]
        for c in feats:
            x = pd.to_numeric(d["g_ZEB1"], errors="coerce")
            y = pd.to_numeric(d[c], errors="coerce")
            ok = x.notna() & y.notna()
            if int(ok.sum()) < 6:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            rows.append({"feature": c, "rho": float(rho), "p": float(p), "n": int(ok.sum())})
        out = pd.DataFrame(rows)
        if len(out):
            out["fdr"] = fdr_bh(out["p"].to_numpy())
            out.sort_values("p").to_csv(TAB / "M5_r2_ZEB1_vs_mapped.tsv", sep="\t", index=False)
    if len(cp) and "g_ZEB1" in d.columns:
        # coarse bins for composition (descriptive; n=10)
        d2 = d.dropna(subset=["g_ZEB1"]).copy()
        if len(d2) >= 6:
            d2["zeb_bin"] = pd.qcut(d2["g_ZEB1"], 2, labels=["ZEB1_low", "ZEB1_high"], duplicates="drop")
            coarse = {
                "DN1": "DN-like", "DN2": "DN-like", "DN3": "DN-like", "ISP": "DN-like",
                "DP_CD3neg": "DP-like", "DP_CD3pos": "DP-like",
                "CD4SP": "SP-like", "CD8SP": "SP-like",
            }
            cp2 = cp.merge(d2[["patient", "zeb_bin"]], on="patient", how="inner")
            if "nearest_stage" in cp2.columns:
                cp2["coarse"] = cp2["nearest_stage"].map(coarse).fillna("other")
                agg = cp2.groupby(["zeb_bin", "coarse"], observed=True)["n"].sum().unstack(fill_value=0)
                agg.to_csv(TAB / "M5_r2_composition_by_ZEB1bin.tsv", sep="\t")
                try:
                    chi, p, dof, _ = chi2_contingency(agg.to_numpy())
                    note = f"chi2={chi:.3f} dof={dof} p={p:.4g} n_patients={d2.patient.nunique()}\n"
                except Exception as e:
                    note = f"chi2 failed: {e}\n"
                (TAB / "M5_r2_composition_note.txt").write_text(note, encoding="utf-8")
                log(note.strip())
    return d


def cite_adt():
    """Patient-level ADT vs ZEB1 RNA. Combined RNA mtx is truncated; ADT+metadata used."""
    meta_p = DATA / "GSE248287/suppl/GSE248287_10x_metadata_all.txt.gz"
    adt_m = DATA / "GSE248287/suppl/GSE248287_10x_ADT_matrix.mtx.gz"
    adt_f = DATA / "GSE248287/suppl/GSE248287_10x_ADT_features.tsv.gz"
    bcs_p = DATA / "GSE248287/suppl/GSE248287_10x_barcodes.tsv.gz"
    if not all(p.exists() for p in (meta_p, adt_m, adt_f, bcs_p)):
        log("CITE combined files missing")
        return
    meta = pd.read_csv(meta_p, sep="\t")
    log(f"CITE meta {meta.shape} cols={list(meta.columns)[:25]}")
    meta.head(2).to_csv(TAB / "M5_r2_cite_meta_head.tsv", sep="\t", index=False)
    (TAB / "M5_r2_cite_meta_cols.txt").write_text("\n".join(map(str, meta.columns)), encoding="utf-8")
    feat = pd.read_csv(adt_f, sep="\t", header=None)
    bcs = pd.read_csv(bcs_p, sep="\t", header=None)[0].astype(str)
    log(f"ADT features {len(feat)} barcodes {len(bcs)}")
    try:
        mat = mmread(gzip.open(adt_m, "rb")).tocsr()
    except Exception as e:
        log(f"ADT mtx failed (likely truncated): {e}")
        (TAB / "M5_r2_cite_note.txt").write_text(f"ADT mtx unreadable: {e}\n", encoding="utf-8")
        return
    names = (feat[1] if feat.shape[1] >= 2 else feat[0]).astype(str).tolist()
    if mat.shape == (len(names), len(bcs)):
        x = mat.T
    elif mat.shape == (len(bcs), len(names)):
        x = mat
    else:
        log(f"ADT shape mismatch {mat.shape} vs {len(bcs)} x {len(names)}")
        (TAB / "M5_r2_cite_note.txt").write_text(
            f"ADT shape {mat.shape}; barcodes {len(bcs)}; features {len(names)}\n", encoding="utf-8"
        )
        return
    a = ad.AnnData(X=x.tocsr())
    a.var_names = pd.Index(names).astype(str)
    a.obs_names = bcs.astype(str)
    a.var_names_make_unique()
    # map metadata
    bc_col = meta.columns[0]
    meta["_bc"] = meta[bc_col].astype(str)
    imap_cols = [c for c in meta.columns if str(c).lower() in
                 ("sample", "patient", "orig.ident", "hash", "hto", "timepoint", "celltypes_all",
                  "celltypes_initial", "type", "disease", "ncount_rna") or
                 any(k in str(c).lower() for k in ("patient", "sample", "celltype", "hto", "hash"))]
    a.obs["_bc"] = a.obs_names.astype(str)
    for c in imap_cols:
        a.obs[str(c)] = a.obs["_bc"].map(dict(zip(meta["_bc"], meta[c])))
    # identify patient column
    pcol = None
    for c in a.obs.columns:
        cl = str(c).lower()
        if "patient" in cl or cl in ("orig.ident", "sample", "hash"):
            nun = a.obs[c].nunique(dropna=True)
            if 3 <= nun <= 80:
                pcol = c
                break
    if pcol is None:
        # try PedTALL from barcodes/metadata strings
        for c in a.obs.columns:
            s = a.obs[c].astype(str)
            if s.str.contains("PedTALL|P\\d+", regex=True).mean() > 0.2:
                pcol = c
                break
    log(f"ADT patient col={pcol} n_obs={a.n_obs}")
    if pcol is None:
        (TAB / "M5_r2_cite_note.txt").write_text("no patient column in ADT metadata\n", encoding="utf-8")
        return
    a.obs["patient"] = a.obs[pcol].astype(str)
    # keep Dx-like if possible
    keep = a.obs["patient"].astype(str)
    mask = keep.str.contains("PedTALL|Dx|_P\\d+|P\\d+", regex=True, na=False)
    if mask.mean() > 0.05:
        a = a[mask].copy()
    markers = [g for g in ["CD1A", "CD1a", "CD7", "CD3", "CD3E", "CD4", "CD8", "CD8A",
                           "CD34", "CD5", "CD2", "CD38", "CD44", "CD45"] if g in a.var_names]
    # also case-insensitive
    lower = {str(v).upper(): v for v in a.var_names}
    for g in ["CD1A", "CD7", "CD3", "CD4", "CD8", "CD34", "CD5", "CD2"]:
        if g not in a.var_names and g in lower:
            markers.append(lower[g])
    markers = list(dict.fromkeys(markers))
    log(f"ADT markers {markers}")
    rows = []
    for p, sub in a.obs.groupby("patient"):
        idx = sub.index
        rec = {"patient": p, "n_cells": len(sub)}
        for g in markers:
            vec = a[idx, g].X
            rec[f"ADT_{g}"] = float(np.asarray(vec.mean()).ravel()[0]) if vec.size else np.nan
        rows.append(rec)
    adt_pb = pd.DataFrame(rows)
    # normalize patient id to P1..P15
    adt_pb["patient_std"] = adt_pb["patient"].astype(str).str.extract(r"(P\d+)", expand=False)
    adt_pb.to_csv(TAB / "M5_r2_cite_ADT_patient.tsv", sep="\t", index=False)
    rna = TAB / "M5_5.42_patient_pseudobulk.tsv"
    if rna.exists() and "patient_std" in adt_pb.columns:
        r = pd.read_csv(rna, sep="\t")
        r["patient_std"] = r["patient"].astype(str)
        m = adt_pb.merge(r, on="patient_std", how="inner")
        m.to_csv(TAB / "M5_r2_cite_ADT_RNA_merged.tsv", sep="\t", index=False)
        assoc = []
        if "g_ZEB1" in m.columns:
            for c in [c for c in m.columns if c.startswith("ADT_")]:
                x = pd.to_numeric(m["g_ZEB1"], errors="coerce")
                y = pd.to_numeric(m[c], errors="coerce")
                ok = x.notna() & y.notna()
                if int(ok.sum()) < 6:
                    continue
                rho, p = spearmanr(x[ok], y[ok])
                assoc.append({"feature": c, "rho": float(rho), "p": float(p), "n": int(ok.sum())})
        pd.DataFrame(assoc).to_csv(TAB / "M5_r2_cite_ZEB1_vs_ADT.tsv", sep="\t", index=False)
        log(f"CITE merged patients {len(m)}")
    (TAB / "M5_r2_cite_note.txt").write_text(
        "WNN skipped: muon not installed. Combined RNA mtx is truncated. "
        "Patient-level ADT vs ZEB1 RNA used instead.\n",
        encoding="utf-8",
    )


def main():
    pos, neg = write_cores()
    pb = score_h5ad(pos, neg)
    merge_mapping(pb)
    try:
        cite_adt()
    except Exception as e:
        log(f"CITE failed: {type(e).__name__}: {e}")
        (TAB / "M5_r2_cite_note.txt").write_text(f"CITE failed: {e}\n", encoding="utf-8")
    log("M5 R2 PREPARE DONE")


if __name__ == "__main__":
    main()
