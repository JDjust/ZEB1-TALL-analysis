#!/usr/bin/env python3
"""Module 10 (optional): GSE260697 PDX L-IC (CD7+CD1a-) vs DP (CD7+CD1a+).

ZEB1 continuous. n is sorted populations within REL-13 / REL-14, not a large
patient series. Do not upgrade occupancy or LMO2 into causality.
"""
from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE260697")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10")
PROC, TAB = ROOT / "processed", ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

KEYS = ["ZEB1", "ZEB2", "LMO2", "LYL1", "TAL1", "TCF7", "GATA3", "BCL11B",
        "CD1A", "CD7", "CD34", "IL7R", "MEF2C", "DNTT", "RAG1", "MKI67"]
ENS = {
    "ENSG00000148516": "ZEB1", "ENSG00000169554": "ZEB2", "ENSG00000135363": "LMO2",
    "ENSG00000104903": "LYL1", "ENSG00000162367": "TAL1", "ENSG00000081059": "TCF7",
    "ENSG00000107485": "GATA3", "ENSG00000127152": "BCL11B", "ENSG00000168685": "IL7R",
    "ENSG00000174059": "CD34", "ENSG00000158473": "CD1A", "ENSG00000173762": "CD7",
    "ENSG00000081189": "MEF2C", "ENSG00000107447": "DNTT", "ENSG00000152804": "RAG1",
    "ENSG00000148773": "MKI67", "ENSG00000164330": "EBF1", "ENSG00000066336": "SPI1",
    "ENSG00000122026": "KIT", "ENSG00000126583": "CD44", "ENSG00000170312": "CDK1",
    "ENSG00000131747": "TOP2A", "ENSG00000132646": "PCNA", "ENSG00000134057": "CCNB1",
    "ENSG00000145386": "CCNA2",
}
ETP = ["CD34", "KIT", "IL7R", "LMO2", "LYL1", "MEF2C", "SPI1", "CD44"]
CYCLE = ["MKI67", "PCNA", "TOP2A", "CDK1", "CCNB1", "CCNA2"]


def log(m):
    print(m, flush=True)


def read_table(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene"})
    drop = [c for c in df.columns if c.lower() in ("type", "description", "gene_biotype", "name") and c != "gene"]
    if drop:
        df = df.drop(columns=drop, errors="ignore")
    df["gene"] = df["gene"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df["gene"] = df["gene"].str.split("|").str[0].str.strip()
    return df


def classify(col: str) -> str:
    s = col.upper().replace("-", "_").replace(":", "_").replace(" ", "_")
    if any(k in s for k in ("CD1AN", "CD1A_N", "CD1A_MINUS", "CD1N", "CD1A_NEG", "CD1A-")):
        return "LIC"
    if any(k in s for k in ("CD1AP", "CD1A_P", "CD1A_PLUS", "CD1P", "CD1A_POS", "CD1A+")):
        return "DP"
    return "other"


def main():
    files = sorted(DATA.rglob("*COUNT*.txt.gz")) + sorted(DATA.rglob("*FPKM*.txt.gz"))
    log("files " + ", ".join(p.name for p in files))
    rows = []
    for fp in files:
        df = read_table(fp)
        kind = "count" if "COUNT" in fp.name.upper() else "fpkm"
        patient = "REL13" if "REL13" in fp.name.upper() else ("REL14" if "REL14" in fp.name.upper() else fp.stem)
        log(f"{fp.name} genes={len(df)} samples={df.shape[1]-1} cols={list(df.columns[1:12])}")
        meta = pd.DataFrame({"sample": df.columns[1:]})
        meta["group"] = [classify(c) for c in meta["sample"]]
        meta["patient"] = patient
        meta["kind"] = kind
        meta.to_csv(TAB / f"M10_{patient}_{kind}_meta.tsv", sep="\t", index=False)
        log(str(meta.group.value_counts().to_dict()))
        gdf = df.set_index("gene")
        gdf = gdf.apply(pd.to_numeric, errors="coerce")
        gdf = gdf.groupby(gdf.index).mean()
        mapped = pd.Index([ENS.get(str(i).split(".")[0], str(i).upper()) for i in gdf.index])
        gdf.index = mapped
        gdf = gdf.groupby(gdf.index).mean()
        log(f"  genes head={list(gdf.index[:8])} ZEB1_in={int('ZEB1' in gdf.index)}")
        for g in KEYS:
            if g not in gdf.index:
                continue
            v = pd.to_numeric(gdf.loc[g], errors="coerce")
            for samp, val in v.items():
                rows.append({
                    "file": fp.name, "patient": patient, "kind": kind,
                    "sample": samp, "group": classify(str(samp)), "gene": g, "value": float(val),
                })
        # signature scores on log1p FPKM/counts
        mat = np.log1p(gdf.apply(pd.to_numeric, errors="coerce").fillna(0))
        for name, genes in [("ETP", ETP), ("cycle", CYCLE)]:
            present = [g for g in genes if g in mat.index]
            if len(present) < 3:
                continue
            z = mat.loc[present]
            z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0)
            sc = z.mean(axis=0)
            for samp, val in sc.items():
                rows.append({
                    "file": fp.name, "patient": patient, "kind": kind,
                    "sample": samp, "group": classify(str(samp)), "gene": f"sig_{name}",
                    "value": float(val) if np.isfinite(val) else np.nan,
                })
    long = pd.DataFrame(rows)
    long.to_csv(TAB / "M10_10.1_keygene_values.tsv", sep="\t", index=False)
    # LIC vs DP within patient, FPKM preferred
    tests = []
    if long.empty or "kind" not in long.columns:
        log("no long table; check gene IDs")
        pd.DataFrame(rows).to_csv(TAB / "M10_10.1_keygene_values.tsv", sep="\t", index=False)
        return
    sub = long[(long["kind"] == "fpkm") & (long["group"].isin(["LIC", "DP"]))]
    for (patient, gene), dd in sub.groupby(["patient", "gene"]):
        a = dd.loc[dd.group.eq("LIC"), "value"].to_numpy()
        b = dd.loc[dd.group.eq("DP"), "value"].to_numpy()
        if len(a) < 2 or len(b) < 2:
            p = np.nan
        else:
            try:
                p = ttest_ind(a, b, equal_var=False).pvalue
            except Exception:
                p = np.nan
        tests.append({
            "patient": patient, "gene": gene, "n_LIC": len(a), "n_DP": len(b),
            "mean_LIC": float(np.nanmean(a)) if len(a) else np.nan,
            "mean_DP": float(np.nanmean(b)) if len(b) else np.nan,
            "delta_LIC_minus_DP": float(np.nanmean(a) - np.nanmean(b)) if len(a) and len(b) else np.nan,
            "p": p,
        })
    testdf = pd.DataFrame(tests)
    testdf.to_csv(TAB / "M10_10.2_LIC_vs_DP.tsv", sep="\t", index=False)
    log("M10 prepare DONE")


if __name__ == "__main__":
    main()
