#!/usr/bin/env python3
"""Module 3 revision: independent T-ALL subtype cohorts.

GSE62156 / GSE26713 (GPL570 series matrix) + GSE110636 RNA-seq + TARGET labels.
Uses gzip -dc because NCBI series_matrix files have trailing garbage that
breaks Python 3.12 gzip.
"""
from __future__ import annotations

import gzip
import re
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
AN = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
ROOT = AN / "module3"
PROC = ROOT / "processed"
M9 = AN / "module9" / "processed"
PROC.mkdir(parents=True, exist_ok=True)

GENES = ["ZEB1", "ZEB2", "LMO2", "TLX1", "TLX3", "TAL1", "LYL1", "MEF2C",
         "GATA3", "TCF7", "BCL11B", "NKX2-1"]
FALLBACK = {
    "ZEB1": {"212764_at", "210875_s_at", "212758_s_at", "239952_at"},
    "ZEB2": {"203603_s_at", "228337_at"},
    "LMO2": {"204249_s_at", "204248_at"},
    "TLX1": {"214455_at"},
    "TLX3": {"219905_at"},
    "TAL1": {"206283_s_at", "216954_x_at"},
    "LYL1": {"210423_s_at"},
    "MEF2C": {"209199_s_at", "207968_s_at"},
    "GATA3": {"209602_s_at", "209604_s_at"},
    "TCF7": {"205254_x_at"},
    "BCL11B": {"222895_s_at"},
    "NKX2-1": {"206632_s_at"},
}


def log(msg: str) -> None:
    print(msg, flush=True)


def gzip_dc_lines(path: Path):
    p = subprocess.Popen(
        ["gzip", "-dc", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert p.stdout is not None
    for line in p.stdout:
        yield line
    p.wait()


def parse_gpl570() -> dict[str, str]:
    path = DATA / "GSE13159" / "GPL570.annot.gz"
    probe2gene: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = None
        idx_id, idx_sym = 0, None
        for line in f:
            if line.startswith("#") or line.startswith("^") or (
                line.startswith("!") and "platform_table_begin" not in line
            ):
                continue
            if line.startswith("!platform_table_begin"):
                continue
            if header is None:
                parts = line.rstrip("\n").split("\t")
                header = parts
                for i, col in enumerate(parts):
                    pl = col.lower().strip()
                    if pl == "id":
                        idx_id = i
                    if pl == "gene symbol":
                        idx_sym = i
                continue
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_id, idx_sym or 0):
                continue
            pid = parts[idx_id].strip().strip('"')
            sym = parts[idx_sym].strip().strip('"') if idx_sym is not None else ""
            if not pid or not sym or sym in {"---", ""}:
                continue
            probe2gene[pid] = re.split(r"///|;|,", sym)[0].strip()
    log(f"GPL570 mapped probes {len(probe2gene)}")
    return probe2gene


def wanted_probes(probe2gene: dict[str, str]) -> dict[str, set[str]]:
    out = {g: set() for g in GENES}
    alias = {
        "ZEB1": {"ZEB1", "TCF8"},
        "ZEB2": {"ZEB2", "SIP1", "ZFHX1B"},
        "LMO2": {"LMO2"},
        "TLX1": {"TLX1", "HOX11"},
        "TLX3": {"TLX3", "HOX11L2"},
        "TAL1": {"TAL1", "SCL"},
        "LYL1": {"LYL1"},
        "MEF2C": {"MEF2C"},
        "GATA3": {"GATA3"},
        "TCF7": {"TCF7"},
        "BCL11B": {"BCL11B"},
        "NKX2-1": {"NKX2-1", "NKX2.1", "TTF1", "TTF-1"},
    }
    for probe, gene in probe2gene.items():
        gu = gene.upper()
        for g, als in alias.items():
            if gu in als:
                out[g].add(probe)
    for g, ps in FALLBACK.items():
        if len(out[g]) == 0:
            out[g] |= ps
            log(f"{g} using fallback probes")
        log(f"{g} probes n={len(out[g])}")
    return out


def parse_series_meta(path: Path) -> pd.DataFrame:
    titles, gsm, chars = [], [], []
    for line in gzip_dc_lines(path):
        if line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            gsm = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
        elif line.startswith("!series_matrix_table_begin"):
            break
    n = max(len(titles), len(gsm), max((len(c) for c in chars), default=0))
    df = pd.DataFrame({"gsm": (gsm + [""] * n)[:n], "title": (titles + [""] * n)[:n]})
    for i, row in enumerate(chars):
        df[f"char{i + 1}"] = (row + [""] * n)[:n]
    return df


def extract_gene_expr(path: Path, probe_map: dict[str, set[str]]) -> pd.DataFrame:
    keep = set().union(*probe_map.values())
    probe_rows = {}
    gsms = None
    in_table = False
    for line in gzip_dc_lines(path):
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if gsms is None:
            gsms = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            log(f"  matrix samples {len(gsms)}")
            continue
        pid = line.split("\t", 1)[0].strip().strip('"')
        if pid in keep:
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            probe_rows[pid] = np.array(
                [float(v) if v not in {"", "NA", "null"} else np.nan for v in vals],
                dtype=float,
            )
    log(f"  extracted probes {len(probe_rows)}")
    gene_mat = {}
    for g, probes in probe_map.items():
        mats = [probe_rows[p] for p in probes if p in probe_rows]
        gene_mat[g] = np.nanmean(np.vstack(mats), axis=0) if mats else np.full(len(gsms), np.nan)
    expr = pd.DataFrame(gene_mat, index=gsms)
    expr.index.name = "gsm"
    return expr


def map_gse62156(raw: str) -> tuple[str, str]:
    s = str(raw).lower()
    s = s.replace("molecular subgroup:", "").strip()
    level2 = s
    if "tlx1" in s or "tlx3" in s:
        return "TLX", level2
    if "tal" in s:
        return "TAL", level2
    if "hoxa13" in s or s.startswith("hox") or "mll-t" in s or "calm-af10" in s:
        return "HOXA", level2
    if "immature" in s:
        return "Immature", level2
    return "Other", level2


def map_gse26713(raw: str) -> tuple[str, str, str]:
    s = str(raw).replace("cytogenetics:", "").strip()
    sl = s.lower()
    if sl in {"bm"}:
        return "NormalBM", s, "normal"
    if sl in {"tlx1", "tlx3"}:
        return "TLX", s, "T-ALL"
    if sl in {"tal1", "tal2"}:
        return "TAL", s, "T-ALL"
    if sl == "tal2/lmo1":
        return "TAL", s, "T-ALL"
    if sl in {"lmo2", "lmo1"}:
        return "Immature", s, "T-ALL"
    if sl == "hoxa":
        return "HOXA", s, "T-ALL"
    if sl in {"nkx2-1", "nkx2.1"}:
        return "NKX", s, "T-ALL"
    return "Other", s, "T-ALL"


def map_gse110636(sg: str) -> tuple[str, str]:
    s = str(sg).strip().upper()
    if s == "IMM":
        return "Immature", s
    if s in {"TLX", "TAL", "HOXA"}:
        return s, s
    return "Other", s


def map_target(sg: str) -> tuple[str, str]:
    s = str(sg).strip()
    if s in {"LMO2/LYL1"}:
        return "Immature", s
    if s in {"TLX", "TAL", "HOXA"}:
        return s, s
    if s in {"NKX2_1", "NKX2-1"}:
        return "NKX", s
    return "Other", s


def prepare_gse62156(probe_map) -> pd.DataFrame:
    sm = DATA / "GSE62156/matrix/GSE62156_series_matrix.txt.gz"
    log("GSE62156")
    meta = parse_series_meta(sm)
    expr = extract_gene_expr(sm, probe_map)
    df = meta.merge(expr, left_on="gsm", right_index=True, how="left")
    raw = df.filter(like="char").apply(
        lambda r: next((str(x) for x in r if "molecular subgroup" in str(x).lower()), ""),
        axis=1,
    )
    mapped = raw.map(map_gse62156)
    df["level2"] = [x[1] for x in mapped]
    df["level1"] = [x[0] for x in mapped]
    df["disease"] = "T-ALL"
    df["cohort"] = "GSE62156"
    df["platform"] = "Affymetrix HG-U133 Plus 2.0"
    log(str(df["level1"].value_counts()))
    return df


def prepare_gse26713(probe_map) -> pd.DataFrame:
    sm = DATA / "GSE26713/matrix/GSE26713_series_matrix.txt.gz"
    log("GSE26713")
    meta = parse_series_meta(sm)
    expr = extract_gene_expr(sm, probe_map)
    df = meta.merge(expr, left_on="gsm", right_index=True, how="left")
    raw = df.filter(like="char").apply(
        lambda r: next((str(x) for x in r if "cytogenetics:" in str(x).lower()), ""),
        axis=1,
    )
    mapped = raw.map(map_gse26713)
    df["level1"] = [x[0] for x in mapped]
    df["level2"] = [x[1] for x in mapped]
    df["disease"] = [x[2] for x in mapped]
    df["cohort"] = "GSE26713"
    df["platform"] = "Affymetrix HG-U133 Plus 2.0"
    log(str(df.groupby(["disease", "level1"]).size()))
    return df


def prepare_gse110636() -> pd.DataFrame:
    log("GSE110636")
    cold = pd.read_csv(M9 / "gse110636_symbol_coldata.tsv", sep="\t")
    counts = pd.read_csv(M9 / "gse110636_symbol_counts.tsv.gz", sep="\t", index_col=0)
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    log2 = np.log2(cpm + 1)
    rows = []
    for _, r in cold.iterrows():
        gsm = r["gsm"]
        if gsm not in log2.columns:
            continue
        lv1, lv2 = map_gse110636(r.get("subgroup", ""))
        rec = {
            "gsm": gsm,
            "title": r.get("title", ""),
            "level1": lv1,
            "level2": lv2,
            "disease": "T-ALL",
            "cohort": "GSE110636",
            "platform": "Illumina RNA-seq",
        }
        for g in GENES:
            rec[g] = float(log2.loc[g, gsm]) if g in log2.index else np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    log(str(df["level1"].value_counts()))
    return df


def prepare_target() -> pd.DataFrame:
    t = pd.read_csv(PROC / "target_tall_subtype.tsv", sep="\t")
    t = t.copy()
    mapped = t["subtype"].map(map_target)
    t["level1"] = [x[0] for x in mapped]
    t["level2"] = [x[1] for x in mapped]
    t["disease"] = "T-ALL"
    t["cohort"] = "TARGET-ALL-P2"
    t["platform"] = "Illumina RNA-seq"
    t["gsm"] = t.get("usi", t.index.astype(str))
    for g in GENES:
        if g in t.columns:
            t[g] = np.log2(pd.to_numeric(t[g], errors="coerce") + 1)
        elif g.replace("-", "_") in t.columns:
            t[g] = np.log2(pd.to_numeric(t[g.replace("-", "_")], errors="coerce") + 1)
        else:
            t[g] = np.nan
    log(f"TARGET n={len(t)}")
    log(str(t["level1"].value_counts()))
    return t


def main():
    p2g = parse_gpl570()
    pmap = wanted_probes(p2g)
    a = prepare_gse62156(pmap)
    b = prepare_gse26713(pmap)
    c = prepare_gse110636()
    d = prepare_target()
    keep = ["cohort", "platform", "gsm", "title", "disease", "level1", "level2"] + GENES
    parts = []
    for df in (a, b, c, d):
        for col in keep:
            if col not in df.columns:
                df[col] = np.nan if col in GENES else ""
        parts.append(df[keep].copy())
    out = pd.concat(parts, ignore_index=True)
    out.to_csv(PROC / "independent_keygenes.tsv", sep="\t", index=False)
    counts = out[out["disease"].eq("T-ALL")].groupby(["cohort", "level1"]).size().unstack(fill_value=0)
    log("\nT-ALL Level1 counts:\n" + counts.to_string())
    counts.to_csv(PROC / "independent_level1_counts.tsv", sep="\t")
    log("wrote independent_keygenes.tsv n=" + str(len(out)))


if __name__ == "__main__":
    main()
