#!/usr/bin/env python3
"""Module 6 prepare: ZEB1 locus chromatin (GSE110637 peaks, GSE243915 contacts)
and promoter methylation (GSE69954). Do not stream-extract the 23 GB tar to disk;
count ZEB1/LMO2-overlapping contacts while reading members from tar.
"""
from __future__ import annotations

import gzip
import io
import re
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6")
PROC = ROOT / "processed"
TAB = ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

# generous windows, hg38 (GSE243915 2023) and hg19 (450k / older ChIP)
LOCI = {
    "ZEB1_hg38": ("chr10", 31318447, 31529814),
    "ZEB1_prom_hg38": ("chr10", 31316447, 31328447),
    "LMO2_hg38": ("chr11", 33858410, 33862300),
    "GATA3_hg38": ("chr10", 8053438, 8075198),
    "ZEB1_hg19": ("chr10", 31607511, 31818878),
    "ZEB1_prom_hg19": ("chr10", 31605511, 31617511),
}


def log(m):
    print(m, flush=True)


def overlaps(chrom, pos, name):
    c, a, b = LOCI[name]
    return chrom == c and a <= pos <= b


def extract_110637_peaks():
    tar_p = DATA / "GSE110637/suppl/GSE110637_RAW.tar"
    out = PROC / "gse110637"
    out.mkdir(parents=True, exist_ok=True)
    want = ("ATAC", "H3K4me1_peaks.bed", "H3K4me3_peaks.bed")
    with tarfile.open(tar_p, "r") as tf:
        for m in tf.getmembers():
            base = Path(m.name).name
            if not any(k in base for k in want):
                continue
            dest = out / base
            if dest.exists() and dest.stat().st_size > 100:
                continue
            src = tf.extractfile(m)
            dest.write_bytes(src.read())
            log(f"extracted {base}")
    rows = []
    for bed in sorted(out.glob("*.bed.gz")):
        kind = "ATAC" if "ATAC" in bed.name else ("H3K4me1" if "H3K4me1" in bed.name else "H3K4me3")
        n_tot = n_prom = n_gene = 0
        with gzip.open(bed, "rt") as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                p = line.split()
                if len(p) < 3:
                    continue
                chrom, start, end = p[0], int(p[1]), int(p[2])
                n_tot += 1
                mid = (start + end) // 2
                # try hg19 then hg38
                if overlaps(chrom, mid, "ZEB1_prom_hg19") or overlaps(chrom, mid, "ZEB1_prom_hg38"):
                    n_prom += 1
                if overlaps(chrom, mid, "ZEB1_hg19") or overlaps(chrom, mid, "ZEB1_hg38"):
                    n_gene += 1
        rows.append({"file": bed.name, "mark": kind, "n_peaks": n_tot, "ZEB1_promoter": n_prom, "ZEB1_genebody": n_gene})
        log(f"{bed.name} peaks={n_tot} ZEB1prom={n_prom} gene={n_gene}")
    pd.DataFrame(rows).to_csv(TAB / "M6_6.1_6.12_ZEB1_peak_overlap.tsv", sep="\t", index=False)


def count_243915_contacts():
    tar_p = DATA / "GSE243915/suppl/GSE243915_RAW.tar"
    meta = [
        ("GSM7798256_TALL-p01_filtered.reg.txt.gz", "TALL", "p01"),
        ("GSM7798257_TALL-p03_filtered.reg.txt.gz", "TALL", "p03"),
        ("GSM7798258_ETP-p04_filtered.reg.txt.gz", "ETP", "p04"),
        ("GSM7798259_TALL-p06_filtered.reg.txt.gz", "TALL", "p06"),
        ("GSM7798260_TALL-p11_filtered.reg.txt.gz", "TALL", "p11"),
        ("GSM7798261_ETP-p15_filtered.reg.txt.gz", "ETP", "p15"),
        ("GSM7798262_ETP-p16_filtered.reg.txt.gz", "ETP", "p16"),
        ("GSM7798263_CD34-2190_filtered.reg.txt.gz", "CD34", "2190"),
        ("GSM7798264_CD34-2583_filtered.reg.txt.gz", "CD34", "2583"),
        ("GSM7798265_THY-2123_filtered.reg.txt.gz", "THY", "2123"),
        ("GSM7798266_THY-2124_filtered.reg.txt.gz", "THY", "2124"),
    ]
    # skip the two giant thymus files first; do them last
    meta = [x for x in meta if x[1] != "THY"] + [x for x in meta if x[1] == "THY"]
    recs = []
    with tarfile.open(tar_p, "r") as tf:
        names = {Path(m.name).name: m for m in tf.getmembers()}
        for fname, group, sid in meta:
            if fname not in names:
                log(f"missing {fname}")
                continue
            log(f"stream {fname}")
            n = n_z = n_l = n_g = n_zz = 0
            fh = gzip.GzipFile(fileobj=tf.extractfile(names[fname]))
            wrap = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
            for line in wrap:
                p = line.split()
                if len(p) < 8:
                    continue
                try:
                    c1, pos1, c2, pos2 = p[1], int(p[3]), p[5], int(p[7])
                except (ValueError, IndexError):
                    continue
                n += 1
                z1 = overlaps(c1, pos1, "ZEB1_prom_hg38") or overlaps(c1, pos1, "ZEB1_hg38")
                z2 = overlaps(c2, pos2, "ZEB1_prom_hg38") or overlaps(c2, pos2, "ZEB1_hg38")
                l1 = overlaps(c1, pos1, "LMO2_hg38")
                l2 = overlaps(c2, pos2, "LMO2_hg38")
                g1 = overlaps(c1, pos1, "GATA3_hg38")
                g2 = overlaps(c2, pos2, "GATA3_hg38")
                if z1 or z2:
                    n_z += 1
                if l1 or l2:
                    n_l += 1
                if g1 or g2:
                    n_g += 1
                if (z1 and l2) or (z2 and l1):
                    n_zz += 1  # ZEB1-LMO2 contacts
                if n % 5000000 == 0:
                    log(f"  {fname} {n/1e6:.1f}M contacts ZEB1={n_z}")
            recs.append({
                "file": fname, "group": group, "sample": sid, "n_contacts": n,
                "ZEB1_locus": n_z, "LMO2_locus": n_l, "GATA3_locus": n_g,
                "ZEB1_LMO2_pairs": n_zz,
                "ZEB1_per_M": 1e6 * n_z / n if n else np.nan,
            })
            log(f"  done n={n} ZEB1={n_z} LMO2={n_l} ZEB1-LMO2={n_zz}")
            pd.DataFrame(recs).to_csv(TAB / "M6_6.13_6.16_ZEB1_contacts.tsv", sep="\t", index=False)
    pd.DataFrame(recs).to_csv(TAB / "M6_6.13_6.16_ZEB1_contacts.tsv", sep="\t", index=False)


def methylation_69954():
    """ZEB1 450k promoter probes. GEO matrix.gz is truncated; fall back to RAW.tar."""
    probes = [
        "cg01472637", "cg02764521", "cg03707168", "cg08856250", "cg16209936",
        "cg17152776", "cg18081940", "cg23210950", "cg24678886", "cg25626449",
    ]
    mat = DATA / "GSE69954/suppl/GSE69954_matrix_raw.txt.gz"
    raw = DATA / "GSE69954/suppl/GSE69954_RAW.tar"
    try:
        header = pd.read_csv(mat, sep="\t", nrows=0)
        beta_cols = [c for c in header.columns if "AVG_Beta" in str(c)]
        use = ["ID_REF"] + beta_cols
        log(f"GSE69954 beta columns {len(beta_cols)}")
        d = pd.read_csv(mat, sep="\t", usecols=lambda c: c in use)
    except Exception as e:
        log(f"GSE69954 matrix.gz unreadable ({e}); listing RAW.tar")
        if raw.exists():
            with tarfile.open(raw, "r") as tf:
                names = tf.getnames()
            pd.DataFrame({"raw_member": names}).to_csv(TAB / "M6_6.17_GSE69954_RAW_members.tsv", sep="\t", index=False)
            log(f"RAW members {len(names)}")
        return
    d["ID_REF"] = d["ID_REF"].astype(str).str.replace('"', "", regex=False)
    sub = d[d["ID_REF"].isin(probes)].copy()
    if sub.empty:
        log("curated ZEB1 probes not found")
        return
    long = sub.melt(id_vars="ID_REF", var_name="sample", value_name="beta")
    long["sample"] = long["sample"].str.replace(".AVG_Beta", "", regex=False).str.replace('"', "", regex=False)
    long.to_csv(TAB / "M6_6.17_ZEB1_promoter_beta.tsv", sep="\t", index=False)
    summ = long.groupby("ID_REF")["beta"].agg(["mean", "median", "std", "count"]).reset_index()
    summ.to_csv(TAB / "M6_6.17_ZEB1_promoter_beta_summary.tsv", sep="\t", index=False)
    log(f"wrote ZEB1 promoter beta probes={sub.shape[0]} samples={sub.shape[1]-1}")


def main():
    try:
        extract_110637_peaks()
    except Exception as e:
        log(f"110637 failed: {e}")
    try:
        methylation_69954()
    except Exception as e:
        log(f"69954 failed: {e}")
    try:
        count_243915_contacts()
    except Exception as e:
        log(f"243915 failed: {e}")
    log("MODULE6 PREPARE DONE")


if __name__ == "__main__":
    main()
