#!/usr/bin/env python3
"""Module 6 round-2: LMO2/TAL1/LDB1/GATA2 occupancy, ALL-SIL H3K27ac,
scATAC ZEB1 peaks (GSE280728; no ArchR), mouse thymocyte ATAC at Zeb1.
Do not interpret occupancy as repression.
"""
from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6")
M9 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module9/processed")
PROC, TAB = ROOT / "processed", ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

ZEB1_HG38 = ("chr10", 31_318_447, 31_529_814)
ZEB1_PROM38 = ("chr10", 31_316_447, 31_328_447)
ZEB1_HG19 = ("chr10", 31_607_511, 31_818_878)
ZEB1_PROM19 = ("chr10", 31_605_511, 31_617_511)
ZEB1_MM10 = ("chr18", 5_200_000, 5_560_000)
ZEB1_MM39 = ("chr18", 5_600_000, 5_850_000)


def log(m):
    print(m, flush=True)


def chrom_ok(c):
    c = str(c)
    if not c.startswith("chr"):
        c = "chr" + c
    return c


def overlaps(chrom, start, end, window):
    wchr, a, b = window
    chrom = chrom_ok(chrom)
    return chrom == wchr and start < b and end > a


def extract_member(tar_p: Path, dest: Path, pred):
    dest.mkdir(parents=True, exist_ok=True)
    out = []
    with tarfile.open(tar_p, "r") as tf:
        for m in tf:
            if not m.isfile():
                continue
            base = Path(m.name).name
            if not pred(base):
                continue
            target = dest / base
            if target.exists() and target.stat().st_size > 100:
                out.append(target)
                continue
            src = tf.extractfile(m)
            if src is None:
                continue
            target.write_bytes(src.read())
            log(f"extracted {base} {target.stat().st_size}")
            out.append(target)
    return out


def peak_overlap(path: Path, windows: dict):
    rec = {"file": path.name, "n_peaks": 0}
    for k in windows:
        rec[k] = 0
    hits = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line[0] in "#b":
                continue
            p = line.split()
            if len(p) < 3:
                continue
            try:
                chrom, start, end = p[0], int(float(p[1])), int(float(p[2]))
            except ValueError:
                continue
            rec["n_peaks"] += 1
            for k, w in windows.items():
                if overlaps(chrom, start, end, w):
                    rec[k] += 1
                    if k.startswith("ZEB1") and len(hits) < 12:
                        hits.append(f"{chrom}:{start}-{end}")
    rec["hits"] = ";".join(hits)
    return rec


def gse70734():
    tar = DATA / "GSE70734/suppl/GSE70734_RAW.tar"
    files = extract_member(tar, PROC / "gse70734", lambda n: "H3K27ac" in n or n.endswith(".txt.gz"))
    rows = []
    wins = {"ZEB1_gene_hg19": ZEB1_HG19, "ZEB1_prom_hg19": ZEB1_PROM19,
            "ZEB1_gene_hg38": ZEB1_HG38, "ZEB1_prom_hg38": ZEB1_PROM38}
    for f in files:
        rows.append(peak_overlap(f, wins))
        log(f"70734 {rows[-1]}")
    pd.DataFrame(rows).to_csv(TAB / "M6_r2_GSE70734_H3K27ac_ZEB1.tsv", sep="\t", index=False)


def gse154675_occupancy():
    try:
        import pyBigWig
    except Exception as e:
        log(f"pyBigWig missing: {e}")
        return
    want = []
    for line in ["ARR", "DU528", "HSB2", "CCRFCEM"]:
        for ab in ["aLMO2", "aTAL1", "aLDB1", "aGATA2"]:
            want.append(f"{line}_{ab}.bigwig")
    dest = PROC / "gse154675_bw"
    # reuse M9 extracts if present
    m9 = M9 / "GSE154675_bw"
    if m9.exists():
        for f in m9.glob("*.bigwig"):
            t = dest / f.name
            if not t.exists():
                dest.mkdir(parents=True, exist_ok=True)
                try:
                    t.symlink_to(f)
                except Exception:
                    t.write_bytes(f.read_bytes())
    tar = DATA / "GSE154675/suppl/GSE154675_RAW.tar"
    files = extract_member(tar, dest, lambda n: n.endswith(".bigwig") and any(n.endswith(w) or w in n for w in want))
    windows = {
        "hg38_gene": ZEB1_HG38, "hg38_prom": ZEB1_PROM38,
        "hg19_gene": ZEB1_HG19, "hg19_prom": ZEB1_PROM19,
    }
    rows = []
    for bw_p in sorted(dest.glob("*.bigwig")):
        if not any(k in bw_p.name for k in ("aLMO2", "aTAL1", "aLDB1", "aGATA2")):
            continue
        if not any(k in bw_p.name for k in ("ARR", "DU528", "HSB2", "CCRFCEM")):
            continue
        try:
            bw = pyBigWig.open(str(bw_p))
        except Exception as e:
            log(f"open fail {bw_p.name}: {e}")
            continue
        rec = {"file": bw_p.name}
        m = re.search(r"(ARR|DU528|HSB2|CCRFCEM)_a(LMO2|TAL1|LDB1|GATA2)", bw_p.name)
        rec["line"] = m.group(1) if m else "unk"
        rec["antibody"] = m.group(2) if m else "unk"
        chroms = set(bw.chroms().keys())
        for tag, (chrom, a, b) in windows.items():
            val = np.nan
            for c in (chrom, chrom.replace("chr", "")):
                if c in chroms:
                    try:
                        end = min(b, int(bw.chroms()[c]))
                        val = bw.stats(c, a, end, type="mean")[0]
                    except Exception:
                        val = np.nan
                    break
            rec[tag] = val
        bw.close()
        rows.append(rec)
        log(f"occupancy {rec}")
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "M6_r2_GSE154675_occupancy.tsv", sep="\t", index=False)
    if len(out):
        # presence: promoter mean > 0 in preferred build
        pref = "hg38_prom" if out["hg38_prom"].notna().any() else "hg19_prom"
        gene = "hg38_gene" if out["hg38_gene"].notna().any() else "hg19_gene"
        out["prom_signal"] = pd.to_numeric(out[pref], errors="coerce")
        out["gene_signal"] = pd.to_numeric(out[gene], errors="coerce")
        out["prom_present"] = (out["prom_signal"].fillna(0) > 0).astype(int)
        wide = out.pivot_table(index="line", columns="antibody", values="prom_present", aggfunc="max")
        wide.to_csv(TAB / "M6_r2_GSE154675_cooccupancy_promoter.tsv", sep="\t")
        sig = out.pivot_table(index="line", columns="antibody", values="prom_signal", aggfunc="mean")
        sig.to_csv(TAB / "M6_r2_GSE154675_promoter_signal.tsv", sep="\t")


def gse280728_scatac():
    tar = DATA / "GSE280728/suppl/GSE280728_RAW.tar"
    dest = PROC / "gse280728"
    beds = extract_member(tar, dest, lambda n: n.endswith("peaks.bed.gz") or n.endswith("-barcodes.tsv.gz"))
    # phenotype from series matrix
    sm = DATA / "GSE280728/matrix/GSE280728_series_matrix.txt.gz"
    meta_rows = []
    if sm.exists():
        with gzip.open(sm, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession") or "characteristics" in line.lower():
                    meta_rows.append(line.strip()[:2000])
        (TAB / "M6_r2_GSE280728_series_head.txt").write_text("\n".join(meta_rows[:40]), encoding="utf-8")
    wins = {"ZEB1_gene_hg38": ZEB1_HG38, "ZEB1_prom_hg38": ZEB1_PROM38,
            "ZEB1_gene_hg19": ZEB1_HG19, "ZEB1_prom_hg19": ZEB1_PROM19}
    rows = []
    for bed in sorted(dest.glob("*peaks.bed.gz")):
        rec = peak_overlap(bed, wins)
        rec["sample"] = re.sub(r"^GSM\d+_", "", bed.name).replace("-peaks.bed.gz", "")
        rows.append(rec)
        log(f"scATAC {rec['sample']} {rec}")
    pd.DataFrame(rows).to_csv(TAB / "M6_r2_GSE280728_ZEB1_peaks.tsv", sep="\t", index=False)
    (TAB / "M6_r2_scATAC_note.txt").write_text(
        "ArchR not installed. GSE280728 provides peak-barcode matrices, not fragments. "
        "Analysis is sample-level ZEB1-overlapping peaks (n=7 libraries). "
        "chromVAR/footprinting/Peak2Gene were not run.\n",
        encoding="utf-8",
    )


def gse287757_mouse_atac():
    tar = DATA / "GSE287757/suppl/GSE287757_RAW.tar"
    dest = PROC / "gse287757"
    files = extract_member(tar, dest, lambda n: n.endswith(".narrowPeak.gz"))
    wins = {"Zeb1_mm10": ZEB1_MM10, "Zeb1_mm39": ZEB1_MM39}
    rows = []
    for f in files:
        rec = peak_overlap(f, wins)
        rec["sample"] = re.sub(r"^GSM\d+_", "", f.name).replace(".narrowPeak.gz", "")
        rec["stage"] = "unk"
        for st in ["ETP_trans", "ETP", "DN2a", "DN2b", "DN3", "DN4", "DN2"]:
            if st in rec["sample"]:
                rec["stage"] = st
                break
        rec["rep"] = (re.search(r"rep(\d+)", rec["sample"], re.I) or type("x", (), {"group": lambda *_: ""})()).group(1)
        rows.append(rec)
        log(f"287757 {rec['sample']} Zeb1_mm39={rec['Zeb1_mm39']} mm10={rec['Zeb1_mm10']}")
    pd.DataFrame(rows).to_csv(TAB / "M6_r2_GSE287757_Zeb1_peaks.tsv", sep="\t", index=False)


def main():
    for label, fn in [
        ("GSE70734", gse70734),
        ("GSE154675", gse154675_occupancy),
        ("GSE280728", gse280728_scatac),
        ("GSE287757", gse287757_mouse_atac),
    ]:
        try:
            fn()
            log(f"OK {label}")
        except Exception as e:
            log(f"FAIL {label}: {type(e).__name__}: {e}")
    log("M6 R2 PREPARE DONE")


if __name__ == "__main__":
    main()
