"""Preregistered replication: ZEB1/ZEB2 balance vs daunorubicin LC50 in GSE2351.

Discovery anchor is St. Jude pharmacotype T-ALL (Lee 2023), n=77,
Spearman rho=-0.281 for log2((ZEB1+eps)/(ZEB2+eps)) vs daunorubicin LC50.

GSE2351 is the Dutch/German pediatric ALL microarray cohort (Lugthart /
Holleman), not Total Therapy XV-XVII. Expression is HG-U133A. The matched
quantity on a log-scale array is the probe difference ZEB1 - ZEB2.
"""
import csv
import gzip
import math
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(r"d:/_bioinformation/ZEB1/total/data/validation")
MATRIX = ROOT / "GSE2351_series_matrix.txt.gz"
CLIN = ROOT / "gse2351_clinical.csv"
OUT = ROOT / "gse2351_zeb_ratio_daunorubicin.tsv"

ZEB1 = ["212764_at", "212758_s_at", "210875_s_at"]
ZEB2 = "203603_s_at"
PROBES = set(ZEB1 + [ZEB2])


def spearman(x, y):
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def main():
    samples = None
    expr = {}
    with gzip.open(MATRIX, "rt", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                next(handle)
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            pid = line.split("\t", 1)[0].strip().strip('"')
            if pid not in PROBES:
                continue
            vals = []
            for token in line.rstrip("\n").split("\t")[1:]:
                token = token.strip().strip('"')
                vals.append(np.nan if token in ("", "null", "NA") else float(token))
            expr[pid] = np.array(vals, dtype=float)

    clin = {}
    with CLIN.open(newline="") as handle:
        for row in csv.DictReader(handle):
            clin[row["gsm"]] = row

    rows = []
    for i, gsm in enumerate(samples):
        meta = clin.get(gsm)
        if meta is None:
            continue
        dnr = meta["LC50 DNR"].strip()
        if dnr == "":
            continue
        rec = {"gsm": gsm, "code": meta["code"], "lc50_dnr": float(dnr)}
        for pid in ZEB1 + [ZEB2]:
            rec[pid] = float(expr[pid][i])
        # Series-matrix values are linear (about 10-1400), so the discovery
        # quantity log2((ZEB1+eps)/(ZEB2+eps)) is a difference of logs.
        eps = 1e-6
        rec["zeb1_mean"] = float(np.mean([rec[p] for p in ZEB1]))
        rec["ratio_primary"] = math.log2((rec["212764_at"] + eps) / (rec[ZEB2] + eps))
        rec["ratio_mean"] = math.log2((rec["zeb1_mean"] + eps) / (rec[ZEB2] + eps))
        rows.append(rec)

    lc50 = np.array([r["lc50_dnr"] for r in rows])
    print(f"n_joined={len(rows)}")
    print(f"expression_range_212764 {np.nanmin(expr['212764_at']):.2f} {np.nanmax(expr['212764_at']):.2f}")
    print(f"lc50_dnr min {lc50.min()} median {np.median(lc50)} max {lc50.max()}")
    print(f"lc50_at_floor_0.002 {(lc50 <= 0.002).sum()} at_ceiling_ge_2 {(lc50 >= 2).sum()}")

    tests = [
        ("log2_ratio_212764", "ratio_primary"),
        ("log2_ratio_meanZEB1", "ratio_mean"),
        ("ZEB1_212764_at", "212764_at"),
        ("ZEB1_212758_s_at", "212758_s_at"),
        ("ZEB1_210875_s_at", "210875_s_at"),
        ("ZEB1_mean", "zeb1_mean"),
        ("ZEB2_203603_s_at", "203603_s_at"),
    ]
    print("feature\tn\trho\tp")
    for name, key in tests:
        x = np.array([r[key] for r in rows])
        ok = np.isfinite(x) & np.isfinite(lc50)
        rho, p = spearman(x[ok], lc50[ok])
        print(f"{name}\t{ok.sum()}\t{rho:.4f}\t{p:.4g}")

    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
