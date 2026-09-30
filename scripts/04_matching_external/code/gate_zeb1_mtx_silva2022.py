"""External ZEB1-MTX gate on Canevarolo 2022 (GSE218348).

Phenotype was locked in silva2022_mtx_ic50.tsv before expression was joined.
Expression is the author-deposited GPL570 log2 matrix in the GEO family SOFT
(42450 probesets per sample). Gene value = mean of probesets annotated to
that symbol. Sensitivity: the probeset with the highest mean across all 13
lines, chosen without using IC50.
Primary test: T-ALL only, n=7, Spearman of ZEB1 vs 96 h MTX IC50.
"""
import gzip
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

val = Path(r"d:\_bioinformation\ZEB1\total\data\validation")
probe_map = {}
for line in (val / "hgu133plus2_gate_probes.tsv").read_text().splitlines()[1:]:
    sym, pid = line.split("\t")
    probe_map.setdefault(pid, sym)

samples = []  # list of (title, {probe: float})
current_title = None
in_table = False
header_seen = False
bag = None
with gzip.open(val / "GSE218348_family.soft.gz", "rt", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            current_title = line.split("=", 1)[1].strip()
        elif line.lower().startswith("!sample_table_begin"):
            in_table = True
            header_seen = False
            bag = {}
        elif line.lower().startswith("!sample_table_end"):
            samples.append((current_title, bag))
            in_table = False
        elif in_table:
            if not header_seen:
                header_seen = True
                continue
            pid, value = line.rstrip("\n").split("\t")[:2]
            if pid in probe_map:
                bag[pid] = float(value)

title_to_line = {
    "BCP-ALL cell line 697": "697",
    "BCP-ALL cell line Nalm6": "Nalm6",
    "BCP-ALL cell line REH": "REH",
    "BCP-ALL cell line RS4,11": "RS4;11",
    "BCP-ALL cell line Nalm16": "Nalm16",
    "BCP-ALL cell line Nalm30": "Nalm30",
    "T-ALL cell line ALL-SIL": "ALL-SIL",
    "T-ALL cell line CCRF-CEM": "CCRF-CEM",
    "T-ALL cell line HPB-ALL": "HPB-ALL",
    "T-ALL cell line Jurkat": "Jurkat",
    "T-ALL cell line Molt4": "Molt-4",
    "T-ALL cell line P12-ICHIKAWA": "P12-ICHIKAWA",
    "T-ALL cell line TALL-1": "TALL-1",
}
assert len(samples) == 13
lines = [title_to_line[t] for t, _ in samples]

genes = ["ZEB1", "FPGS", "SLC19A1", "GGH", "DHFR", "TYMS", "SHMT1", "SHMT2"]
by_gene = {g: [] for g in genes}
for pid, g in probe_map.items():
    by_gene[g].append(pid)

present = {g: [] for g in genes}
for g, pids in by_gene.items():
    for pid in pids:
        if all(pid in bag for _, bag in samples):
            present[g].append(pid)

def col(pid):
    return np.array([bag[pid] for _, bag in samples], dtype=float)

gene_mean = {}
gene_best = {}
notes = []
for g in genes:
    pids = present[g]
    if not pids:
        raise SystemExit(f"no probes for {g}")
    mats = np.vstack([col(pid) for pid in pids])
    gene_mean[g] = mats.mean(axis=0)
    means = mats.mean(axis=1)
    order = np.argsort(-means)
    gene_best[g] = mats[order[0]]
    notes.append(
        g + "\t" + ";".join(f"{pids[i]}:{means[i]:.3f}" for i in order)
    )
(val / "silva2022_probe_means.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")

pheno = {}
for rec in (val / "silva2022_mtx_ic50.tsv").read_text().splitlines()[1:]:
    a = rec.split("\t")
    pheno[a[0]] = (a[1], float(a[2]), float(a[3]))
gsh = {}
gsh_header = (val / "silva2022_basal_metabolites.tsv").read_text().splitlines()
gsh_cols = gsh_header[0].split("\t")
for rec in gsh_header[1:]:
    a = rec.split("\t")
    gsh[a[0]] = {c: float(v) for c, v in zip(gsh_cols[1:], a[1:])}

rows_out = ["cell_line\tlineage\tic50_48h_nM\tic50_96h_nM\tGlutathione\t" + "\t".join(genes) + "\t" + "\t".join(g + "_bestprobe" for g in genes)]
for i, line in enumerate(lines):
    lineage, ic48, ic96 = pheno[line]
    vals = [f"{gene_mean[g][i]:.6f}" for g in genes]
    best = [f"{gene_best[g][i]:.6f}" for g in genes]
    rows_out.append("\t".join([
        line, lineage, str(ic48), str(ic96), f"{gsh[line]['Glutathione']:.6f}",
        *vals, *best,
    ]))
(val / "silva2022_zeb1_mtx_joined.tsv").write_text("\n".join(rows_out) + "\n", encoding="utf-8")

def spearman(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    rho, p = spearmanr(x[ok], y[ok])
    return int(ok.sum()), float(rho), float(p)

def partial_spearman(x, y, z):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    rx = np.argsort(np.argsort(x[ok])).astype(float)
    ry = np.argsort(np.argsort(y[ok])).astype(float)
    rz = np.argsort(np.argsort(z[ok])).astype(float)
    xr = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    yr = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rho, p = pearsonr(xr, yr)
    return int(ok.sum()), float(rho), float(p)

records = {line: i for i, line in enumerate(lines)}

def take(cohort_lines, key):
    idx = [records[ln] for ln in cohort_lines]
    if key == "ic50_96h":
        return np.array([pheno[ln][2] for ln in cohort_lines])
    if key == "ic50_48h":
        return np.array([pheno[ln][1] for ln in cohort_lines])
    if key == "Glutathione":
        return np.array([gsh[ln]["Glutathione"] for ln in cohort_lines])
    if key.endswith("_bestprobe"):
        return gene_best[key.replace("_bestprobe", "")][idx]
    return gene_mean[key][idx]

tall = [ln for ln in lines if pheno[ln][0] == "T-ALL"]
bcp = [ln for ln in lines if pheno[ln][0] == "BCP-ALL"]
covars = ["FPGS", "SLC19A1", "GGH", "DHFR", "TYMS", "SHMT1", "SHMT2", "Glutathione"]

gate_rows = ["cohort\ty\tx\tn\trho\tp"]
def add(cohort, yname, xname, stat):
    n, rho, p = stat
    gate_rows.append(f"{cohort}\t{yname}\t{xname}\t{n}\t{rho:.6f}\t{p:.6g}")

for cohort, members in (("T-ALL", tall), ("BCP-ALL", bcp)):
    y = take(members, "ic50_96h")
    add(cohort, "ic50_96h", "ZEB1", spearman(take(members, "ZEB1"), y))
    add(cohort, "ic50_96h", "ZEB1_bestprobe", spearman(take(members, "ZEB1_bestprobe"), y))
    add(cohort, "ic50_48h", "ZEB1", spearman(take(members, "ZEB1"), take(members, "ic50_48h")))
    for g in covars:
        add(cohort, "ic50_96h", g, spearman(take(members, g), y))
        add(cohort, "ZEB1", g, spearman(take(members, "ZEB1"), take(members, g)))
        add(cohort, "ic50_96h_partial", f"ZEB1|{g}", partial_spearman(take(members, "ZEB1"), y, take(members, g)))
    for drop in members:
        keep = [ln for ln in members if ln != drop]
        st = spearman(take(keep, "ZEB1"), take(keep, "ic50_96h"))
        add(f"{cohort}_loo_{drop}", "ic50_96h", "ZEB1", st)

(val / "silva2022_zeb1_mtx_gate.tsv").write_text("\n".join(gate_rows) + "\n", encoding="utf-8")
print("\n".join(notes))
print("--- T-ALL values ---")
for ln in tall:
    i = records[ln]
    print(f"{ln}\tZEB1 {gene_mean['ZEB1'][i]:.3f}\tbest {gene_best['ZEB1'][i]:.3f}\tIC50_96 {pheno[ln][2]}\tGSH {gsh[ln]['Glutathione']:.3f}")
print("--- gate ---")
for row in gate_rows:
    if row.startswith("T-ALL\t") or row.startswith("BCP-ALL\t") or row.startswith("cohort"):
        print(row)
