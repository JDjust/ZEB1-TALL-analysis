"""GSE146901: ETP vs non-ETP ZEB1/ZEB2 expression and merged-loop polarity.

Expression is per sample (8 ETP, 10 non-ETP, 4 normal T cells).
Loops are the author merged calls, one file per group, so chromatin
is a group-level comparison, not a per-patient test.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
from scipy import stats

RAW = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\data\validation\gse146901\raw")
OUT = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\data\validation\gse146901")
OUT.mkdir(parents=True, exist_ok=True)

# Same hg38 spans already used for the HiChIP comparison.
# hg19 spans are Ensembl GRCh37 gene bodies, used only if hg38 hits nothing.
BUILDS = {
    "hg38": {
        "ZEB1": ("chr10", 31318447, 31529814),
        "ZEB2": ("chr2", 144384735, 145277958),
    },
    "hg19": {
        "ZEB1": ("chr10", 31606184, 31819137),
        "ZEB2": ("chr2", 145141942, 145277958),
    },
}

ETP = {"076", "077", "093", "097", "098", "107", "108", "115"}
NON = {"102", "103", "116", "117", "118", "121", "122", "123", "124", "132"}
NORMAL = {"N1", "N2", "N3", "N4"}


def group_of(sample: str) -> str | None:
    key = sample.replace("_RNAseq", "").replace("_HiC", "")
    if key in ETP:
        return "ETP"
    if key in NON:
        return "non-ETP"
    if key in NORMAL:
        return "normal T"
    return None


def load_fpkm() -> pd.DataFrame:
    import gzip
    import io
    path = RAW / "GSE146901_TALL_count_FPKM.xls.gz"
    # Workbook has a count sheet and an FPKM sheet. Counts differ ~9-fold in
    # library size, so the expression comparison uses FPKM.
    raw = gzip.open(path, "rb").read()
    df = pd.read_excel(io.BytesIO(raw), sheet_name="FPKM", engine="xlrd")
    df.columns = [str(c) for c in df.columns]
    print("FPKM columns:", list(df.columns[:8]), "nrow", len(df), flush=True)
    return df


def pick_gene_rows(df: pd.DataFrame) -> pd.DataFrame:
    cols = {c.lower(): c for c in df.columns}
    symbol_col = None
    for key in ("gene_name", "gene symbol", "symbol", "genename", "gene"):
        if key in cols:
            symbol_col = cols[key]
            break
    if symbol_col is None:
        # many GEO xls files put the symbol in the first column
        symbol_col = df.columns[0]
    sub = df[df[symbol_col].astype(str).str.upper().isin(["ZEB1", "ZEB2"])].copy()
    if sub.empty:
        # Ensembl-style "ZEB1" inside a longer id
        mask = df[symbol_col].astype(str).str.contains(r"(?:^|;)ZEB[12](?:;|$)", regex=True)
        sub = df[mask].copy()
    print("gene rows", len(sub), "via", symbol_col, flush=True)
    print(sub.iloc[:, :6].to_string(), flush=True)
    sub["symbol"] = sub[symbol_col].astype(str).str.upper().str.extract(r"(ZEB[12])", expand=False)
    return sub


def expression_table(sub: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [c for c in sub.columns if group_of(str(c)) is not None or group_of(str(c) + "_RNAseq") is not None]
    if not sample_cols:
        raise SystemExit(f"no sample columns among {list(sub.columns)}")
    long = sub.melt(id_vars=["symbol"], value_vars=sample_cols, var_name="sample", value_name="fpkm")
    long["fpkm"] = pd.to_numeric(long["fpkm"], errors="coerce")
    wide = long.pivot_table(index="sample", columns="symbol", values="fpkm", aggfunc="mean")
    wide = wide.reset_index()
    wide["group"] = wide["sample"].map(lambda s: group_of(str(s)))
    # log2(FPKM+1), then z inside the samples that have both genes
    for g in ("ZEB1", "ZEB2"):
        wide[f"log2_{g}"] = (wide[g] + 1).map(lambda x: __import__("math").log2(x))
        mu, sd = wide[f"log2_{g}"].mean(), wide[f"log2_{g}"].std(ddof=1)
        wide[f"z_{g}"] = (wide[f"log2_{g}"] - mu) / sd
    wide["balance"] = wide["z_ZEB1"] - wide["z_ZEB2"]
    wide["log2_ratio"] = wide["log2_ZEB1"] - wide["log2_ZEB2"]
    return wide.sort_values(["group", "sample"])


def wilcox(df: pd.DataFrame, col: str) -> dict:
    a = df.loc[df.group == "ETP", col].dropna()
    b = df.loc[df.group == "non-ETP", col].dropna()
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "metric": col,
        "n_ETP": int(a.size),
        "n_nonETP": int(b.size),
        "median_ETP": float(a.median()),
        "median_nonETP": float(b.median()),
        "median_normal": float(df.loc[df.group == "normal T", col].median()),
        "mannwhitney_p": float(res.pvalue),
    }


def overlap(chrom, start, end, locus) -> bool:
    c, a, b = locus
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom == c and int(start) <= b and int(end) >= a


def count_loops(path: Path, loci: dict) -> dict:
    n = 0
    hits = {k: 0 for k in loci}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.lower().startswith("chrom") or line.startswith("chr1\tx"):
                # keep a real header skip only for non-coordinate headers
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0].replace("chr", "").split("_")[0][:1].isdigit() and not parts[0].startswith("chr"):
                    if n == 0 and not parts[0].startswith("chr") and not parts[0][:1].isdigit():
                        continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                c1, s1, e1, c2, s2, e2 = parts[0], int(parts[1]), int(parts[2]), parts[3], int(parts[4]), int(parts[5])
            except ValueError:
                continue
            n += 1
            for name, locus in loci.items():
                if overlap(c1, s1, e1, locus) or overlap(c2, s2, e2, locus):
                    hits[name] += 1
    return {"n_loops": n, **hits}


def main() -> None:
    sub = pick_gene_rows(load_fpkm())
    expr = expression_table(sub)
    expr.to_csv(OUT / "zeb_expression_by_sample.tsv", sep="\t", index=False)
    tests = pd.DataFrame([wilcox(expr, c) for c in ("log2_ZEB1", "log2_ZEB2", "balance", "log2_ratio")])
    tests.to_csv(OUT / "zeb_expression_etp_vs_nonetp.tsv", sep="\t", index=False)
    print(expr.to_string(index=False), flush=True)
    print(tests.to_string(index=False), flush=True)

    # hg38 gene bodies, plus the same span extended 1 Mb each side.
    # Called loops use bare chromosome names. ZEB2's promoter neighborhood
    # has almost no anchors; the 1 Mb flank is what sits next to the gene.
    loci = {
        "ZEB1_body": ("chr10", 31318447, 31529814),
        "ZEB2_body": ("chr2", 144384735, 145277958),
        "ZEB1_locus_1Mb": ("chr10", 30318447, 32529814),
        "ZEB2_locus_1Mb": ("chr2", 143384735, 146277958),
    }
    expected = {
        "GSE146901_ETP_merged_loops.bedpe.gz": 1306573,
        "GSE146901_non_ETP_merged_loops.bedpe.gz": 1432687,
        "GSE146901_normal_merged_loops.bedpe.gz": 1366914,
    }
    missing = [name for name, size in expected.items() if not (RAW / name).exists() or (RAW / name).stat().st_size != size]
    if missing:
        print("loop files not complete yet:", missing, flush=True)
        return

    rows = []
    for label, fname in (
        ("ETP", "GSE146901_ETP_merged_loops.bedpe.gz"),
        ("non-ETP", "GSE146901_non_ETP_merged_loops.bedpe.gz"),
        ("normal T", "GSE146901_normal_merged_loops.bedpe.gz"),
    ):
        rec = count_loops(RAW / fname, loci)
        n = rec["n_loops"]
        row = {"build": "hg38", "group": label, "n_loops": n}
        for name in loci:
            row[name] = rec[name]
            row[name + "_per_1000"] = 1000 * rec[name] / n
        rows.append(row)
    loops = pd.DataFrame(rows)
    loops.to_csv(OUT / "zeb_merged_loop_rates.tsv", sep="\t", index=False)
    print(loops.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
