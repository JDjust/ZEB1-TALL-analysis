#!/usr/bin/env python3
"""GSE69954 ZEB1 promoter 450k betas. CIMP+/- is on-array phenotype (n=65).
No matched RNA in this series; CIMP x MRD is not in GEO metadata.
ZEB1 is minus-strand; TSS200/TSS1500/5'UTR come from HM450 hg19 annotation,
not the TES-adjacent window used earlier as 'promoter'.
"""
from __future__ import annotations

import gzip
import io
import math
import shutil
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

LOCAL_DATA = Path(r"d:\_bioinformation\local_data_staging\revision_tier1\GSE69954")
HULU_DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE69954")
DATA = LOCAL_DATA if (LOCAL_DATA / "matrix" / "GSE69954_series_matrix.txt.gz").exists() else HULU_DATA

if (Path(r"d:\_bioinformation\modules\module6")).exists():
    ROOT = Path(r"d:\_bioinformation\modules\module6")
else:
    ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6")

TAB = ROOT / "tables"
FIG = ROOT / "figures"
ANN = ROOT / "processed"
for d in (TAB, FIG, ANN, ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

# minus-strand TSS ~ hg19 chr10:31818742; keep a genomic window for DMR-style mean
ZEB1_TSS_HG19 = 31818742
ZEB1_PROM_HG19 = ("chr10", ZEB1_TSS_HG19 - 2000, ZEB1_TSS_HG19 + 2000)
ZEB1_GENE_HG19 = ("chr10", 31608101, 31818742)

# fallback TSS/UTR probes if annotation download fails (Illumina HM450 UCSC ZEB1)
FALLBACK_PROM = [
    "cg01472637", "cg02764521", "cg03707168", "cg08856250", "cg16209936",
    "cg17152776", "cg18081940", "cg23210950", "cg24678886", "cg25626449",
    "cg00188947", "cg06278081", "cg08362783", "cg10565456", "cg14048100",
    "cg17792490", "cg23242182", "cg24732722", "cg27533199",
]

ANN_URLS = [
    "https://zwdzwd.github.io/InfiniumAnnotation/current/HM450/HM450.hg19.manifest.tsv.gz",
    "https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg19.manifest.tsv.gz",
]


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
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[np.flatnonzero(mask)[order]] = np.clip(q, 0, 1)
    return out


def mannwhitney(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    try:
        from scipy.stats import mannwhitneyu
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        # exact U via rank; two-sided normal approximation
        n1, n2 = len(a), len(b)
        ranks = pd.Series(np.concatenate([a, b])).rank().to_numpy()
        u = ranks[:n1].sum() - n1 * (n1 + 1) / 2
        mu = n1 * n2 / 2
        sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        z = abs(u - mu) / sigma if sigma else 0
        return float(math.erfc(z / math.sqrt(2)))


def parse_series_meta(path: Path) -> pd.DataFrame:
    title = acc = cimp = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip().strip('"') for v in vals]
            if key == "!Sample_title":
                title = vals
            elif key == "!Sample_geo_accession":
                acc = vals
            elif key == "!Sample_characteristics_ch1" and vals and vals[0].lower().startswith("cimp"):
                cimp = vals
    if not acc:
        raise RuntimeError("no sample accessions")
    d = pd.DataFrame({"geo": acc, "title": title or acc, "cimp_raw": cimp or [""] * len(acc)})
    d["sample_id"] = d["title"].str.replace("T-ALL Sample ", "", regex=False)
    d["cimp"] = np.where(d["cimp_raw"].str.contains("CIMP\\+", regex=True), "CIMP+",
                         np.where(d["cimp_raw"].str.contains("CIMP-", regex=False), "CIMP-", "NA"))
    log(f"samples {len(d)} CIMP+ {(d.cimp=='CIMP+').sum()} CIMP- {(d.cimp=='CIMP-').sum()}")
    return d


def load_gpl_csv(path: Path) -> pd.DataFrame:
    import csv
    rows = []
    header = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        for line in fh:
            if header is None:
                if line.startswith("IlmnID"):
                    header = next(csv.reader([line]))
                continue
            if "ZEB1" not in line:
                continue
            rec = next(csv.reader([line]))
            d = dict(zip(header, rec))
            gene = d.get("UCSC_RefGene_Name", "")
            if "ZEB1" not in gene.split(";"):
                continue
            rows.append({
                "probe": d.get("IlmnID") or d.get("Name"),
                "chr": "chr" + str(d.get("CHR", "")).replace("chr", ""),
                "pos": d.get("MAPINFO"),
                "gene": gene,
                "group": d.get("UCSC_RefGene_Group", ""),
                "cgi": d.get("Relation_to_UCSC_CpG_Island", ""),
            })
    return pd.DataFrame(rows).drop_duplicates("probe")


def load_annotation(want: set[str] | None = None) -> pd.DataFrame:
    cache = ANN / "HM450_ZEB1_probes.tsv"
    gpl = ANN / "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"
    if cache.exists() and cache.stat().st_size > 500:
        d = pd.read_csv(cache, sep="\t")
        if len(d) >= 10:
            log(f"load cached annotation n={len(d)}")
            return d
    if gpl.exists() and gpl.stat().st_size > 1_000_000:
        log(f"parse local GPL {gpl.name}")
        d = load_gpl_csv(gpl)
        if len(d):
            d.to_csv(cache, sep="\t", index=False)
            log(f"ZEB1-locus probes {len(d)}")
            return d
    if cache.exists() and cache.stat().st_size > 200:
        d = pd.read_csv(cache, sep="\t")
        log(f"load cached annotation n={len(d)}")
        return d
    rows = []
    for url in ANN_URLS:
        log(f"download annotation {url}")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ZEB1-TALL-analysis"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                raw = gzip.GzipFile(fileobj=resp)
                wrap = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                header = wrap.readline().rstrip("\n").split("\t")
                idx = {n: i for i, n in enumerate(header)}
                # sesame/zhou columns vary: Probe_ID / CpG_chrm / CpG_beg / genesUniq / UCSC_RefGene_Name
                id_col = next((c for c in ("Probe_ID", "probeID", "IlmnID", "Name") if c in idx), header[0])
                chr_col = next((c for c in ("CpG_chrm", "chr", "CHR", "CpG_chr") if c in idx), None)
                pos_col = next((c for c in ("CpG_beg", "MAPINFO", "pos", "CpG_pos") if c in idx), None)
                gene_col = next((c for c in ("genesUniq", "UCSC_RefGene_Name", "gene", "Gene") if c in idx), None)
                grp_col = next((c for c in ("UCSC_RefGene_Group", "featureUniq", "group") if c in idx), None)
                for line in wrap:
                    p = line.rstrip("\n").split("\t")
                    pid = p[idx[id_col]] if idx[id_col] < len(p) else ""
                    gene = p[idx[gene_col]] if gene_col and idx[gene_col] < len(p) else ""
                    chrom = p[idx[chr_col]] if chr_col and idx[chr_col] < len(p) else ""
                    try:
                        pos = int(float(p[idx[pos_col]])) if pos_col and idx[pos_col] < len(p) and p[idx[pos_col]] else None
                    except ValueError:
                        pos = None
                    grp = p[idx[grp_col]] if grp_col and idx[grp_col] < len(p) else ""
                    hit_gene = "ZEB1" in gene.split(";") or gene == "ZEB1" or "ZEB1" in gene
                    in_win = chrom.replace("chr", "") in {"10"} and pos is not None and ZEB1_GENE_HG19[1] - 3000 <= pos <= ZEB1_GENE_HG19[2] + 3000
                    if hit_gene or in_win:
                        rows.append({"probe": pid, "chr": chrom if str(chrom).startswith("chr") else f"chr{chrom}",
                                     "pos": pos, "gene": gene, "group": grp})
            break
        except Exception as e:
            log(f"annotation fail {e}")
            rows = []
    if not rows:
        log("annotation download failed; using fallback promoter probe list")
        d = pd.DataFrame({"probe": FALLBACK_PROM, "chr": "chr10", "pos": np.nan, "gene": "ZEB1", "group": "TSS1500"})
    else:
        d = pd.DataFrame(rows).drop_duplicates("probe")
    d.to_csv(cache, sep="\t", index=False)
    log(f"ZEB1-locus probes {len(d)}")
    return d


def classify_region(row) -> str:
    """Illumina UCSC_RefGene_Group is the source of truth. A probe can map to
    several ZEB1 isoforms; promoter wins if any isoform is TSS/UTR/exon1.
    3'UTR is gene body, not promoter.
    """
    g = str(row.get("group", "") or "")
    tokens = [t.strip() for t in g.replace(",", ";").split(";") if t.strip()]
    prom_tok = {"TSS200", "TSS1500", "5'UTR", "5UTR", "1stExon"}
    body_tok = {"Body", "3'UTR", "3UTR"}
    if any(t in prom_tok for t in tokens):
        return "promoter"
    if any(t in body_tok for t in tokens):
        return "gene_body"
    pos = row.get("pos")
    if pd.notna(pos):
        x = float(pos)
        # UCSC/Illumina TSS cluster ~31607-31611 kb (hg19)
        if 31605000 <= x <= 31612000:
            return "promoter"
        if ZEB1_GENE_HG19[1] <= x <= ZEB1_GENE_HG19[2]:
            return "gene_body"
    return "other"


def extract_betas(path: Path, probes: set[str]) -> pd.DataFrame:
    """Stream GEO series matrix; keep only requested probes."""
    recs = []
    gsm = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                hdr = next(fh)
                gsm = [x.strip().strip('"') for x in hdr.rstrip("\n").split("\t")[1:]]
                break
        if gsm is None:
            raise RuntimeError("no matrix table")
        for line in fh:
            if line.startswith("!series_matrix_table_end"):
                break
            p = line.rstrip("\n").split("\t")
            pid = p[0].strip().strip('"')
            if pid not in probes:
                continue
            vals = []
            for x in p[1:]:
                x = x.strip().strip('"')
                try:
                    vals.append(float(x) if x not in ("", "NA", "null") else np.nan)
                except ValueError:
                    vals.append(np.nan)
            recs.append((pid, vals))
    if not recs:
        raise RuntimeError("no ZEB1 probes in series matrix")
    mat = pd.DataFrame({pid: vals for pid, vals in recs}, index=gsm).T
    log(f"extracted probes {mat.shape[0]} samples {mat.shape[1]}")
    return mat


def hedges_g(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    na, nb = len(a), len(b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    sp = math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)) if na + nb > 2 else np.nan
    if not sp or not np.isfinite(sp) or sp == 0:
        return np.nan
    g = (a.mean() - b.mean()) / sp
    J = 1 - 3 / (4 * (na + nb) - 9)
    return float(g * J)


def main():
    sm = DATA / "matrix" / "GSE69954_series_matrix.txt.gz"
    if not sm.exists():
        raise SystemExit(f"missing {sm}")
    meta = parse_series_meta(sm)
    meta.to_csv(TAB / "M6_6.17_GSE69954_sample_meta.tsv", sep="\t", index=False)

    ann = load_annotation()
    ann["region"] = ann.apply(classify_region, axis=1)
    ann.to_csv(TAB / "M6_6.17_ZEB1_probe_annotation.tsv", sep="\t", index=False)
    probes = set(ann["probe"].astype(str))
    mat = extract_betas(sm, probes)
    keep = [p for p in mat.index if p in probes]
    mat = mat.loc[keep]
    long = mat.T.reset_index().rename(columns={"index": "geo"})
    long = long.melt(id_vars="geo", var_name="probe", value_name="beta")
    long = long.merge(meta[["geo", "sample_id", "cimp"]], on="geo", how="left")
    long = long.merge(ann[["probe", "chr", "pos", "group", "region"]], on="probe", how="left")
    long["region"] = long["region"].fillna("promoter")
    long.to_csv(TAB / "M6_6.17_ZEB1_promoter_beta.tsv", sep="\t", index=False)

    # per-probe CIMP+ vs CIMP-
    rows = []
    for probe, g in long.groupby("probe"):
        a = g.loc[g.cimp == "CIMP+", "beta"]
        b = g.loc[g.cimp == "CIMP-", "beta"]
        pv = mannwhitney(a, b)
        rec = g.iloc[0]
        rows.append({
            "probe": probe, "chr": rec.get("chr"), "pos": rec.get("pos"),
            "group": rec.get("group"), "region": rec.get("region"),
            "n_CIMP_pos": int(a.notna().sum()), "n_CIMP_neg": int(b.notna().sum()),
            "mean_CIMP_pos": float(a.mean()) if len(a) else np.nan,
            "mean_CIMP_neg": float(b.mean()) if len(b) else np.nan,
            "mean_all": float(g["beta"].mean()),
            "median_all": float(g["beta"].median()),
            "delta_pos_minus_neg": float(a.mean() - b.mean()) if len(a) and len(b) else np.nan,
            "hedges_g": hedges_g(a, b),
            "p_mannwhitney": pv,
        })
    stats = pd.DataFrame(rows)
    stats["q_bh"] = fdr_bh(stats["p_mannwhitney"].values)
    stats = stats.sort_values(["region", "pos", "probe"])
    stats.to_csv(TAB / "M6_6.17_ZEB1_probe_CIMP_stats.tsv", sep="\t", index=False)

    # sample-level mean beta in promoter / gene body
    samp = (long[long.region.isin(["promoter", "gene_body"])]
            .groupby(["geo", "sample_id", "cimp", "region"], dropna=False)["beta"]
            .mean().reset_index(name="mean_beta"))
    samp.to_csv(TAB / "M6_6.17_ZEB1_sample_region_mean.tsv", sep="\t", index=False)

    prom = samp[samp.region == "promoter"]
    a = prom.loc[prom.cimp == "CIMP+", "mean_beta"]
    b = prom.loc[prom.cimp == "CIMP-", "mean_beta"]
    summary = pd.DataFrame([{
        "n_total": int(meta.shape[0]),
        "n_CIMP_pos": int((meta.cimp == "CIMP+").sum()),
        "n_CIMP_neg": int((meta.cimp == "CIMP-").sum()),
        "n_promoter_probes": int((ann.region == "promoter").sum()),
        "n_genebody_probes": int((ann.region == "gene_body").sum()),
        "promoter_mean_CIMP_pos": float(a.mean()) if len(a) else np.nan,
        "promoter_mean_CIMP_neg": float(b.mean()) if len(b) else np.nan,
        "promoter_delta": float(a.mean() - b.mean()) if len(a) and len(b) else np.nan,
        "promoter_hedges_g": hedges_g(a, b),
        "promoter_p": mannwhitney(a, b),
        "n_probe_FDR05": int((stats.q_bh < 0.05).sum()),
        "note": "GSE69954 NOPHO T-ALL 450k, CIMP+/- in GEO; no RNA and no MRD in GEO. Promoter = Illumina TSS200/TSS1500/5UTR/1stExon.",
    }])
    summary.to_csv(TAB / "M6_6.17_ZEB1_promoter_beta_summary.tsv", sep="\t", index=False)
    log(summary.to_string(index=False))
    log("GSE69954 methylation DONE")


if __name__ == "__main__":
    main()
