"""A3-2/A3-3: ZEB loci occupancy and H3K27ac HiChIP contacts.

Primary ChIP question: reproducible BCL11B peak at ZEB1/ZEB2 promoter ±2 kb or gene body.
HiChIP is H3K27ac, not BCL11B HiChIP. DND41 replicates are collapsed.
"""
from pathlib import Path
import gzip
import pandas as pd
import numpy as np

ROOT = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a3_gse165209")
CHIP = ROOT / "gse165016"
HIC = ROOT / "gse165207"
OUT = ROOT / "locus"
OUT.mkdir(parents=True, exist_ok=True)

# GRCh38 / Ensembl
LOCI = {
    "ZEB1": dict(chrom="chr10", start=31318495, end=31529814, strand="+", tss=31318495),
    "ZEB2": dict(chrom="chr2", start=144364364, end=144521057, strand="-", tss=144521057),
}
PROM_FLANK = 2000
CONTACT_FLANK = 5000

BCL11B_PEAKS = {
    "DND41": CHIP / "GSM5024537_DND41_BCL11B_peaks.narrowPeak.gz",
    "SJTALL005006_rep1": CHIP / "GSM5024540_SJTALL005006_BCL11B_rep1_peaks.narrowPeak.gz",
    "SJTALL005006_rep2": CHIP / "GSM5024541_SJTALL005006_BCL11B_rep2_peaks.narrowPeak.gz",
    "SJALL068666": CHIP / "GSM5265324_BCL11B_E15804_peaks.narrowPeak.gz",
    "SJALL067672": CHIP / "GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz",
    "SJALL067671": CHIP / "GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz",
}
H3K27AC_PEAKS = {
    "DND41": CHIP / "GSM5024536_DND41_H3K27ac_peaks.narrowPeak.gz",
    "SJTALL005006": CHIP / "GSM5024539_SJTALL005006_H3K27ac_peaks.narrowPeak.gz",
}
# 5 primary BCL11B-related + collapsed DND41 + Jurkat + CD34
BEDPE = {
    "SJALL068279": dict(path=HIC / "GSM5028224_H3K27ac_HiChIP_loops_SJALL068279.bedpe.gz", group="BCL11B_primary"),
    "SJAUL068292": dict(path=HIC / "GSM5028225_H3K27ac_HiChIP_loops_SJAUL068292.bedpe.gz", group="BCL11B_primary"),
    "SJMPAL011914": dict(path=HIC / "GSM5028226_H3K27ac_HiChIP_loops_SJMPAL011914.bedpe.gz", group="BCL11B_primary"),
    "SJTALL005006": dict(path=HIC / "GSM5028227_H3K27ac_HiChIP_loops_SJTALL005006.bedpe.gz", group="BCL11B_primary"),
    "SJMPAL011911": dict(path=HIC / "GSM5028228_H3K27ac_HiChIP_loops_SJMPAL011911.bedpe.gz", group="BCL11B_primary"),
    "DND41_rep1": dict(path=HIC / "GSM5028229_H3K27ac_HiChIP_loops_DND41_rep1.bedpe.gz", group="DND41"),
    "DND41_rep2": dict(path=HIC / "GSM5028230_H3K27ac_HiChIP_loops_DND41_rep2.bedpe.gz", group="DND41"),
    "Jurkat": dict(path=HIC / "GSM5028231_H3K27ac_HiChIP_loops_Jurkat.bedpe.gz", group="Tlineage_other"),
    "CD34": dict(path=HIC / "GSM5028232_H3K27ac_HiChIP_loops_CD34_5M.bedpe.gz", group="CD34"),
}


def norm_chrom(x):
    x = str(x)
    return x if x.startswith("chr") else "chr" + x


def read_narrowpeak(path):
    cols = ["chrom", "start", "end", "name", "score", "strand", "signal", "pval", "qval", "peak"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, compression="gzip")
    df["chrom"] = df["chrom"].map(norm_chrom)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def overlaps(a0, a1, b0, b1):
    return (a0 < b1) and (a1 > b0)


def locus_windows(gene):
    loc = LOCI[gene]
    prom0 = loc["tss"] - PROM_FLANK
    prom1 = loc["tss"] + PROM_FLANK
    return loc, (prom0, prom1), (loc["start"], loc["end"])


def nearest_and_overlap(peaks, gene):
    loc, prom, body = locus_windows(gene)
    sub = peaks.loc[peaks.chrom == loc["chrom"]].copy()
    if sub.empty:
        return dict(n_promoter=0, n_body=0, nearest_bp=np.nan, nearest_signal=np.nan)
    mid = (sub.start + sub.end) // 2
    dist = np.minimum(np.abs(mid - loc["tss"]), np.minimum(np.abs(sub.start - loc["tss"]), np.abs(sub.end - loc["tss"])))
    # distance 0 if overlaps gene body
    body_hit = (sub.start < body[1]) & (sub.end > body[0])
    dist = np.where(body_hit, 0, dist)
    i = int(np.argmin(dist))
    n_prom = int(((sub.start < prom[1]) & (sub.end > prom[0])).sum())
    n_body = int(body_hit.sum())
    return dict(n_promoter=n_prom, n_body=n_body, nearest_bp=int(dist.iloc[i] if hasattr(dist, "iloc") else dist[i]),
                nearest_signal=float(sub.iloc[i]["signal"]))


def read_bedpe(path):
    raw = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    raw = raw.rename(columns={0: "c1", 1: "s1", 2: "e1", 3: "c2", 4: "s2", 5: "e2"})
    raw["c1"] = raw["c1"].map(norm_chrom)
    raw["c2"] = raw["c2"].map(norm_chrom)
    for col in ("s1", "e1", "s2", "e2"):
        raw[col] = raw[col].astype(int)
    # count column: last numeric column if present
    num = raw.select_dtypes(include=[np.number])
    raw["count"] = num.iloc[:, -1] if num.shape[1] else 1.0
    return raw


def promoter_centered(bedpe, gene):
    loc = LOCI[gene]
    w0, w1 = loc["tss"] - CONTACT_FLANK, loc["tss"] + CONTACT_FLANK
    a = (bedpe.c1 == loc["chrom"]) & (bedpe.s1 < w1) & (bedpe.e1 > w0)
    b = (bedpe.c2 == loc["chrom"]) & (bedpe.s2 < w1) & (bedpe.e2 > w0)
    hit = bedpe.loc[a | b]
    total = float(bedpe["count"].sum())
    raw = float(hit["count"].sum()) if len(hit) else 0.0
    cpm = 1e6 * raw / total if total > 0 else np.nan
    return raw, cpm, int(len(hit)), total


def main():
    chip_rows = []
    for sample, path in BCL11B_PEAKS.items():
        if not path.exists() or path.stat().st_size < 50000:
            print("SKIP missing/bad", path.name, flush=True)
            continue
        peaks = read_narrowpeak(path)
        for gene in LOCI:
            d = nearest_and_overlap(peaks, gene)
            chip_rows.append(dict(assay="BCL11B_ChIP", sample_id=sample, gene=gene, n_peaks_genome=len(peaks), **d))
    for sample, path in H3K27AC_PEAKS.items():
        peaks = read_narrowpeak(path)
        for gene in LOCI:
            d = nearest_and_overlap(peaks, gene)
            chip_rows.append(dict(assay="H3K27ac", sample_id=sample, gene=gene, n_peaks_genome=len(peaks), **d))
    chip = pd.DataFrame(chip_rows)
    chip["promoter_hit"] = chip["n_promoter"] > 0
    chip["body_hit"] = chip["n_body"] > 0
    chip.to_csv(OUT / "a3_chip_zeb_loci.tsv", sep="\t", index=False)

    # reproducibility: exclude DND41 from the "independent primary" count
    bcl = chip.loc[chip.assay == "BCL11B_ChIP"].copy()
    bcl["independent"] = ~bcl.sample_id.str.startswith("DND41")
    # collapse SJTALL005006 reps: hit if either rep hits
    rec = []
    for gene in LOCI:
        g = bcl.loc[bcl.gene == gene]
        prim = g.loc[g.independent]
        # unique cases among files that actually loaded
        cases = {
            "SJTALL005006": bool(prim.loc[prim.sample_id.str.startswith("SJTALL005006"), "promoter_hit"].any()) if prim.sample_id.str.startswith("SJTALL005006").any() else None,
            "SJALL068666": bool(prim.loc[prim.sample_id.eq("SJALL068666"), "promoter_hit"].any()) if prim.sample_id.eq("SJALL068666").any() else None,
            "SJALL067672": bool(prim.loc[prim.sample_id.eq("SJALL067672"), "promoter_hit"].any()) if prim.sample_id.eq("SJALL067672").any() else None,
            "SJALL067671": bool(prim.loc[prim.sample_id.eq("SJALL067671"), "promoter_hit"].any()) if prim.sample_id.eq("SJALL067671").any() else None,
        }
        body_cases = {
            "SJTALL005006": bool(prim.loc[prim.sample_id.str.startswith("SJTALL005006"), "body_hit"].any()) if prim.sample_id.str.startswith("SJTALL005006").any() else None,
            "SJALL068666": bool(prim.loc[prim.sample_id.eq("SJALL068666"), "body_hit"].any()) if prim.sample_id.eq("SJALL068666").any() else None,
            "SJALL067672": bool(prim.loc[prim.sample_id.eq("SJALL067672"), "body_hit"].any()) if prim.sample_id.eq("SJALL067672").any() else None,
            "SJALL067671": bool(prim.loc[prim.sample_id.eq("SJALL067671"), "body_hit"].any()) if prim.sample_id.eq("SJALL067671").any() else None,
        }
        rec.append(dict(
            gene=gene,
            n_independent_cases=int(sum(v is not None for v in cases.values())),
            n_promoter_independent=int(sum(bool(v) for v in cases.values() if v is not None)),
            n_body_independent=int(sum(bool(v) for v in body_cases.values() if v is not None)),
            DND41_promoter=bool(g.loc[g.sample_id.eq("DND41"), "promoter_hit"].any()),
            DND41_body=bool(g.loc[g.sample_id.eq("DND41"), "body_hit"].any()),
            promoter_cases=",".join([k for k, v in cases.items() if v]) or "none",
            body_cases=",".join([k for k, v in body_cases.items() if v]) or "none",
        ))
    rec = pd.DataFrame(rec)
    rec.to_csv(OUT / "a3_chip_reproducibility.tsv", sep="\t", index=False)

    hic_rows = []
    for sample, meta in BEDPE.items():
        if not meta["path"].exists() or meta["path"].stat().st_size < 50000:
            print("SKIP missing/bad", meta["path"].name, flush=True)
            continue
        bedpe = read_bedpe(meta["path"])
        for gene in LOCI:
            raw, cpm, nloop, total = promoter_centered(bedpe, gene)
            hic_rows.append(dict(sample=sample, group=meta["group"], gene=gene,
                                 n_loops=nloop, raw_count=raw, cpm=cpm, total_count=total))
    hic = pd.DataFrame(hic_rows)
    hic.to_csv(OUT / "a3_hichip_zeb_promoter_contacts.tsv", sep="\t", index=False)
    if hic.empty:
        print("No valid BEDPE yet; ChIP-only Gate C this round.")
        rec.to_csv(OUT / "a3_chip_reproducibility.tsv", sep="\t", index=False)
        print(chip.to_string(index=False))
        print(rec.to_string(index=False))
        return

    # collapse DND41 reps
    wide = []
    for sample in hic.sample.unique():
        z1 = hic.loc[(hic["sample"] == sample) & (hic.gene == "ZEB1"), "cpm"].iloc[0]
        z2 = hic.loc[(hic["sample"] == sample) & (hic.gene == "ZEB2"), "cpm"].iloc[0]
        grp = hic.loc[hic["sample"] == sample, "group"].iloc[0]
        wide.append(dict(sample=sample, group=grp, ZEB1_cpm=z1, ZEB2_cpm=z2,
                         zeb2_over_zeb1=z2 / z1 if z1 else np.nan))
    wide = pd.DataFrame(wide)
    dnd = wide.loc[wide.group == "DND41", ["ZEB1_cpm", "ZEB2_cpm", "zeb2_over_zeb1"]].mean()
    collapsed = pd.concat([
        wide.loc[wide.group != "DND41"],
        pd.DataFrame([dict(sample="DND41_collapsed", group="DND41", **dnd.to_dict())]),
    ], ignore_index=True)
    collapsed.to_csv(OUT / "a3_hichip_sample_units.tsv", sep="\t", index=False)

    prim = collapsed.loc[collapsed.group == "BCL11B_primary"]
    other = collapsed.loc[collapsed.group.isin(["Tlineage_other", "CD34"])]
    def mw(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        if len(a) < 2 or len(b) < 1:
            return np.nan
        try:
            from scipy.stats import mannwhitneyu
            return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:
            return np.nan

    summary = pd.DataFrame([
        dict(item="ZEB1_promoter_independent_BCL11B", value=int(rec.loc[rec.gene == "ZEB1", "n_promoter_independent"].iloc[0])),
        dict(item="ZEB2_promoter_independent_BCL11B", value=int(rec.loc[rec.gene == "ZEB2", "n_promoter_independent"].iloc[0])),
        dict(item="ZEB1_body_independent_BCL11B", value=int(rec.loc[rec.gene == "ZEB1", "n_body_independent"].iloc[0])),
        dict(item="ZEB2_body_independent_BCL11B", value=int(rec.loc[rec.gene == "ZEB2", "n_body_independent"].iloc[0])),
        dict(item="direct_occupancy_support",
             value="YES" if (rec.loc[rec.gene == "ZEB2", "n_promoter_independent"].iloc[0] >= 2) else "NO"),
        dict(item="BCL11B_primary_n", value=int(len(prim))),
        dict(item="BCL11B_primary_ZEB2_cpm_median", value=float(prim.ZEB2_cpm.median())),
        dict(item="BCL11B_primary_ZEB1_cpm_median", value=float(prim.ZEB1_cpm.median())),
        dict(item="BCL11B_primary_ratio_median", value=float(prim.zeb2_over_zeb1.median())),
        dict(item="other_ZEB2_cpm_median", value=float(other.ZEB2_cpm.median()) if len(other) else np.nan),
        dict(item="other_ZEB1_cpm_median", value=float(other.ZEB1_cpm.median()) if len(other) else np.nan),
        dict(item="other_ratio_median", value=float(other.zeb2_over_zeb1.median()) if len(other) else np.nan),
        dict(item="ratio_MW_p_primary_vs_other", value=mw(prim.zeb2_over_zeb1, other.zeb2_over_zeb1)),
        dict(item="ZEB2_cpm_MW_p_primary_vs_other", value=mw(prim.ZEB2_cpm, other.ZEB2_cpm)),
        dict(item="note", value="H3K27ac HiChIP is chromatin consistency, not BCL11B looping. DND41 collapsed and kept separate."),
    ])
    summary.to_csv(OUT / "a3_locus_gateC_summary.tsv", sep="\t", index=False)
    print(chip.to_string(index=False))
    print(rec.to_string(index=False))
    print(collapsed.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
