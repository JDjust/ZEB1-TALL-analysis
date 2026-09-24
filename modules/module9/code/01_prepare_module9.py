#!/usr/bin/env python3
"""Module 9 prepare: TLX1 / LMO2 perturbation counts, patient TLX ranks, ChIP overlap.

Do not impose LMO2→ZEB1 causality. Extract matrices and metadata only.
"""
from __future__ import annotations

import gzip
import io
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
AN = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis")
ROOT = AN / "module9"
PROC = ROOT / "processed"
TAB = ROOT / "tables"
for d in (PROC, TAB, ROOT / "figures", ROOT / "logs"):
    d.mkdir(parents=True, exist_ok=True)

KEY = [
    "ZEB1", "ZEB2", "LMO2", "TLX1", "TLX3", "TAL1", "LYL1", "GATA3", "TCF7",
    "BCL11B", "MEF2C", "RUNX1", "MYC", "IL7R", "CD34", "CD1A", "NOTCH1",
    "SPI1", "RAG1", "DNTT", "LEF1", "SOX4",
]
MOUSE_KEY = {g.capitalize() if g not in {"MYC", "RAG1", "DNTT", "IL7R", "CD34", "CD1A", "LEF1", "SOX4"} else g.title() if g in {"MYC"} else g.capitalize(): g for g in KEY}
MOUSE_KEY.update({
    "Zeb1": "ZEB1", "Zeb2": "ZEB2", "Lmo2": "LMO2", "Tlx1": "TLX1", "Tlx3": "TLX3",
    "Tal1": "TAL1", "Lyl1": "LYL1", "Gata3": "GATA3", "Tcf7": "TCF7", "Bcl11b": "BCL11B",
    "Mef2c": "MEF2C", "Runx1": "RUNX1", "Myc": "MYC", "Il7r": "IL7R", "Cd34": "CD34",
    "Cd1a": "CD1A", "Notch1": "NOTCH1", "Spi1": "SPI1", "Rag1": "RAG1", "Dntt": "DNTT",
    "Lef1": "LEF1", "Sox4": "SOX4",
})
# mm10 / mm39 Zeb1 windows (generous)
ZEB1_MM = {"chr18": (5_600_000, 5_850_000), "18": (5_600_000, 5_850_000)}
ZEB1_HG38 = {"chr10": (31_318_447, 31_529_814), "10": (31_318_447, 31_529_814)}
ZEB1_HG19 = {"chr10": (31_607_511, 31_818_878), "10": (31_607_511, 31_818_878)}
ZEB1_PROM_HG38 = {"chr10": (31_316_447, 31_328_447), "10": (31_316_447, 31_328_447)}
ZEB1_PROM_HG19 = {"chr10": (31_605_511, 31_617_511), "10": (31_605_511, 31_617_511)}


def log(m):
    print(m, flush=True)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    titles, gsm, chars = [], [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        for line in f:
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
        df[f"char{i+1}"] = (row + [""] * n)[:n]
    return df


def extract_tar_members(tar_p: Path, dest: Path, pred) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    out = []
    with tarfile.open(tar_p, "r") as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            name = Path(m.name).name
            if not pred(name):
                continue
            target = dest / name
            if not target.exists() or target.stat().st_size < 20:
                src = tf.extractfile(m)
                target.write_bytes(src.read())
            out.append(target)
    return sorted(out)


def read_star_rpg(path: Path) -> pd.Series:
    genes, vals = [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("N_") or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            gid = p[0].split(".")[0]
            try:
                v = float(p[1])
            except ValueError:
                continue
            genes.append(gid)
            vals.append(v)
    return pd.Series(vals, index=pd.Index(genes, name="gene_id"), dtype="float64")


def build_ensembl_map() -> pd.DataFrame:
    star_dir = DATA / "TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts"
    files = sorted(star_dir.glob("*.tsv"))
    if not files:
        log("WARNING no TARGET STAR files for gene map")
        return pd.DataFrame(columns=["gene_id", "symbol", "gene_type"])
    rows = []
    with files[0].open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or line.startswith("gene_id") or line.startswith("N_"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            rows.append((p[0].split(".")[0], p[1], p[2]))
    mp = pd.DataFrame(rows, columns=["gene_id", "symbol", "gene_type"]).drop_duplicates("gene_id")
    mp.to_csv(PROC / "ensembl_to_symbol.tsv", sep="\t", index=False)
    log(f"ensembl map {len(mp)} from {files[0].name}")
    return mp


def collapse_symbol(counts: pd.DataFrame, mp: pd.DataFrame) -> pd.DataFrame:
    """Map Ensembl rows to symbols; keep max-count gene if many-to-one."""
    idx = counts.index.to_series().astype(str).str.split(".").str[0]
    counts = counts.copy()
    counts.index = idx
    m = mp.drop_duplicates("gene_id").set_index("gene_id")["symbol"]
    sym = idx.map(m)
    keep = sym.notna() & (sym.astype(str) != "") & ~sym.astype(str).str.startswith("ENS")
    sub = counts.loc[keep.values]
    sub.index = pd.Index(sym[keep].astype(str), name="gene")
    # sum split Ensembl of the same symbol (rare)
    sub = sub.groupby(level=0).sum()
    return sub


def write_counts(name: str, counts: pd.DataFrame, coldata: pd.DataFrame):
    counts.to_csv(PROC / f"{name}_counts.tsv.gz", sep="\t")
    coldata.to_csv(PROC / f"{name}_coldata.tsv", sep="\t", index=False)
    log(f"wrote {name} genes={counts.shape[0]} samples={counts.shape[1]}")
    log(str(coldata))


def prepare_star_kd(acc: str, classify_title):
    meta = parse_series_matrix(DATA / acc / "matrix" / f"{acc}_series_matrix.txt.gz")
    files = extract_tar_members(
        DATA / acc / "suppl" / f"{acc}_RAW.tar",
        PROC / acc,
        lambda n: n.endswith("ReadsPerGene.out.tab.gz"),
    )
    cols = {}
    rows = []
    for p in files:
        m = re.search(r"(GSM\d+)", p.name)
        gsm = m.group(1) if m else p.name
        hit = meta[meta["gsm"].eq(gsm)]
        title = hit["title"].iloc[0] if len(hit) else p.name
        group, sirna, rep = classify_title(title)
        cols[gsm] = read_star_rpg(p)
        rows.append({"gsm": gsm, "title": title, "group": group, "sirna": sirna, "rep": rep, "file": p.name})
    counts = pd.DataFrame(cols).fillna(0).astype(np.int64)
    coldata = pd.DataFrame(rows)
    return counts, coldata, meta


def class_110635(title: str):
    t = title.lower()
    rep = re.search(r"repl\s*(\d)", t)
    r = int(rep.group(1)) if rep else np.nan
    if "scrambled" in t or "ntc" in t:
        return "control", "NTC", r
    if "s46" in t:
        return "siTLX1_1", "s46", r
    if "s47" in t:
        return "siTLX1_2", "s47", r
    return "other", "unk", r


def class_110632(title: str):
    t = title.lower()
    rep = re.search(r"rep(\d)", t)
    r = int(rep.group(1)) if rep else np.nan
    if "ntc" in t:
        return "control", "NTC", r
    if "sirna1" in t:
        return "siTLX1_1", "siRNA1", r
    if "sirna2" in t:
        return "siTLX1_2", "siRNA2", r
    return "other", "unk", r


def parse_geo_matrix_table(path: Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    in_table = False
    header = None
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = p
                continue
            if p and p[0]:
                rows.append(p)
    if not header or not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows, columns=header)
    df = df.set_index(df.columns[0])
    return df.apply(pd.to_numeric, errors="coerce")


def parse_agilent_fe(path: Path) -> pd.Series:
    opener = gzip.open if str(path).endswith(".gz") else open
    header = None
    genes, vals = [], []
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if header is None:
                if line.startswith("FEATURES"):
                    header = line.rstrip("\n").split("\t")
                continue
            if not line.startswith("DATA"):
                continue
            p = line.rstrip("\n").split("\t")
            rec = {header[i]: p[i] if i < len(p) else "" for i in range(len(header))}
            probe = rec.get("ProbeName") or rec.get("SystematicName") or rec.get("GeneName")
            sig = rec.get("gProcessedSignal") or rec.get("gMeanSignal") or rec.get("gMedianSignal")
            if not probe or not sig:
                continue
            try:
                v = float(sig)
            except ValueError:
                continue
            genes.append(probe)
            vals.append(v)
    if not genes:
        raise ValueError(f"no FEATURES/gProcessedSignal in {path.name}")
    return pd.Series(vals, index=pd.Index(genes, name="probe")).groupby(level=0).mean()


def prepare_gse6214x():
    for acc in ["GSE62143", "GSE62141"]:
        tar = DATA / acc / "suppl" / f"{acc}_RAW.tar"
        extract_tar_members(
            tar, PROC / acc,
            lambda n: n.endswith(".txt.gz") and (n.startswith("GSM") or n.startswith("GPL")),
        )
        matrix = DATA / acc / "matrix" / f"{acc}_series_matrix.txt.gz"
        if acc == "GSE62141" and not matrix.exists():
            meta = parse_series_matrix(DATA / acc / "soft" / "GSE62141_family.soft.gz") if (DATA / acc / "soft" / "GSE62141_family.soft.gz").exists() else pd.DataFrame()
            expr = pd.DataFrame()
        else:
            meta = parse_series_matrix(matrix) if matrix.exists() else pd.DataFrame()
            expr = parse_geo_matrix_table(matrix) if matrix.exists() else pd.DataFrame()
        log(f"{acc} series table {expr.shape}")
        if expr.empty:
            files = sorted((PROC / acc).glob("GSM*.txt.gz"))
            cols = {}
            for p in files:
                gsm = re.search(r"(GSM\d+)", p.name).group(1)
                log(f"  Agilent FE {p.name}")
                cols[gsm] = parse_agilent_fe(p)
            expr = pd.DataFrame(cols)
            log(f"{acc} parsed FE {expr.shape}")
        rows = []
        for gsm in expr.columns:
            title = ""
            if len(meta) and "gsm" in meta.columns:
                hit = meta[meta["gsm"].astype(str).eq(str(gsm))]
                if len(hit):
                    title = str(hit["title"].iloc[0])
            files = list((PROC / acc).glob(f"{gsm}*.txt.gz"))
            fname = files[0].name if files else ""
            blob = (title + " " + fname).lower()
            if "empty" in blob:
                group = "empty"
            elif "tlx1_vector" in blob or "tlx1 vector" in blob:
                group = "TLX1_OE"
            elif "scrambled" in blob:
                group = "control"
            elif "sirna1" in blob:
                group = "siTLX1_1"
            elif "sirna2" in blob:
                group = "siTLX1_2"
            else:
                group = "other"
            rows.append({"sample": str(gsm), "gsm": str(gsm), "title": title, "file": fname, "group": group})
        expr.index.name = "gene"
        expr.to_csv(PROC / f"{acc.lower()}_expr.tsv.gz", sep="\t")
        pd.DataFrame(rows).to_csv(PROC / f"{acc.lower()}_coldata.tsv", sep="\t", index=False)
        log(f"{acc} expr {expr.shape} groups={pd.DataFrame(rows)['group'].value_counts().to_dict()}")


def prepare_gse110636(mp: pd.DataFrame):
    meta = parse_series_matrix(DATA / "GSE110636/matrix/GSE110636_series_matrix.txt.gz")
    files = extract_tar_members(
        DATA / "GSE110636/suppl/GSE110636_RAW.tar",
        PROC / "GSE110636",
        lambda n: n.endswith("ReadsPerGene.out.tab.gz"),
    )
    cols, rows = {}, []
    for p in files:
        gsm = re.search(r"(GSM\d+)", p.name).group(1)
        hit = meta[meta["gsm"].eq(gsm)]
        title = hit["title"].iloc[0] if len(hit) else p.name
        subgroup = ""
        for c in hit.columns:
            val = str(hit.iloc[0][c])
            if val.lower().startswith("subgroup:"):
                subgroup = val.split(":", 1)[1].strip()
        cols[gsm] = read_star_rpg(p)
        rows.append({"gsm": gsm, "title": title, "subgroup": subgroup})
    counts = pd.DataFrame(cols).fillna(0).astype(np.int64)
    coldata = pd.DataFrame(rows)
    write_counts("gse110636", counts, coldata)
    sym = collapse_symbol(counts, mp)
    write_counts("gse110636_symbol", sym, coldata)


def prepare_gse186943():
    p = DATA / "GSE186943/suppl/GSE186943_FPKM_and_read_count_matrix_data_file.csv.gz"
    df = pd.read_csv(p)
    log(f"GSE186943 cols {len(df.columns)}: {list(df.columns)[:40]} ... {list(df.columns)[-8:]}")
    Path(PROC / "gse186943_colnames.txt").write_text("\n".join(map(str, df.columns)), encoding="utf-8")
    gene = "gene_symbol" if "gene_symbol" in df.columns else df.columns[2]
    exp_cols = [c for c in df.columns if str(c).startswith("exp_")]
    cnt_cols = [c for c in df.columns if str(c).startswith("cnt_") or str(c).startswith("count_")]
    fpkm = df.set_index(gene)[exp_cols].apply(pd.to_numeric, errors="coerce")
    fpkm = fpkm[~fpkm.index.duplicated(keep="first")]
    fpkm.index.name = "gene"
    fpkm.to_csv(PROC / "gse186943_fpkm.tsv.gz", sep="\t")
    coldata = []
    for c in exp_cols:
        s = str(c).replace("exp_", "")
        gt = "TG" if s.startswith("TG_") else ("WT" if s.startswith("WT_") else "other")
        stage = re.sub(r"^(TG|WT)_", "", s)
        stage = re.sub(r"_\d+N?$", "", stage)
        stage = re.sub(r"\d+N?$", "", stage)
        coldata.append({"sample": c, "short": s, "genotype": gt, "stage": stage})
    cd = pd.DataFrame(coldata)
    cd.to_csv(PROC / "gse186943_coldata.tsv", sep="\t", index=False)
    log(f"GSE186943 FPKM {fpkm.shape} stages={cd['stage'].value_counts().to_dict()}")
    if cnt_cols:
        cnt = df.set_index(gene)[cnt_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        cnt = cnt[~cnt.index.duplicated(keep="first")].round().astype(np.int64)
        cnt.to_csv(PROC / "gse186943_counts.tsv.gz", sep="\t")
        log(f"GSE186943 counts {cnt.shape}")
    # ChIP peaks vs Zeb1
    peaks = extract_tar_members(
        DATA / "GSE186943/suppl/GSE186943_RAW.tar",
        PROC / "GSE186943_chip",
        lambda n: n.endswith("narrowPeak.gz") or n.endswith(".narrowPeak"),
    )
    rows = []
    for pk in peaks:
        n = n_zeb = 0
        hits = []
        opener = gzip.open if str(pk).endswith(".gz") else open
        with opener(pk, "rt", errors="replace") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.split()
                if len(p) < 3:
                    continue
                chrom, start, end = p[0], int(p[1]), int(p[2])
                n += 1
                win = ZEB1_MM.get(chrom)
                if win and not (end < win[0] or start > win[1]):
                    n_zeb += 1
                    hits.append(f"{chrom}:{start}-{end}")
        rows.append({"file": pk.name, "n_peaks": n, "n_Zeb1_window": n_zeb, "hits": ";".join(hits[:12])})
        log(f"  peak {pk.name} n={n} Zeb1={n_zeb}")
    pd.DataFrame(rows).to_csv(PROC / "gse186943_chip_zeb1.tsv", sep="\t", index=False)


def prepare_gse144035():
    p = DATA / "GSE144035/suppl/GSE144035_All_counts_data.txt.gz"
    df = pd.read_csv(p, sep="\t")
    log(f"GSE144035 shape {df.shape} cols={list(df.columns)[:20]}")
    df.columns = [str(c).strip() for c in df.columns]
    gene = "Gene.Name" if "Gene.Name" in df.columns else df.columns[2]
    # precomputed DE if present
    keep_de = [c for c in df.columns if c in {"Gene.ID", "Gene.Name", "Biotype", "Dependent -Dox", "Independent -Dox", "FDR", "AveExpr", "P value"} or "Dox" in c or c in {"FDR", "AveExpr"}]
    if keep_de:
        de = df[keep_de].copy()
        de.to_csv(PROC / "gse144035_precomputed_de.tsv.gz", sep="\t", index=False)
    # sample count columns: typically VTL-*
    samp = [c for c in df.columns if str(c).startswith("VTL") or re.search(r"_S\d+$", str(c))]
    log(f"GSE144035 sample cols n={len(samp)} {samp[:12]}")
    if samp:
        mat = df.set_index(gene)[samp].apply(pd.to_numeric, errors="coerce").fillna(0)
        mat = mat[~mat.index.duplicated(keep="first")]
        # split integer-like vs fpkm-like: if max>1000 treat as counts
        if mat.max().max() > 1000:
            mat = mat.round().astype(np.int64)
        mat.to_csv(PROC / "gse144035_counts.tsv.gz", sep="\t")
        cd = pd.DataFrame({"sample": samp})
        cd.to_csv(PROC / "gse144035_coldata.tsv", sep="\t", index=False)
    soft = DATA / "GSE144035/soft/GSE144035_family.soft.gz"
    if soft.exists():
        with gzip.open(soft, "rt", errors="replace") as f:
            txt = f.read(20000)
        (PROC / "gse144035_soft_head.txt").write_text(txt[:15000], encoding="utf-8")


def prepare_gse188225():
    p = DATA / "GSE188225/suppl/GSE188225_Lmo2_transgenic_thymocyte_CellPop_GenewiseCounts.txt.gz"
    df = pd.read_csv(p, sep="\t")
    log(f"GSE188225 {df.shape} {list(df.columns)}")
    gene = df.columns[0]
    mat = df.set_index(gene)
    if "Length" in mat.columns:
        mat = mat.drop(columns=["Length"])
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0).round().astype(np.int64)
    mat.to_csv(PROC / "gse188225_counts.tsv.gz", sep="\t")
    meta = parse_series_matrix(DATA / "GSE188225/matrix/GSE188225_series_matrix.txt.gz")
    meta.to_csv(PROC / "gse188225_series_meta.tsv", sep="\t", index=False)
    cd = []
    for c in mat.columns:
        gen = "TG" if ".TG." in c or "_TG." in c or ".TG" in c else ("WT" if ".WT." in c or "_WT." in c else "unk")
        stage = "DN2" if "DN2" in c else ("DN3" if "DN3" in c else "other")
        age = "Old" if "Old" in c else ("Young" if "Young" in c else "unk")
        cd.append({"sample": c, "genotype": gen, "stage": stage, "age": age})
    pd.DataFrame(cd).to_csv(PROC / "gse188225_coldata.tsv", sep="\t", index=False)


def prepare_gse287751():
    tar = DATA / "GSE287751/suppl/GSE287751_RAW.tar"
    names = []
    with tarfile.open(tar, "r") as tf:
        for m in tf.getmembers():
            if m.isfile():
                names.append(f"{m.name}\t{m.size}")
    (PROC / "gse287751_tar_members.txt").write_text("\n".join(names), encoding="utf-8")
    log("GSE287751 members:\n" + "\n".join(names[:40]))
    files = extract_tar_members(
        tar, PROC / "GSE287751",
        lambda n: any(n.endswith(ext) for ext in (".txt.gz", ".tsv.gz", ".csv.gz", ".txt", ".tsv", ".csv", ".counts.txt.gz"))
        and "HTO" not in n and "cDNA" not in n.lower(),
    )
    log(f"GSE287751 extracted {len(files)}")
    for p in files[:8]:
        log(f"  {p.name} {p.stat().st_size}")


def prepare_m4_cores():
    st = pd.read_csv(AN / "module4/processed/target_zeb1_gene_stats.tsv", sep="\t")
    pos = st[(st["fdr_rho"] < 0.05) & (st["rho_zeb1"] > 0.3)]["gene"].astype(str)
    neg = st[(st["fdr_rho"] < 0.05) & (st["rho_zeb1"] < -0.3)]["gene"].astype(str)
    (PROC / "m4_zeb1_pos_core.txt").write_text("\n".join(pos), encoding="utf-8")
    (PROC / "m4_zeb1_neg_core.txt").write_text("\n".join(neg), encoding="utf-8")
    log(f"M4 cores pos={len(pos)} neg={len(neg)}")
    lab = pd.read_csv(AN / "module3/processed/target_tall_subtype.tsv", sep="\t")
    expr_p = AN / "module4/processed/target_log2tpm_filtered.tsv.gz"
    if expr_p.exists():
        expr = pd.read_csv(expr_p, sep="\t", index_col=0)
        common = [c for c in expr.columns if c in set(lab["usi"].astype(str))]
        lab = lab[lab["usi"].astype(str).isin(common)].copy()
        tlx = set(lab.loc[lab["subtype"].astype(str).eq("TLX"), "usi"].astype(str))
        rest = set(lab.loc[~lab["subtype"].astype(str).eq("TLX") & ~lab["subtype"].astype(str).eq("Unknown"), "usi"].astype(str))
        from scipy.stats import mannwhitneyu
        rows = []
        for g in expr.index:
            a = expr.loc[g, list(tlx)].astype(float)
            b = expr.loc[g, list(rest)].astype(float)
            a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
            if len(a) < 8 or len(b) < 8:
                continue
            try:
                u = mannwhitneyu(a, b, alternative="two-sided")
                p = u.pvalue
            except Exception:
                p = np.nan
            lfc = float(a.mean() - b.mean())
            rows.append((g, lfc, p, a.mean(), b.mean()))
        out = pd.DataFrame(rows, columns=["gene", "lfc_TLX_minus_rest", "p", "mean_TLX", "mean_rest"])
        out["signed_stat"] = np.sign(out["lfc_TLX_minus_rest"]) * (-np.log10(out["p"].clip(lower=1e-300)))
        out.to_csv(PROC / "target_tlx_vs_rest.tsv.gz", sep="\t", index=False)
        log(f"TARGET TLX vs rest genes={len(out)}")


def chip_human_zeb1():
    """Mean bigWig signal at ZEB1 if pyBigWig can be imported or installed."""
    try:
        import pyBigWig  # noqa: F401
    except Exception:
        import subprocess, sys
        r = subprocess.run([sys.executable, "-m", "pip", "install", "--user", "pyBigWig"], capture_output=True, text=True)
        log("pip pyBigWig rc=" + str(r.returncode) + " " + (r.stdout[-200:] if r.stdout else "") + (r.stderr[-300:] if r.stderr else ""))
        try:
            import pyBigWig
        except Exception as e:
            log(f"pyBigWig unavailable: {e}")
            pd.DataFrame([{"status": "skipped", "reason": str(e)}]).to_csv(PROC / "gse154675_chip_zeb1.tsv", sep="\t", index=False)
            return
    tar = DATA / "GSE154675/suppl/GSE154675_RAW.tar"
    want = ("ARR_aLMO2", "DU528_aLMO2", "HSB2_aLMO2", "CCRFCEM_aLMO2",
            "ARR_aTAL1", "DU528_aTAL1", "HSB2_aTAL1", "CCRFCEM_aTAL1")
    files = extract_tar_members(tar, PROC / "GSE154675_bw", lambda n: n.endswith(".bigwig") and any(k in n for k in want))
    rows = []
    for bw_p in files:
        try:
            bw = pyBigWig.open(str(bw_p))
        except Exception as e:
            log(f"open fail {bw_p.name}: {e}")
            continue
        chroms = list(bw.chroms().keys())[:8]
        log(f"  {bw_p.name} chroms={chroms}")
        rec = {"file": bw_p.name}
        for build, windows, tag in [("hg38", ZEB1_HG38, "gene"), ("hg19", ZEB1_HG19, "gene"),
                                    ("hg38p", ZEB1_PROM_HG38, "prom"), ("hg19p", ZEB1_PROM_HG19, "prom")]:
            chrom_map = windows
            val = np.nan
            for chrom, (a, b) in chrom_map.items():
                if chrom in bw.chroms():
                    try:
                        val = bw.stats(chrom, a, min(b, bw.chroms()[chrom]), type="mean")[0]
                    except Exception:
                        val = np.nan
                    break
            rec[f"{build}_{tag}"] = val
        bw.close()
        rows.append(rec)
        log(f"  occupancy {rec}")
    pd.DataFrame(rows).to_csv(PROC / "gse154675_chip_zeb1.tsv", sep="\t", index=False)


def safe(label, fn):
    try:
        fn()
        log(f"OK {label}")
    except Exception as e:
        log(f"FAIL {label}: {type(e).__name__}: {e}")


def main():
    mp = build_ensembl_map()

    def star_kd():
        if not (PROC / "gse110635_symbol_counts.tsv.gz").exists():
            c635, d635, _ = prepare_star_kd("GSE110635", class_110635)
            write_counts("gse110635", c635, d635)
            write_counts("gse110635_symbol", collapse_symbol(c635, mp), d635)
        else:
            log("skip GSE110635, already prepared")
        if not (PROC / "gse110632_symbol_counts.tsv.gz").exists():
            c632, d632, _ = prepare_star_kd("GSE110632", class_110632)
            write_counts("gse110632", c632, d632)
            write_counts("gse110632_symbol", collapse_symbol(c632, mp), d632)
        else:
            log("skip GSE110632, already prepared")

    safe("STAR KD", star_kd)
    safe("GSE6214x Agilent", prepare_gse6214x)
    safe("GSE110636", lambda: prepare_gse110636(mp))
    safe("GSE186943", prepare_gse186943)
    safe("GSE144035", prepare_gse144035)
    safe("GSE188225", prepare_gse188225)
    safe("GSE287751", prepare_gse287751)
    safe("M4 cores", prepare_m4_cores)
    safe("GSE154675 ChIP", chip_human_zeb1)
    log("PREPARE DONE")


if __name__ == "__main__":
    main()
