#!/usr/bin/env python3
"""Prepare Module 1 expression + labels from MILE, TARGET, Pharmacotype."""
from __future__ import annotations

import gzip
import json
import math
import re
import ssl
import time
import urllib.request
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1")
PROC = ROOT / "processed"
PROC.mkdir(parents=True, exist_ok=True)
CTX = ssl.create_default_context()
UA = {"User-Agent": "ZEB1-TALL-analysis/1.0 (academic)"}
GENES = ["ZEB1", "ZEB2", "LMO2"]
ALIASES = {
    "ZEB1": {"ZEB1", "TCF8", "ZEB-1", "AREB6"},
    "ZEB2": {"ZEB2", "SIP1", "ZFHX1B", "SMADIP1"},
    "LMO2": {"LMO2", "RBTN2", "RHOM2", "TTG2"},
}


def log(msg: str) -> None:
    print(msg, flush=True)


def hedges_g(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float]:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return (np.nan, np.nan, np.nan, np.nan)
    s1, s2 = a.var(ddof=1), b.var(ddof=1)
    sp = math.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if sp == 0:
        return (0.0, 0.0, 0.0, 0.0)
    d = (a.mean() - b.mean()) / sp
    J = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    g = J * d
    se = math.sqrt((n1 + n2) / (n1 * n2) + (d * d) / (2.0 * (n1 + n2))) * J
    return g, se, g - 1.96 * se, g + 1.96 * se


def classify_mile(raw: str) -> str:
    s = raw.strip()
    sl = s.lower()
    if s == "T-ALL" or sl == "t-all":
        return "T-ALL"
    if "non-leukemia" in sl or "healthy" in sl:
        return "Normal BM"
    if "aml" in sl:
        return "AML"
    if sl == "cml" or "cml" == sl:
        return "CML"
    if sl == "cll":
        return "CLL"
    if sl == "mds" or sl.startswith("mds"):
        return "MDS"
    if "all" in sl:
        return "B-ALL"
    return "Other"


def parse_mile_metadata() -> pd.DataFrame:
    brief = DATA / "GSE13159" / "GSE13159_gsm_brief.txt"
    rows = []
    cur = {"gsm": None, "title": "", "source": "", "sample_type": "", "leukemia_class": ""}
    with brief.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur["gsm"]:
                    rows.append(cur)
                cur = {
                    "gsm": line.split("=", 1)[-1].strip(),
                    "title": "",
                    "source": "",
                    "sample_type": "",
                    "leukemia_class": "",
                }
            elif line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[-1].strip()
            elif line.startswith("!Sample_source_name_ch1"):
                cur["source"] = line.split("=", 1)[-1].strip()
            elif line.startswith("!Sample_characteristics_ch1"):
                v = line.split("=", 1)[-1].strip()
                if v.lower().startswith("sample type:"):
                    cur["sample_type"] = v.split(":", 1)[-1].strip()
                elif v.lower().startswith("leukemia class:"):
                    cur["leukemia_class"] = v.split(":", 1)[-1].strip()
        if cur["gsm"]:
            rows.append(cur)
    df = pd.DataFrame(rows)
    df["disease"] = df["leukemia_class"].map(classify_mile)
    df["cohort"] = "GSE13159"
    df["platform"] = "Affymetrix HG-U133 Plus 2.0"
    log(f"MILE samples {len(df)}")
    log(df["leukemia_class"].value_counts().to_string())
    log(df["disease"].value_counts().to_string())
    return df


def parse_gpl570() -> dict[str, str]:
    """probe -> gene symbol (Entrez Gene symbol column)."""
    path = DATA / "GSE13159" / "GPL570.annot.gz"
    probe2gene: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        header = None
        idx_id, idx_sym = 0, None
        for line in f:
            if line.startswith("#") or line.startswith("^") or (line.startswith("!") and "platform_table_begin" not in line):
                continue
            if line.startswith("!platform_table_begin"):
                continue
            if header is None:
                parts = line.rstrip("\n").split("\t")
                header = parts
                for i, p in enumerate(parts):
                    pl = p.lower().strip()
                    if pl == "id":
                        idx_id = i
                    if pl == "gene symbol":
                        idx_sym = i
                if idx_sym is None:
                    raise SystemExit(f"Gene symbol column not found in {header[:12]}")
                log(f"GPL570 header symbol col={header[idx_sym]!r} idx={idx_sym}")
                continue
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_id, idx_sym):
                continue
            pid = parts[idx_id].strip().strip('"')
            sym = parts[idx_sym].strip().strip('"')
            if not pid or not sym or sym in {"---", ""}:
                continue
            gene = re.split(r"///|;|,", sym)[0].strip()
            probe2gene[pid] = gene
    log(f"GPL570 mapped probes {len(probe2gene)}")
    return probe2gene


def wanted_probes(probe2gene: dict[str, str]) -> dict[str, set[str]]:
    out = {g: set() for g in GENES}
    for probe, gene in probe2gene.items():
        gu = gene.upper()
        for g, als in ALIASES.items():
            if gu in als:
                out[g].add(probe)
    # canonical HG-U133 Plus 2.0 fallbacks
    fallback = {
        "ZEB1": {"212764_at", "210875_s_at", "239948_at"},
        "ZEB2": {"203603_s_at", "228337_at"},
        "LMO2": {"204249_s_at", "204248_at"},
    }
    for g, ps in fallback.items():
        if len(out[g]) == 0:
            out[g] |= ps
            log(f"{g} using fallback probes")
    for g, ps in out.items():
        log(f"{g} probes n={len(ps)} {sorted(ps)}")
    return out


def extract_mile_expr(probe_map: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sm = DATA / "GSE13159" / "GSE13159_series_matrix.txt.gz"
    keep = set().union(*probe_map.values())
    probe_rows = {}
    gsms = None
    with gzip.open(sm, "rt", encoding="utf-8", errors="replace") as f:
        in_table = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                header = next(f)
                gsms = [x.strip().strip('"') for x in header.rstrip("\n").split("\t")[1:]]
                log(f"MILE matrix samples {len(gsms)}")
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            pid = line.split("\t", 1)[0].strip().strip('"')
            if pid in keep:
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                probe_rows[pid] = np.array([float(v) if v not in {"", "NA", "null"} else np.nan for v in vals], dtype=float)
    log(f"extracted probes {len(probe_rows)}")
    gene_mat = {}
    for g, probes in probe_map.items():
        mats = [probe_rows[p] for p in probes if p in probe_rows]
        if not mats:
            raise SystemExit(f"no probes extracted for {g}")
        gene_mat[g] = np.nanmean(np.vstack(mats), axis=0)
    expr = pd.DataFrame(gene_mat, index=gsms)
    expr.index.name = "gsm"
    probe_df = pd.DataFrame(probe_rows, index=gsms)
    probe_df.index.name = "gsm"
    return expr, probe_df


def gdc_map_rna_files(file_names: list[str]) -> pd.DataFrame:
    url = "https://api.gdc.cancer.gov/files"
    rows = []
    chunk = 80
    for i in range(0, len(file_names), chunk):
        names = file_names[i : i + chunk]
        payload = {
            "filters": {"op": "in", "content": {"field": "file_name", "value": names}},
            "fields": ",".join(
                [
                    "file_id",
                    "file_name",
                    "cases.submitter_id",
                    "cases.samples.submitter_id",
                    "cases.samples.sample_type",
                    "cases.samples.tissue_type",
                    "cases.project.project_id",
                    "cases.diagnoses.primary_diagnosis",
                    "cases.diagnoses.age_at_diagnosis",
                ]
            ),
            "size": chunk,
        }
        req = urllib.request.Request(
            url,
            data=json.dumps(payload).encode(),
            headers={**UA, "Content-Type": "application/json"},
            method="POST",
        )
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, context=CTX, timeout=60) as resp:
                    obj = json.loads(resp.read().decode())
                break
            except Exception as e:
                log(f"GDC retry {attempt} {e}")
                time.sleep(2 * (attempt + 1))
        else:
            raise SystemExit("GDC API failed")
        hits = obj.get("data", {}).get("hits", [])
        log(f"GDC mapped {i}-{i+len(names)} hits={len(hits)}")
        for h in hits:
            case0 = (h.get("cases") or [{}])[0]
            samples = case0.get("samples") or [{}]
            diag = (case0.get("diagnoses") or [{}])[0]
            # prefer RNA aliquot-like sample if multiple
            samp = samples[0]
            for s in samples:
                sid = str(s.get("submitter_id") or "")
                if sid.endswith("R") or "-01R" in sid or "-09A" in sid or "-03A" in sid:
                    samp = s
                    break
            rows.append(
                {
                    "file_id": h.get("file_id"),
                    "file_name": h.get("file_name"),
                    "case_submitter_id": case0.get("submitter_id"),
                    "sample_submitter_id": samp.get("submitter_id"),
                    "sample_type": samp.get("sample_type"),
                    "tissue_type": samp.get("tissue_type"),
                    "project_id": (case0.get("project") or {}).get("project_id"),
                    "primary_diagnosis": diag.get("primary_diagnosis"),
                    "age_at_diagnosis": diag.get("age_at_diagnosis"),
                }
            )
    return pd.DataFrame(rows)


def extract_star_genes(path: Path, genes=GENES) -> dict[str, float]:
    want = set(genes)
    out = {}
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or line.startswith("gene_id") or line.startswith("N_"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            name = parts[1]
            if name in want:
                try:
                    out[name] = float(parts[6])  # tpm_unstranded
                except ValueError:
                    out[name] = np.nan
                if len(out) == len(want):
                    break
    return out


def usi_from_barcode(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x)
    m = re.search(r"TARGET-\d+-([A-Z0-9]+)", s)
    if m:
        return m.group(1)
    if re.fullmatch(r"[A-Z0-9]{5,8}", s):
        return s
    return s


def prepare_target() -> pd.DataFrame:
    rna_dir = DATA / "TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts"
    files = sorted(p for p in rna_dir.iterdir() if p.is_file() and p.name.endswith(".tsv"))
    log(f"TARGET RNA files {len(files)}")
    map_path = PROC / "target_gdc_file_map.tsv"
    if map_path.exists():
        fmap = pd.read_csv(map_path, sep="\t")
        log(f"reuse GDC map {len(fmap)}")
    else:
        fmap = gdc_map_rna_files([p.name for p in files])
        fmap.to_csv(map_path, sep="\t", index=False)
    recs = []
    fmap_idx = fmap.drop_duplicates("file_name").set_index("file_name")
    for i, p in enumerate(files, 1):
        expr = extract_star_genes(p)
        rec = {"file_name": p.name, **expr}
        if p.name in fmap_idx.index:
            rec.update(fmap_idx.loc[p.name].to_dict())
        recs.append(rec)
        if i % 50 == 0:
            log(f"TARGET extracted {i}/{len(files)}")
    df = pd.DataFrame(recs)
    clin = pd.read_csv(PROC / "target_clinical_combined.tsv", sep="\t")
    clin["usi"] = clin["TARGET USI"].map(usi_from_barcode)

    df["usi"] = df["sample_submitter_id"].map(usi_from_barcode)
    df["usi2"] = df["case_submitter_id"].map(usi_from_barcode)
    df["usi"] = np.where(df["usi"].astype(str).isin(["", "nan", "None"]), df["usi2"], df["usi"])
    keep_cols = [c for c in ["usi", "Protocol", "Gender", "Age at Diagnosis in Days", "Cell of Origin"] if c in clin.columns]
    prot = clin.drop_duplicates("usi")[keep_cols]
    df = df.merge(prot, on="usi", how="left")
    tall_sup = PROC / "target_t_all_supplement.tsv"
    tall_usi = set()
    if tall_sup.exists():
        sup = pd.read_csv(tall_sup, sep="\t")
        if "TARGET USI" in sup.columns:
            tall_usi = set(sup["TARGET USI"].map(usi_from_barcode))

    def lineage_row(r) -> str:
        usi = str(r.get("usi") or "")
        proto = str(r.get("Protocol") or "")
        coo = str(r.get("Cell of Origin") or "").lower()
        diag = str(r.get("primary_diagnosis") or "").lower()
        if "t cell" in coo or coo.strip() in {"t", "t-all", "t all"}:
            return "T-ALL"
        if "b cell" in coo or "b-precursor" in coo or "b precursor" in coo or coo.strip() in {"b", "b-all"}:
            return "B-ALL"
        if usi in tall_usi:
            return "T-ALL"
        if "AALL0434" in proto:
            return "T-ALL"
        if any(x in proto for x in ["AALL0232", "AALL0331", "P9906", "9906", "AALL03B1", "AALL07P4", "AALL08P1"]):
            return "B-ALL"
        if "t lymphoblastic" in diag or "t-cell" in diag or "t cell" in diag:
            return "T-ALL"
        if "b lymphoblastic" in diag or "b-cell" in diag or "precursor b" in diag:
            return "B-ALL"
        return "Unknown"

    df["disease"] = df.apply(lineage_row, axis=1)
    df["cohort"] = "TARGET-ALL-P2"
    df["platform"] = "RNA-seq STAR TPM"
    log(df["disease"].value_counts(dropna=False).to_string())
    log(df["sample_type"].value_counts(dropna=False).head(10).to_string())
    return df


def _norm_id(s) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def _pharm_disease(x: str) -> str:
    xu = str(x).upper()
    if re.search(r"T[- ]?ALL|\bT\b.*LINEAGE|T-LINEAGE|^T$", xu):
        return "T-ALL"
    if re.search(r"B[- ]?ALL|\bB\b.*LINEAGE|B-LINEAGE|^B$", xu):
        return "B-ALL"
    if xu.strip() in {"T", "TALL"}:
        return "T-ALL"
    if xu.strip() in {"B", "BALL"}:
        return "B-ALL"
    return "Other"


def prepare_pharmacotype() -> pd.DataFrame:
    zpath = DATA / "StJude_ALL_Pharmacotype/pharmacotyping_ped_rnaseq_fpkm.zip"
    with zipfile.ZipFile(zpath) as z:
        names = z.namelist()
        log(f"FPKM zip entries {names[:10]}")
        info = max(z.infolist(), key=lambda x: x.file_size)
        log(f"using {info.filename} size={info.file_size}")
        with z.open(info.filename) as fh:
            peek = fh.read(2000).decode("utf-8", "replace")
        log("FPKM peek: " + peek[:400].replace("\n", " | "))
        first = peek.splitlines()[0]
        sep = "\t" if first.count("\t") > first.count(",") else ","
        with z.open(info.filename) as fh:
            mat = pd.read_csv(fh, sep=sep)
    log(f"FPKM matrix {mat.shape} cols {list(mat.columns)[:8]}")
    gene_col = mat.columns[0]
    for c in mat.columns[:6]:
        vals = mat[c].astype(str).str.upper()
        if vals.isin(GENES).any():
            gene_col = c
            break
    mat = mat.rename(columns={gene_col: "gene"})
    mat["gene_u"] = mat["gene"].astype(str).str.upper().str.replace(r"\.\d+$", "", regex=True)
    sub = mat[mat["gene_u"].isin(GENES)].copy()
    if sub.empty:
        for c in mat.columns:
            cand = mat[c].astype(str).str.upper()
            if cand.isin(GENES).any():
                mat["gene_u"] = cand
                sub = mat[mat["gene_u"].isin(GENES)].copy()
                break
    if sub.empty:
        raise SystemExit("ZEB1/ZEB2/LMO2 not found in pharmacotype FPKM")
    keep = [c for c in sub.columns if c not in {"gene", "gene_u"} and pd.api.types.is_numeric_dtype(sub[c])]
    expr = sub.set_index("gene_u")[keep].T
    expr.index.name = "sample"
    expr = expr.reindex(columns=GENES)
    t1 = pd.read_csv(PROC / "pharmacotype_supp_table1.tsv", sep="\t")
    t2 = pd.read_csv(PROC / "pharmacotype_supp_table2.tsv", sep="\t")
    log(f"pharm labels {t1.shape} immuno={t1['Immunophenotype'].value_counts().to_dict() if 'Immunophenotype' in t1.columns else None}")
    lab = t1.copy()
    lab["lineage_raw"] = lab["Immunophenotype"].astype(str) if "Immunophenotype" in lab.columns else ""
    lab["disease"] = lab["lineage_raw"].map(_pharm_disease)
    lab["k_patient"] = lab["Patient ID"].map(_norm_id)
    if "Sample ID" in t2.columns:
        t2["k_sample"] = t2["Sample ID"].map(_norm_id)
        t2["k_patient"] = t2["Patient ID"].map(_norm_id)
        lab = lab.merge(t2[["k_patient", "k_sample", "Sample ID"]], on="k_patient", how="left")
    else:
        lab["k_sample"] = lab["k_patient"]
    expr = expr.reset_index()
    expr["k"] = expr["sample"].map(_norm_id)
    merged = expr.merge(lab, left_on="k", right_on="k_patient", how="left")
    n1 = merged["disease"].isin(["T-ALL", "B-ALL"]).sum()
    if n1 < 50 and "k_sample" in lab.columns:
        merged2 = expr.merge(lab, left_on="k", right_on="k_sample", how="left")
        n2 = merged2["disease"].isin(["T-ALL", "B-ALL"]).sum()
        log(f"pharm match patient={n1} sample={n2}")
        if n2 > n1:
            merged = merged2
    merged["cohort"] = "StJude_Pharmacotype"
    merged["platform"] = "RNA-seq FPKM"
    log(merged["disease"].value_counts(dropna=False).to_string())
    log(f"pharmacotype T/B matched {merged['disease'].isin(['T-ALL','B-ALL']).sum()} / {len(merged)}")
    return merged


def main() -> None:
    mile_path = PROC / "gse13159_keygenes.tsv"
    if mile_path.exists():
        log(f"reuse {mile_path}")
    else:
        meta = parse_mile_metadata()
        meta.to_csv(PROC / "gse13159_sample_metadata.tsv", sep="\t", index=False)
        p2g = parse_gpl570()
        pmap = wanted_probes(p2g)
        expr, probe_df = extract_mile_expr(pmap)
        mile = meta.merge(expr.reset_index(), on="gsm", how="inner")
        mile.to_csv(mile_path, sep="\t", index=False)
        probe_df.to_csv(PROC / "gse13159_keygene_probes.tsv", sep="\t")
        log(f"MILE merged {mile.shape}")

    tgt = prepare_target()
    tgt.to_csv(PROC / "target_keygenes.tsv", sep="\t", index=False)

    ph = prepare_pharmacotype()
    ph.to_csv(PROC / "pharmacotype_keygenes.tsv", sep="\t", index=False)
    log("PREPARE DONE")


if __name__ == "__main__":
    main()
