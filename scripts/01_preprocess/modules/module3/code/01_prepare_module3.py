#!/usr/bin/env python3
"""Prepare Module 3: ZEB1 vs T-ALL subtype / ETP / CIMP / mutation context."""
from __future__ import annotations

import gzip
import json
import re
import time
import urllib.request
import ssl
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
M1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3")
PROC = ROOT / "processed"
PROC.mkdir(parents=True, exist_ok=True)
CTX = ssl.create_default_context()
UA = {"User-Agent": "ZEB1-TALL-analysis/1.0 (academic)"}
GENES = ["ZEB1", "ZEB2", "LMO2"]
MUT_GENES = ["NOTCH1", "PTEN", "NRAS", "KRAS", "NF1", "JAK1", "JAK3", "IL7R", "FLT3", "FBXW7", "PHF6", "WT1", "ZEB1", "ZEB2", "LMO2"]


def log(msg: str) -> None:
    print(msg, flush=True)


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


def collapse_group(g: str) -> str:
    s = str(g).strip()
    if s in {"TAL1", "TAL2"}:
        return "TAL"
    if s in {"LMO1/2", "LMO2_LYL1"}:
        return "LMO2/LYL1"
    if s in {"TLX1", "TLX3"}:
        return "TLX"
    if s in {"HOXA", "NKX2_1"}:
        return s
    if s in {"Unknown", "", "nan"}:
        return "Unknown"
    return s


def prepare_target() -> pd.DataFrame:
    tgt = pd.read_csv(M1 / "target_keygenes.tsv", sep="\t")
    tgt = tgt[tgt["disease"].eq("T-ALL")].copy()
    tgt = tgt[tgt["sample_type"].astype(str).str.startswith("Primary")]
    tgt["pref"] = np.where(tgt["sample_type"].astype(str).str.contains("Bone Marrow"), 1, 2)
    tgt["usi"] = tgt["usi"].map(usi_from_barcode)
    tgt = tgt.sort_values(["usi", "pref"]).drop_duplicates("usi")
    sup = pd.read_csv(M1 / "target_t_all_supplement.tsv", sep="\t")
    sup["usi"] = sup["TARGET USI"].map(usi_from_barcode)
    keep = ["usi", "Group", "Lesion", "ETP status", "Maturation stage", "BM Blasts at Diagnosis"]
    keep = [c for c in keep if c in sup.columns]
    df = tgt.merge(sup[keep], on="usi", how="left")
    df["subtype"] = df["Group"].map(collapse_group)
    df["etp"] = df["ETP status"].astype(str).replace({"nan": np.nan, "NA": np.nan, "Unevaluable": np.nan})
    df["lmo2_subtype"] = df["subtype"].eq("LMO2/LYL1")
    df["tal_subtype"] = df["subtype"].eq("TAL")
    df["hoxa_subtype"] = df["subtype"].eq("HOXA")
    df["cohort"] = "TARGET-ALL-P2"
    log(f"TARGET T-ALL primary unique {len(df)}")
    log(str(df["subtype"].value_counts(dropna=False)))
    log(str(df["etp"].value_counts(dropna=False)))
    return df


def prepare_pharmacotype() -> pd.DataFrame:
    ph = pd.read_csv(M1 / "pharmacotype_keygenes.tsv", sep="\t")
    ph = ph[ph["disease"].eq("T-ALL")].copy()
    t1 = pd.read_csv(M1 / "pharmacotype_supp_table1.tsv", sep="\t")
    t1.columns = [c.strip() for c in t1.columns]
    log("pharm t1 cols " + str(list(t1.columns)[:15]))
    sub_col = next((c for c in t1.columns if "molecular" in c.lower() and "subtype" in c.lower()), None)
    t1["k_patient"] = t1["Patient ID"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    ph["k"] = ph["sample"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True) if "sample" in ph.columns else ph.iloc[:, 0].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    if "k" in ph.columns and ph["k"].isna().all() is False:
        pass
    if "k_patient" in ph.columns:
        ph["k"] = ph["k_patient"].astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    keep = ["k_patient"]
    if sub_col:
        keep.append(sub_col)
    if "Age at diagnosis (years)" in t1.columns:
        keep.append("Age at diagnosis (years)")
    lab = t1[keep].drop_duplicates("k_patient")
    df = ph.merge(lab, left_on="k", right_on="k_patient", how="left")
    if sub_col and sub_col in df.columns:
        df["subtype"] = df[sub_col].astype(str)
    else:
        df["subtype"] = "Unknown"
    df["etp"] = np.where(df["subtype"].str.contains("ETP", case=False, na=False), "ETP", "notETP")
    df["cohort"] = "StJude_Pharmacotype"
    log(f"Pharmacotype T-ALL {len(df)}")
    log(str(df["subtype"].value_counts(dropna=False).head(20)))
    return df


def prepare_gse272023() -> pd.DataFrame:
    brief = DATA / "GSE272023/GSE272023_gsm_brief.txt"
    rows = []
    cur = {"gsm": None, "cimp": "", "title": ""}
    with brief.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("^SAMPLE"):
                if cur["gsm"]:
                    rows.append(cur)
                cur = {"gsm": line.split("=", 1)[-1].strip(), "cimp": "", "title": ""}
            elif line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[-1].strip()
            elif "cimp subgroup" in line.lower():
                cur["cimp"] = line.split(":", 1)[-1].strip()
        if cur["gsm"]:
            rows.append(cur)
    meta = pd.DataFrame(rows)
    log(f"GSE272023 meta {len(meta)} cimp={meta['cimp'].value_counts().to_dict()}")

    mat_path = DATA / "GSE272023/GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz"
    with gzip.open(mat_path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        header = [h.strip().strip('"') for h in header]
        samples = header[1:]
        found = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            gid = parts[0].strip().strip('"')
            gu = gid.upper().split(".")[0]
            if gu in GENES or gid in GENES:
                found[gu if gu in GENES else gid] = np.array([float(x) if x not in {"", "NA"} else np.nan for x in parts[1:]], float)
            if len(found) == 3:
                break
    if "ZEB1" not in found:
        # ENSG lookup in a second pass
        ensg = {"ENSG00000148516": "ZEB1", "ENSG00000169554": "ZEB2", "ENSG00000135363": "LMO2"}
        with gzip.open(mat_path, "rt", encoding="utf-8", errors="replace") as fh:
            _ = fh.readline()
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                gid = parts[0].strip().strip('"').split(".")[0]
                if gid in ensg and ensg[gid] not in found:
                    found[ensg[gid]] = np.array([float(x) if x not in {"", "NA"} else np.nan for x in parts[1:]], float)
                if len(found) == 3:
                    break
    log(f"GSE272023 genes {list(found)}")
    expr = pd.DataFrame(found, index=samples)
    expr.index.name = "sample"
    expr = expr.reset_index()
    # GSM vs TALL201 mapping: series matrix is tiny; sample IDs in VST are TALL###
    expr["cimp"] = expr["sample"].map(dict(zip(meta["title"], meta["cimp"])))
    if expr["cimp"].isna().all():
        # try gsm_brief title contains sample
        title_map = {}
        for r in meta.itertuples(index=False):
            t = str(r.title)
            for s in expr["sample"]:
                if s in t:
                    title_map[s] = r.cimp
        if title_map:
            expr["cimp"] = expr["sample"].map(title_map)
    # series matrix sample titles
    sm = DATA / "GSE272023/matrix/GSE272023_series_matrix.txt.gz"
    if sm.exists() and expr["cimp"].isna().mean() > 0.5:
        with gzip.open(sm, "rt", encoding="utf-8", errors="replace") as f:
            gsm, titles, cimps = None, None, None
            for line in f:
                if line.startswith("!Sample_geo_accession"):
                    gsm = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                if line.startswith("!Sample_title"):
                    titles = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
                if "cimp" in line.lower() and line.startswith("!"):
                    cimps = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
        if titles and cimps:
            tmap = {t: c for t, c in zip(titles, cimps)}
            expr["cimp"] = expr["sample"].map(tmap).fillna(expr["cimp"])
            gsm_map = {t: g for t, g in zip(titles, gsm)} if gsm else {}
            expr["gsm"] = expr["sample"].map(gsm_map)
        if gsm and cimps and expr["cimp"].isna().mean() > 0.5:
            # characteristics may be "cimp subgroup: CIMP high"
            expr = expr.merge(meta.rename(columns={"gsm": "gsm"}), on="gsm", how="left", suffixes=("", "_m"))
            if "cimp_m" in expr.columns:
                expr["cimp"] = expr["cimp"].fillna(expr["cimp_m"])
    expr["cimp"] = expr["cimp"].astype(str).str.replace(r"^cimp subgroup:\s*", "", regex=True, flags=re.I)
    expr["cohort"] = "GSE272023"
    log(str(expr["cimp"].value_counts(dropna=False)))
    return expr


def gdc_map_maf(file_names: list[str]) -> pd.DataFrame:
    url = "https://api.gdc.cancer.gov/files"
    rows = []
    for i in range(0, len(file_names), 80):
        names = file_names[i : i + 80]
        payload = {
            "filters": {"op": "in", "content": {"field": "file_name", "value": names}},
            "fields": "file_id,file_name,cases.submitter_id,cases.samples.submitter_id",
            "size": 80,
        }
        req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers={**UA, "Content-Type": "application/json"}, method="POST")
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, context=CTX, timeout=60) as resp:
                    obj = json.loads(resp.read().decode())
                break
            except Exception as e:
                log(f"GDC maf retry {attempt} {e}")
                time.sleep(2 * (attempt + 1))
        else:
            continue
        for h in obj.get("data", {}).get("hits", []):
            case0 = (h.get("cases") or [{}])[0]
            samp = (case0.get("samples") or [{}])[0]
            rows.append({"file_name": h.get("file_name"), "case_submitter_id": case0.get("submitter_id"), "sample_submitter_id": samp.get("submitter_id")})
        log(f"GDC maf mapped {i}-{i+len(names)} n={len(rows)}")
    return pd.DataFrame(rows)


def prepare_mutations(tall_usi: set[str]) -> pd.DataFrame:
    maf_dir = DATA / "TARGET-ALL-P2/Masked_Somatic_Mutation/Aliquot_Ensemble_Somatic_Variant_Merging_and_Masking"
    files = sorted(p for p in maf_dir.iterdir() if p.suffixes[-2:] == [".maf", ".gz"] or p.name.endswith(".maf.gz"))
    log(f"MAF files {len(files)}")
    map_path = PROC / "target_maf_gdc_map.tsv"
    if map_path.exists():
        fmap = pd.read_csv(map_path, sep="\t")
    else:
        fmap = gdc_map_maf([p.name for p in files])
        fmap.to_csv(map_path, sep="\t", index=False)
    fmap["usi"] = fmap["case_submitter_id"].map(usi_from_barcode)
    by_name = fmap.drop_duplicates("file_name").set_index("file_name")
    hits = []
    for i, p in enumerate(files, 1):
        usi = by_name.loc[p.name, "usi"] if p.name in by_name.index else ""
        if usi not in tall_usi:
            continue
        try:
            maf = pd.read_csv(p, sep="\t", comment="#", low_memory=False, usecols=lambda c: c in {"Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"})
        except Exception:
            continue
        maf = maf[maf["Hugo_Symbol"].isin(MUT_GENES)]
        maf = maf[~maf["Variant_Classification"].isin(["Silent", "Intron", "IGR", "RNA", "lincRNA", "3'UTR", "5'UTR", "3'Flank", "5'Flank"])]
        for g in maf["Hugo_Symbol"].unique():
            hits.append({"usi": usi, "gene": g, "mut": 1})
        if i % 100 == 0:
            log(f"MAF scanned {i}/{len(files)}")
    if not hits:
        return pd.DataFrame(columns=["usi"] + MUT_GENES)
    long = pd.DataFrame(hits).drop_duplicates()
    wide = long.pivot_table(index="usi", columns="gene", values="mut", aggfunc="max", fill_value=0)
    if isinstance(wide.columns, pd.MultiIndex):
        wide.columns = wide.columns.get_level_values(-1)
    wide = wide.reset_index()
    wide = wide.rename(columns={g: f"{g}_mut" for g in MUT_GENES if g in wide.columns})
    for g in MUT_GENES:
        col = f"{g}_mut"
        if col not in wide.columns:
            wide[col] = 0
    ras = [c for c in ["NRAS_mut", "KRAS_mut", "NF1_mut", "JAK1_mut", "JAK3_mut", "IL7R_mut", "FLT3_mut"] if c in wide.columns]
    wide["RAS_JAK"] = (wide[ras].sum(axis=1) > 0).astype(int) if ras else 0
    log(f"mutation-positive T-ALL patients {len(wide)}")
    log("mut counts " + str({c: int(wide[c].sum()) for c in wide.columns if c.endswith("_mut") or c == "RAS_JAK"}))
    return wide


def prepare_cnv(tall_usi: set[str]) -> pd.DataFrame:
    cnv_dir = DATA / "TARGET-ALL-P2/Gene_Level_Copy_Number/ASCAT2"
    files = sorted(p for p in cnv_dir.rglob("*") if p.is_file() and p.stat().st_size > 1000)
    log(f"CNV gene-level files {len(files)}")
    # peek one
    if not files:
        return pd.DataFrame()
    peek = files[0]
    log(f"CNV peek {peek.name}")
    head = peek.read_text(encoding="utf-8", errors="replace").splitlines()[:5]
    for h in head:
        log("  " + h[:200])
    # GDC gene level: gene_id, gene_name, copy_number, min_copy_number, max_copy_number
    recs = []
    # map filenames via existing RNA map if same case, else GDC
    rna_map = pd.read_csv(M1 / "target_gdc_file_map.tsv", sep="\t")
    # CNV files likely UUID named; query GDC
    map_path = PROC / "target_cnv_gdc_map.tsv"
    if map_path.exists():
        fmap = pd.read_csv(map_path, sep="\t")
    else:
        fmap = gdc_map_maf([p.name for p in files[:800]])
        fmap.to_csv(map_path, sep="\t", index=False)
    fmap["usi"] = fmap["case_submitter_id"].map(usi_from_barcode)
    by_name = fmap.drop_duplicates("file_name").set_index("file_name") if len(fmap) else None
    for i, p in enumerate(files, 1):
        usi = ""
        if by_name is not None and p.name in by_name.index:
            usi = by_name.loc[p.name, "usi"]
        if usi not in tall_usi:
            continue
        try:
            d = pd.read_csv(p, sep="\t")
        except Exception:
            continue
        name_col = "gene_name" if "gene_name" in d.columns else ("Gene Symbol" if "Gene Symbol" in d.columns else None)
        cn_col = "copy_number" if "copy_number" in d.columns else None
        if name_col is None:
            continue
        if cn_col is None:
            for c in d.columns:
                if "copy" in c.lower() or c.lower() in {"value", "segment_mean"}:
                    cn_col = c
                    break
        if cn_col is None:
            continue
        sub = d[d[name_col].isin(["ZEB1", "ZEB2", "LMO2"])]
        rec = {"usi": usi}
        # burden: fraction |log2| or copy != 2
        cn = pd.to_numeric(d[cn_col], errors="coerce")
        rec["cnv_burden"] = float(np.mean(np.abs(cn - 2) > 0.5)) if cn.notna().any() else np.nan
        for _, r in sub.iterrows():
            rec[f"{r[name_col]}_cn"] = r[cn_col]
        recs.append(rec)
        if i % 50 == 0:
            log(f"CNV {i}/{len(files)} kept={len(recs)}")
    log(f"CNV T-ALL matched {len(recs)}")
    return pd.DataFrame(recs).drop_duplicates("usi") if recs else pd.DataFrame()


def main() -> None:
    tgt_path = PROC / "target_tall_subtype.tsv"
    if tgt_path.exists():
        log(f"reuse {tgt_path}")
        tgt = pd.read_csv(tgt_path, sep="\t")
    else:
        tgt = prepare_target()
        mut = prepare_mutations(set(tgt["usi"]))
        if len(mut):
            tgt = tgt.merge(mut, on="usi", how="left")
            for g in [f"{x}_mut" for x in MUT_GENES] + ["RAS_JAK"]:
                if g in tgt.columns:
                    tgt[g] = tgt[g].fillna(0)
        cnv = prepare_cnv(set(tgt["usi"]))
        if len(cnv):
            tgt = tgt.merge(cnv, on="usi", how="left")
        tgt.to_csv(tgt_path, sep="\t", index=False)

    ph = prepare_pharmacotype()
    ph.to_csv(PROC / "pharmacotype_tall_subtype.tsv", sep="\t", index=False)

    cimp = prepare_gse272023()
    cimp.to_csv(PROC / "gse272023_cimp.tsv", sep="\t", index=False)
    log("MODULE3 PREPARE DONE")


if __name__ == "__main__":
    main()
