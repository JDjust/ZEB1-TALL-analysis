#!/usr/bin/env python3
"""Rebuild Module 3 labels: Pharmacotype subtype/ETP and GSE272023 CIMP mapping.

TARGET mutation/CNV table is kept. Pharmacotype is taken from module1 keygenes
(already merged); the previous module3 re-merge blanked Molecular subtype.
CIMP is mapped via Sample_description TALL### in gsm_brief.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

M1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
PROC = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed")
DATA = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data")
PROC.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    print(msg, flush=True)


def fusion_class(lesion) -> str:
    if pd.isna(lesion):
        return "unknown"
    s = str(lesion).strip()
    if s.lower() in {"wt", "wildtype", "wild type", "none", ""}:
        return "no_fusion"
    u = s.upper()
    if "ZEB1" in u:
        return "ZEB1_fusion"
    if "TAL1" in u or "TAL2" in u:
        return "TAL"
    if "TLX1" in u or "TLX3" in u:
        return "TLX"
    if "LMO1" in u or "LMO2" in u:
        return "LMO"
    if "HOXA" in u:
        return "HOXA"
    if "KMT2A" in u or re.search(r"\bMLL\b", u):
        return "KMT2A"
    if "MLLT10" in u:
        return "MLLT10"
    if "NKX2" in u:
        return "NKX2-1"
    return "other_fusion"


def fix_target() -> pd.DataFrame:
    t = pd.read_csv(PROC / "target_tall_subtype.tsv", sep="\t")
    t["fusion_class"] = t["Lesion"].map(fusion_class) if "Lesion" in t.columns else "unknown"
    t["zeb1_fusion"] = t["Lesion"].astype(str).str.contains("ZEB1", case=False, na=False)
    if "Age at Diagnosis in Days" in t.columns:
        t["age_years"] = pd.to_numeric(t["Age at Diagnosis in Days"], errors="coerce") / 365.25
    elif "age_at_diagnosis" in t.columns:
        t["age_years"] = pd.to_numeric(t["age_at_diagnosis"], errors="coerce") / 365.25
    t["lmo2_subtype"] = t["subtype"].eq("LMO2/LYL1")
    t["tal_subtype"] = t["subtype"].eq("TAL")
    t["hoxa_subtype"] = t["subtype"].eq("HOXA")
    t["tlx_subtype"] = t["subtype"].eq("TLX")
    log(f"TARGET n={len(t)}")
    log("fusion\n" + t["fusion_class"].value_counts(dropna=False).to_string())
    log(f"ZEB1 fusion n={int(t['zeb1_fusion'].sum())}")
    log(f"age_years finite {int(t['age_years'].notna().sum())} range {t['age_years'].min():.2f}-{t['age_years'].max():.2f}")
    t.to_csv(PROC / "target_tall_subtype.tsv", sep="\t", index=False)
    return t


def collapse_pharm_subtype(x) -> str:
    s = str(x).strip()
    if s.lower() in {"nan", "none", ""}:
        return "Unknown"
    if s.upper() == "ETP" or "ETP" in s.upper():
        return "ETP"
    if s.upper() in {"T-ALL", "TALL", "T"}:
        return "T-ALL (non-ETP)"
    return s


def fix_pharmacotype() -> pd.DataFrame:
    ph = pd.read_csv(M1 / "pharmacotype_keygenes.tsv", sep="\t")
    ph = ph[ph["disease"].eq("T-ALL")].copy()
    sub_col = "Molecular subtype" if "Molecular subtype" in ph.columns else None
    if sub_col is None:
        raise SystemExit("Molecular subtype missing from pharmacotype_keygenes")
    ph["subtype"] = ph[sub_col].map(collapse_pharm_subtype)
    ph["etp"] = np.where(ph["subtype"].eq("ETP"), "ETP", "notETP")
    ph["age_years"] = pd.to_numeric(ph.get("Age at diagnosis (years)"), errors="coerce")
    ph["cohort"] = "StJude_Pharmacotype"
    keep = [
        "sample", "ZEB1", "ZEB2", "LMO2", "Patient ID", "Immunophenotype",
        "Molecular subtype", "subtype", "etp", "age_years",
        "Day 15 MRD (%)", "Day 42 or 46 MRD (%)", "NCI risk", "Protocol",
        "disease", "cohort", "platform",
    ]
    keep = [c for c in keep if c in ph.columns]
    out = ph[keep].copy()
    log(f"Pharmacotype T-ALL n={len(out)}")
    log("subtype\n" + out["subtype"].value_counts(dropna=False).to_string())
    log("etp\n" + out["etp"].value_counts(dropna=False).to_string())
    out.to_csv(PROC / "pharmacotype_tall_subtype.tsv", sep="\t", index=False)
    return out


def parse_cimp_brief() -> pd.DataFrame:
    brief = DATA / "GSE272023/GSE272023_gsm_brief.txt"
    rows = []
    cur = {"gsm": "", "title": "", "sample": "", "cimp": "", "tissue_type": "", "source": ""}
    with brief.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("^SAMPLE"):
                if cur["gsm"]:
                    rows.append(cur)
                cur = {"gsm": line.split("=", 1)[-1].strip(), "title": "", "sample": "",
                       "cimp": "", "tissue_type": "", "source": ""}
            elif line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[-1].strip()
            elif line.startswith("!Sample_source_name_ch1"):
                cur["source"] = line.split("=", 1)[-1].strip()
            elif line.startswith("!Sample_description"):
                v = line.split("=", 1)[-1].strip()
                if re.fullmatch(r"TALL\d+", v, flags=re.I) or re.fullmatch(r"T-ALL\d+", v, flags=re.I):
                    cur["sample"] = v.upper().replace("-", "")
                elif not cur["sample"] and re.search(r"TALL\d+", v, flags=re.I):
                    cur["sample"] = re.search(r"TALL\d+", v, flags=re.I).group(0).upper()
            elif line.startswith("!Sample_characteristics"):
                v = line.split("=", 1)[-1].strip()
                low = v.lower()
                if "cimp subgroup" in low:
                    cur["cimp"] = v.split(":", 1)[-1].strip()
                elif low.startswith("tissue type:"):
                    cur["tissue_type"] = v.split(":", 1)[-1].strip()
        if cur["gsm"]:
            rows.append(cur)
    meta = pd.DataFrame(rows)
    if meta["sample"].eq("").any():
        meta.loc[meta["sample"].eq(""), "sample"] = meta["title"].str.extract(r"(TALL\d+)", flags=re.I, expand=False)
        meta["sample"] = meta["sample"].astype(str).str.upper()
    n_map = int(meta["sample"].str.match(r"TALL\d+", na=False).sum())
    log(f"CIMP meta n={len(meta)} mapped_sample={n_map}")
    log("cimp\n" + meta["cimp"].value_counts(dropna=False).to_string())
    log("tissue\n" + meta["tissue_type"].value_counts(dropna=False).to_string())
    return meta


def extract_vst_genes(samples_wanted: list[str] | None = None) -> pd.DataFrame:
    mat_path = DATA / "GSE272023/GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz"
    genes = {"ZEB1", "ZEB2", "LMO2"}
    ensg = {"ENSG00000148516": "ZEB1", "ENSG00000169554": "ZEB2", "ENSG00000135363": "LMO2"}
    with gzip.open(mat_path, "rt", encoding="utf-8", errors="replace") as fh:
        header = [h.strip().strip('"') for h in fh.readline().rstrip("\n").split("\t")]
        samples = header[1:]
        found = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            gid = parts[0].strip().strip('"')
            gu = gid.upper().split(".")[0]
            name = None
            if gu in genes:
                name = gu
            elif gid.split(".")[0] in ensg:
                name = ensg[gid.split(".")[0]]
            if name and name not in found:
                found[name] = np.array([float(x) if x not in {"", "NA"} else np.nan for x in parts[1:]], float)
            if len(found) == 3:
                break
    expr = pd.DataFrame(found, index=samples)
    expr.index.name = "sample"
    expr = expr.reset_index()
    log(f"VST genes {list(found)} n_samples={len(expr)}")
    return expr


def fix_cimp() -> pd.DataFrame:
    meta = parse_cimp_brief()
    expr = extract_vst_genes()
    df = expr.merge(meta[["sample", "gsm", "title", "cimp", "tissue_type", "source"]], on="sample", how="left")
    df["cimp"] = df["cimp"].astype(str).str.replace(r"^cimp subgroup:\s*", "", regex=True, flags=re.I)
    df["cimp_g"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df.loc[df["cimp"].str.contains("high", case=False, na=False), "cimp_g"] = "CIMP-high"
    df.loc[df["cimp"].str.contains("low", case=False, na=False), "cimp_g"] = "CIMP-low"
    df["is_tumor"] = df["tissue_type"].astype(str).str.contains("tumor", case=False, na=False)
    df["cohort"] = "GSE272023"
    log("merged cimp_g\n" + df["cimp_g"].value_counts(dropna=False).to_string())
    log("tumor vs not " + str(df["is_tumor"].value_counts(dropna=False).to_dict()))
    df.to_csv(PROC / "gse272023_cimp.tsv", sep="\t", index=False)
    meta.to_csv(PROC / "gse272023_sample_map.tsv", sep="\t", index=False)
    return df


def main() -> None:
    fix_target()
    fix_pharmacotype()
    fix_cimp()
    log("MODULE3 LABELS FIXED")


if __name__ == "__main__":
    main()
