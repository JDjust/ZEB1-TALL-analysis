#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

p = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/total_fig467F/tables/F7_target_surv_merged.tsv")
d = pd.read_csv(p, sep="\t")
print("surv shape", d.shape)
print(d.isna().mean().round(3).to_string())
print("os_event", d.os_event.value_counts(dropna=False).to_dict())
print("efs_event", d.efs_event.value_counts(dropna=False).to_dict())
print("os_time>0", int((d.os_time > 0).sum()), "efs_time>0", int((d.efs_time > 0).sum()))

pts = pd.read_csv(Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/total_fig467F/tables/F7_mtx_points.tsv"), sep="\t")
print("mtx points", pts.shape, pts.columns.tolist())
print(pts.head(3).to_string())
print("lineage", pts.lineage.value_counts(dropna=False).head(10).to_dict())
print("ZEB1 non-na", int(pts.ZEB1.notna().sum()))
print("sample cell_line", pts.cell_line.head(10).tolist())

clin = pd.read_csv("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed/target_clinical_combined.tsv", sep="\t")
print("First Event", clin["First Event"].value_counts(dropna=False).head(20).to_dict())
print("Vital Status", clin["Vital Status"].value_counts(dropna=False).to_dict())

dep = pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/Model.csv", usecols=["ModelID","CellLineName","StrippedCellLineName","OncotreePrimaryDisease","COSMICID","SangerModelID"])
print("model", dep.shape, "COSMIC non-na", dep.COSMICID.notna().sum())
tall = dep[dep.OncotreePrimaryDisease.eq("T-Lymphoblastic Leukemia/Lymphoma")]
print("T-ALL names", tall.StrippedCellLineName.tolist())
print("T-ALL COSMIC", tall.COSMICID.tolist())

gd = pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F/GDSC/GDSC2_fitted_dose_response_24Jul22.csv", nrows=2)
print("gdsc2 cols", list(gd.columns))
mtx = pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F/GDSC/GDSC2_fitted_dose_response_24Jul22.csv", usecols=["CELL_LINE_NAME","DRUG_NAME","COSMIC_ID","SANGER_MODEL_ID","LN_IC50"])
mtx = mtx[mtx.DRUG_NAME.astype(str).str.contains("ethotrexate", case=False, na=False)]
print("mtx n", len(mtx), "unique lines", mtx.CELL_LINE_NAME.nunique())
keys = set(tall.StrippedCellLineName.astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True))
gkeys = mtx.CELL_LINE_NAME.astype(str).str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
print("name overlap T-ALL", len(keys & set(gkeys)))
print("cosmic overlap", len(set(pd.to_numeric(tall.COSMICID, errors="coerce").dropna()) & set(pd.to_numeric(mtx.COSMIC_ID, errors="coerce").dropna())))
print("FINISH_INSPECT_OK")
