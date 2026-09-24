#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
p = Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_CCLE")
m = pd.read_csv(p / "Model.csv", usecols=lambda c: c in {
    "ModelID", "CellLineName", "StrippedCellLineName", "DepmapModelType",
    "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode", "Age"
})
print("TALL", (m.OncotreePrimaryDisease == "T-Lymphoblastic Leukemia/Lymphoma").sum())
print("BALL", (m.OncotreePrimaryDisease == "B-Lymphoblastic Leukemia/Lymphoma").sum())
print(m.loc[m.OncotreePrimaryDisease == "T-Lymphoblastic Leukemia/Lymphoma",
            ["ModelID", "CellLineName", "OncotreeSubtype", "DepmapModelType"]].to_string(index=False))
exp = p / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
with open(exp) as fh:
    header = fh.readline().rstrip("\n").split(",")
print("expr_ncols", len(header), "rowname", header[0][:40])
hits = [h for h in header if h.startswith("ZEB1 ") or h.startswith("ZEB2 ") or h.startswith("LMO2 ")]
print("expr_hits", hits[:10])
cr = p / "CRISPRGeneEffect.csv"
with open(cr) as fh:
    h2 = fh.readline().rstrip("\n").split(",")
print("crispr_ncols", len(h2), "rowname", h2[0][:40])
ch = [h for h in h2 if h.startswith("ZEB1 ") or h.startswith("ZEB2 ") or h.startswith("LMO2 ")]
print("crispr_hits", ch)
print("has_prob", (p / "CRISPRGeneDependency.csv").exists())
print("has_cnv", any(x.name.lower().find("cnv") >= 0 for x in p.iterdir()))
print("has_prism", any("prism" in x.name.lower() or "drug" in x.name.lower() for x in p.iterdir()))
print("files", [x.name for x in p.iterdir()])
