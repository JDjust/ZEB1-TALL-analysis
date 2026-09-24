#!/usr/bin/env python3
"""Prepare Module 7: Pharmacotype T-ALL drug LC50 + ZEB1."""
from pathlib import Path

import numpy as np
import pandas as pd

M1 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/processed")
M3 = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3/processed")
ROOT = Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module7")
PROC = ROOT / "processed"
PROC.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    print(msg, flush=True)


def main() -> None:
    ph = pd.read_csv(M1 / "pharmacotype_keygenes.tsv", sep="\t")
    log(f"pharmacotype all {len(ph)} disease={ph['disease'].value_counts(dropna=False).to_dict()}")
    drug_cols = [c for c in ph.columns if c.endswith(")") and any(u in c for u in ["IU/ml", "nM", "µM", "uM", "μM"])]
    norm_cols = [c for c in ph.columns if c.endswith("_normalized")]
    log("raw LC50 cols " + str(drug_cols))
    log("normalized cols " + str(norm_cols))
    mtx = [c for c in ph.columns if "methotrexate" in c.lower() or c.lower().startswith("mtx") or "mtx" in c.lower()]
    log("MTX-like columns " + str(mtx))

    keep_id = ["sample", "Patient ID", "ZEB1", "ZEB2", "LMO2", "disease", "Immunophenotype",
               "Molecular subtype", "Age at diagnosis (years)", "Day 15 MRD (%)",
               "Day 42 or 46 MRD (%)", "NCI risk", "Protocol", "Sex"]
    keep_id = [c for c in keep_id if c in ph.columns]
    out = ph[keep_id + drug_cols + norm_cols].copy()
    out["ZEB1_log"] = np.log2(out["ZEB1"] + 1)
    out["ZEB2_log"] = np.log2(out["ZEB2"] + 1)
    out["LMO2_log"] = np.log2(out["LMO2"] + 1)
    out["ratio_zeb"] = np.log2((out["ZEB1"] + 1e-6) / (out["ZEB2"] + 1e-6))
    if (M3 / "pharmacotype_tall_subtype.tsv").exists():
        lab = pd.read_csv(M3 / "pharmacotype_tall_subtype.tsv", sep="\t")
        add = [c for c in ["sample", "subtype", "etp"] if c in lab.columns]
        out = out.merge(lab[add].drop_duplicates("sample"), on="sample", how="left")
    out.to_csv(PROC / "pharmacotype_all_drugs.tsv", sep="\t", index=False)
    tall = out[out["disease"].eq("T-ALL")].copy()
    tall.to_csv(PROC / "pharmacotype_tall_drugs.tsv", sep="\t", index=False)
    ball = out[out["disease"].eq("B-ALL")].copy()
    ball.to_csv(PROC / "pharmacotype_ball_drugs.tsv", sep="\t", index=False)
    log(f"T-ALL {len(tall)} B-ALL {len(ball)}")
    n_drug = {c: int(pd.to_numeric(tall[c], errors="coerce").notna().sum()) for c in drug_cols}
    log("T-ALL n with LC50 " + str(n_drug))
    note = PROC / "M7_MTX_note.txt"
    if not mtx:
        note.write_text(
            "MTX is not among the 18 drugs in St. Jude Pharmacotype (PMID 35145305 / 41591_2022_2112).\n"
            "The panel is: asparaginase, bortezomib, CHZ868, cytarabine, dasatinib, daunorubicin,\n"
            "dexamethasone, ibrutinib, mercaptopurine, nelarabine, panobinostat, prednisolone,\n"
            "ruxolitinib, thioguanine, trametinib, venetoclax, vincristine, vorinostat.\n"
            "Closest antimetabolites: mercaptopurine and thioguanine. Nelarabine is T-ALL directed.\n"
            "7.3-7.6 MTX therefore cannot be run on this resource; 7.7 screens the 18-drug panel.\n",
            encoding="utf-8",
        )
        log("wrote MTX absence note")
    log("MODULE7 PREPARE DONE")


if __name__ == "__main__":
    main()
