#!/usr/bin/env python3
import os, subprocess, sys
from pathlib import Path
os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HOST = "hulu-workstation"

pdfs = [
    ("module2/figures", [
        "M2_r2_2.1_umap_facs.pdf", "M2_r2_2.2_umap_ZEB1.pdf", "M2_r2_2.3_umap_LMO2.pdf",
        "M2_r2_2.4_ZEB1_facs_violin.pdf", "M2_r2_2.5_gene_trajectories.pdf",
        "M2_r2_2.6_gene_vs_stage.pdf", "M2_r2_2.7_stage_heatmap.pdf",
        "M2_r2_2.8_umap_206_donor.pdf", "M2_r2_2.9_donor_ZEB1_fraction.pdf",
        "M2_r2_2.10_patient_majority_stage.pdf", "M2_r2_2.11_patient_stage_stack.pdf",
        "M2_r2_2.12_residual_entropy.pdf", "M2_r2_2.13_TALL_ZEB1_by_mapped_stage.pdf",
        "M2_r2_2.14_dpt_density.pdf", "M2_r2_paga_195812.pdf", "M2_r2_paga_206710.pdf",
        "M2_r2_2.1_umap_facs.png", "M2_r2_2.8_umap_206_donor.png",
        "M2_r2_2.9_donor_ZEB1_fraction.png", "M2_r2_2.12_residual_entropy.png",
    ]),
    ("module1/figures", [
        "M1_r2_1.5_ZEB1_ZEB2_phase.pdf", "M1_r2_1.6_celltype_heatmap.pdf",
        "M1_r2_1.1_celltype_rank.pdf", "M1_r2_1.4_Tscore_scatter.pdf",
        "M1_r2_1.9_median_rank.pdf",
    ]),
]
base = Path(r"d:\_bioinformation\modules")
ok = fail = 0
for rel, names in pdfs:
    for name in names:
        loc = base / rel.split("/")[0] / rel.split("/")[1] / name
        loc.parent.mkdir(parents=True, exist_ok=True)
        if loc.exists() and loc.stat().st_size > 200:
            print("skip", name, loc.stat().st_size, flush=True)
            ok += 1
            continue
        rem = f"{HOST}:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/{rel}/{name}"
        print("get", name, flush=True)
        r = subprocess.run(
            ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", rem, str(loc)],
            capture_output=True, timeout=75,
        )
        if r.returncode == 0 and loc.exists() and loc.stat().st_size > 200:
            print("  ok", loc.stat().st_size, flush=True)
            ok += 1
        else:
            print("  FAIL", r.returncode, (r.stderr or b"")[-120:].decode("ascii","replace"), flush=True)
            fail += 1
print("done ok", ok, "fail", fail, flush=True)
