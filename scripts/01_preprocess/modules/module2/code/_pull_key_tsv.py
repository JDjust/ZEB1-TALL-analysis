#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HOST = "hulu-workstation"
files = [
    "M2_r2_GSE206710_donor_fraction.tsv",
    "M2_r2_GSE206710_stage_n.tsv",
    "M2_r2_GSE206710_sample_means.tsv",
    "M2_r2_map_GSE227122_patient.tsv",
    "M2_r2_map_GSE227122_stage_comp.tsv",
    "M2_r2_data_note.txt",
    "M2_r2_GSE195812_ZEB1_by_stage.tsv",
]
tab = Path(r"d:\_bioinformation\modules\module2\tables")
for name in files:
    loc = tab / name
    rem = f"/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/tables/{name}"
    p = subprocess.run(
        ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", f"hulu-workstation:{rem}", str(loc)],
        capture_output=True, timeout=60,
    )
    print(name, "rc", p.returncode, "size", loc.stat().st_size if loc.exists() else 0)
