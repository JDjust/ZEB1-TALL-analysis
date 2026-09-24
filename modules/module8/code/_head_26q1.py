#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
cmd = r"""
head -c 400 /data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv; echo; echo ---
head -c 200 /data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/OmicsCNGeneWGS.csv; echo; echo ---
head -c 200 /data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/CRISPRGeneEffect.csv; echo
echo M2
tail -n 15 /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/logs/01b_prepare_thymus.log
pgrep -af '01b_prepare_thymus' || true
"""
p = subprocess.run(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=25",HOST,cmd],
                   capture_output=True, timeout=90)
print((p.stdout or b"").decode("utf-8","replace").encode("ascii","replace").decode("ascii")[:3500])
print("rc", p.returncode)
