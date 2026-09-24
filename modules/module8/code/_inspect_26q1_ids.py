#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
cmd = r"""
python3 - <<'PY'
import pandas as pd
from pathlib import Path
p=Path("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
with p.open() as f:
    h=f.readline()[:300]
    r=f.readline()[:300]
print("HDR", h)
print("ROW", r)
m=pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/Model.csv", nrows=3, usecols=lambda c: c in ["ModelID","CellLineName","StrippedCellLineName"])
print("MODEL", m.to_string(index=False))
c=pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/CRISPRGeneEffect.csv", nrows=1)
print("CRISPR_COL0", c.columns[0], c.iloc[0,0])
cnv=pd.read_csv("/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1/OmicsCNGeneWGS.csv", nrows=1)
print("CNV_COL0", cnv.columns[0], cnv.iloc[0,0])
print("M1FIG")
PY
ls /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/figures/M1_r2_* | head
echo M2
tail -n 12 /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/logs/01b_prepare_thymus.log
"""
p = subprocess.run(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=25",HOST,cmd],
                   capture_output=True, timeout=90)
Path(r"d:\_bioinformation\modules\module8\tables\_26q1_ids.txt").write_text(
    (p.stdout or b"").decode("utf-8","replace"), encoding="utf-8")
print((p.stdout or b"").decode("utf-8","replace").encode("ascii","replace").decode("ascii")[:4000])
print("rc", p.returncode)
