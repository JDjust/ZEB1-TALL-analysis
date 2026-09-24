#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
# tables
subprocess.call(["scp","-o","BatchMode=yes","-o","ConnectTimeout=25",
    f"{HOST}:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8/tables/M8_r2_*",
    r"d:\_bioinformation\modules\module8\tables\\"])
subprocess.call(["scp","-o","BatchMode=yes","-o","ConnectTimeout=25",
    f"{HOST}:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8/tables/M8_8.1_8.5_line_table.tsv",
    r"d:\_bioinformation\modules\module8\tables\M8_26q1_line_table.tsv"])
subprocess.call(["scp","-o","BatchMode=yes","-o","ConnectTimeout=25",
    f"{HOST}:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8/logs/01b_prepare_depmap26q1.log",
    r"d:\_bioinformation\modules\module8\tables\_26q1_prep.log"])
print("pulled")
print(list(Path(r"d:\_bioinformation\modules\module8\tables").glob("M8_r2*")))
print(list(Path(r"d:\_bioinformation\modules\module8\tables").glob("*26q1*")))
