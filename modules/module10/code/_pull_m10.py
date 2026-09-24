#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
Path(r"d:\_bioinformation\modules\module10\tables").mkdir(parents=True, exist_ok=True)
Path(r"d:\_bioinformation\modules\module10\figures").mkdir(parents=True, exist_ok=True)
subprocess.call(["scp","-o","BatchMode=yes","-o","ConnectTimeout=25",
    f"{HOST}:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/tables/M10_*",
    r"d:\_bioinformation\modules\module10\tables\\"])
subprocess.call(["scp","-o","BatchMode=yes","-o","ConnectTimeout=25",
    f"{HOST}:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/figures/M10_*",
    r"d:\_bioinformation\modules\module10\figures\\"])
print("tables", list(Path(r"d:\_bioinformation\modules\module10\tables").glob("M10_*")))
print("figs", list(Path(r"d:\_bioinformation\modules\module10\figures").glob("M10_*")))
