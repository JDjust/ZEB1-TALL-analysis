#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
loc = r"d:\_bioinformation\modules\module8\code\01b_prepare_depmap26q1.py"
rem = "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8/code/01b_prepare_depmap26q1.py"
Path(loc).write_text(Path(loc).read_text(encoding="utf-8").replace("\r\n", "\n"), encoding="utf-8")
subprocess.check_call(["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", loc, f"{HOST}:{rem}"])
cmd = r"""
sed -i 's/\r$//' /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8/code/01b_prepare_depmap26q1.py
cd /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8
nohup bash -lc 'PY=/opt/bioinfo/envs/omicverse-core/bin/python; RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
$PY code/01b_prepare_depmap26q1.py > logs/01b_prepare_depmap26q1.log 2>&1
echo PREP:$?
tail -n 30 logs/01b_prepare_depmap26q1.log
$RS code/02_analyze_module8.R > logs/02_analyze_module8_26q1.log 2>&1
echo FIG:$?
tail -n 20 logs/02_analyze_module8_26q1.log
ls tables/M8_r2_* 2>/dev/null | wc -l
echo M8_26Q1_DONE
' > logs/run_26q1.log 2>&1 < /dev/null &
echo STARTED:$!
"""
p = subprocess.run(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=25",HOST,cmd],
                   capture_output=True, timeout=90)
text = (p.stdout or b"").decode("utf-8","replace")
print(text.encode("ascii","replace").decode("ascii"))
print("rc", p.returncode)
