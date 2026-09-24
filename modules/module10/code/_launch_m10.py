#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
subprocess.check_call(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=25",HOST,
                       "mkdir -p /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code"])
pairs = [
    (r"d:\_bioinformation\modules\module10\code\01_prepare_module10.py",
     "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/01_prepare_module10.py"),
    (r"d:\_bioinformation\modules\module10\code\02_analyze_module10.R",
     "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/02_analyze_module10.R"),
    (r"d:\_bioinformation\modules\module10\code\run_module10.sh",
     "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/run_module10.sh"),
]
for loc, rem in pairs:
    Path(loc).write_text(Path(loc).read_text(encoding="utf-8").replace("\r\n", "\n"), encoding="utf-8")
    subprocess.check_call(["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", loc, f"{HOST}:{rem}"])
print("scp", Path(loc).name)
cmd = r"""
mkdir -p /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/{code,logs,figures,tables,processed}
sed -i 's/\r$//' /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/*.py \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/*.R \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/*.sh
chmod +x /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/run_module10.sh
cd /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10
nohup bash code/run_module10.sh > logs/run_module10.log 2>&1 < /dev/null &
echo STARTED:$!
"""
p = subprocess.run(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=25",HOST,cmd],
                   capture_output=True, timeout=90)
print((p.stdout or b"").decode("utf-8","replace"))
print("rc", p.returncode)
