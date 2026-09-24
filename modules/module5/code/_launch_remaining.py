#!/usr/bin/env python3
"""Upload remaining round-2 scripts to hulu and launch sequentially."""
import subprocess
from pathlib import Path

HOST = "hulu-workstation"
ROOT = Path(r"d:\_bioinformation\modules")
REMOTE = "/data-b/liangfuhua/projects/ZEB1_TALL_analysis"

pairs = [
    (ROOT / "module5/code/01c_round2.py", f"{REMOTE}/module5/code/01c_round2.py"),
    (ROOT / "module5/code/02c_analyze_r2.R", f"{REMOTE}/module5/code/02c_analyze_r2.R"),
    (ROOT / "module5/code/run_remaining_r2.sh", f"{REMOTE}/module5/code/run_remaining_r2.sh"),
    (ROOT / "module6/code/01c_chip_scatac.py", f"{REMOTE}/module6/code/01c_chip_scatac.py"),
    (ROOT / "module6/code/02c_analyze_r2.R", f"{REMOTE}/module6/code/02c_analyze_r2.R"),
    (ROOT / "module7/code/02c_analyze_r2.R", f"{REMOTE}/module7/code/02c_analyze_r2.R"),
    (ROOT / "module4/code/06_leftover_r2.R", f"{REMOTE}/module4/code/06_leftover_r2.R"),
    (ROOT / "module10/code/01c_scrna.py", f"{REMOTE}/module10/code/01c_scrna.py"),
    (ROOT / "module10/code/02c_analyze_r2.R", f"{REMOTE}/module10/code/02c_analyze_r2.R"),
]

for loc, rem in pairs:
    loc.parent.mkdir(parents=True, exist_ok=True)
    text = loc.read_text(encoding="utf-8").replace("\r\n", "\n")
    loc.write_text(text, encoding="utf-8")
    subprocess.check_call(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", HOST,
                           f"mkdir -p {Path(rem).parent.as_posix()}"])
    subprocess.check_call(["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", str(loc), f"{HOST}:{rem}"])
    print("scp", loc.name, flush=True)

cmd = r"""
sed -i 's/\r$//' /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module*/code/*r2* \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/code/run_remaining_r2.sh \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/code/06_leftover_r2.R \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/code/01c_round2.py \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6/code/01c_chip_scatac.py \
  /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/01c_scrna.py
chmod +x /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/code/run_remaining_r2.sh
mkdir -p /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs
cd /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5
nohup bash code/run_remaining_r2.sh > logs/run_remaining_r2.log 2>&1 < /dev/null &
echo STARTED:$!
"""
p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", HOST, cmd],
                   capture_output=True, timeout=90)
print((p.stdout or b"").decode("utf-8", "replace"))
print("rc", p.returncode, flush=True)
