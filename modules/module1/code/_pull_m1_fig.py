#!/usr/bin/env python3
import os, subprocess, sys
from pathlib import Path
os.environ["PYTHONUNBUFFERED"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
p = subprocess.run(
    ["ssh","-o","BatchMode=yes","-o","ConnectTimeout=20","hulu-workstation",
     "ls /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/figures/M1_r2_*"],
    capture_output=True, timeout=60,
)
names = [Path(x).name for x in (p.stdout or b"").decode("ascii","replace").split() if x]
print("n", len(names), flush=True)
dst = Path(r"d:\_bioinformation\modules\module1\figures")
dst.mkdir(parents=True, exist_ok=True)
for name in names:
    loc = dst / name
    if loc.exists() and loc.stat().st_size > 200:
        print("skip", name, flush=True)
        continue
    print("get", name, flush=True)
    r = subprocess.run(
        ["scp","-o","BatchMode=yes","-o","ConnectTimeout=20",
         f"hulu-workstation:/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/figures/{name}",
         str(loc)], capture_output=True, timeout=75)
    print(" ", "ok" if r.returncode==0 else "FAIL", r.returncode, flush=True)
