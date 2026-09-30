#!/usr/bin/env python3
"""Pull M2/M1 r2 files one-by-one so a dropped glob scp cannot kill the batch."""
import subprocess
import sys
from pathlib import Path

HOST = "hulu-workstation"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def ssh(cmd, timeout=90):
    p = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25",
         "-o", "ServerAliveInterval=15", "-o", "ServerAliveCountMax=8", HOST, cmd],
        capture_output=True, timeout=timeout,
    )
    return p.returncode, (p.stdout or b"").decode("utf-8", "replace")

def scp_one(remote, local, timeout=120):
    Path(local).parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        r = subprocess.run(
            ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", remote, str(local)],
            capture_output=True, timeout=timeout,
        )
        if r.returncode == 0 and Path(local).exists() and Path(local).stat().st_size > 0:
            return True
        print(f"  retry {attempt+1} {Path(remote).name} rc={r.returncode}")
    return False

rc, listing = ssh(
    "ls /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/tables/M2_r2_* "
    "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/figures/M2_r2_* "
    "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1/figures/M1_r2_* 2>/dev/null"
)
print("list rc", rc)
files = [ln.strip() for ln in listing.splitlines() if ln.strip()]
print("n remote", len(files))

# tables first (small), then pdf, then png, skip huge gz until last
def rank(p):
    n = Path(p).name
    if n.endswith(".tsv") or n.endswith(".txt"):
        return 0
    if n.endswith(".pdf"):
        return 1
    if n.endswith(".png"):
        return 2
    return 3

files = sorted(files, key=rank)
ok = fail = 0
for rem in files:
    name = Path(rem).name
    if "/module1/" in rem.replace("\\", "/"):
        loc = Path(r"d:\_bioinformation\modules\module1\figures") / name
    elif "/figures/" in rem.replace("\\", "/"):
        loc = Path(r"d:\_bioinformation\modules\module2\figures") / name
    else:
        loc = Path(r"d:\_bioinformation\modules\module2\tables") / name
    if loc.exists() and loc.stat().st_size > 50:
        print("skip", name)
        ok += 1
        continue
    to = 240 if name.endswith(".gz") else 90
    print("get", name)
    if scp_one(f"{HOST}:{rem}", loc, timeout=to):
        ok += 1
        print("  ok", loc.stat().st_size)
    else:
        fail += 1
        print("  FAIL", name)
print("done ok", ok, "fail", fail)
