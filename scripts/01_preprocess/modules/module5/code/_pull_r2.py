#!/usr/bin/env python3
import subprocess
from pathlib import Path

HOST = "hulu-workstation"
REMOTE = "/data-b/liangfuhua/projects/ZEB1_TALL_analysis"
LOCAL = Path(r"d:\_bioinformation\modules")

def ssh(cmd, timeout=40):
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", HOST, cmd],
                       capture_output=True, timeout=timeout)
    return (p.stdout or b"").decode("utf-8", "replace")

# list remote r2 files
listing = ssh(
    "ls /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module{4,5,6,7,10}/figures/M*_r2_* "
    "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module{4,5,6,7,10}/tables/M*_r2_* "
    "2>/dev/null"
)
names = [ln.strip() for ln in listing.splitlines() if ln.strip()]
print("remote files", len(names), flush=True)
ok = fail = skip = 0
for rem in names:
    rel = rem.split("/ZEB1_TALL_analysis/")[-1]
    loc = LOCAL / rel.replace("/", "\\") if False else LOCAL / Path(*Path(rel).parts)
    loc.parent.mkdir(parents=True, exist_ok=True)
    if loc.exists() and loc.stat().st_size > 200 and loc.suffix.lower() in {".pdf", ".png"}:
        skip += 1
        continue
    r = subprocess.run(
        ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", f"{HOST}:{rem}", str(loc)],
        capture_output=True, timeout=90,
    )
    if r.returncode == 0 and loc.exists() and loc.stat().st_size > 50:
        print("ok", loc.name, loc.stat().st_size, flush=True)
        ok += 1
    else:
        print("FAIL", rem, (r.stderr or b"")[-80:].decode("ascii", "replace"), flush=True)
        fail += 1
print("done ok", ok, "skip", skip, "fail", fail, flush=True)
