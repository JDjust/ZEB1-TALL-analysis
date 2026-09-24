#!/usr/bin/env python3
"""SCP the Fig4/6/7 F pipeline to hulu and start it in the background."""
import subprocess
from pathlib import Path

HOST = "hulu-workstation"
ROOT = Path(r"d:\_bioinformation")
loc = ROOT / "total" / "code" / "hulu_fig467F_pipeline.py"
rem_dir = "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/total_fig467F"
rem = rem_dir + "/hulu_fig467F_pipeline.py"
kg = ROOT / "modules" / "module3" / "tables" / "independent_keygenes.tsv"

text = loc.read_text(encoding="utf-8").replace("\r\n", "\n")
loc.write_text(text, encoding="utf-8", newline="\n")

subprocess.check_call([
    "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", HOST,
    f"mkdir -p {rem_dir}/logs {rem_dir}/tables",
])
subprocess.check_call([
    "scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25",
    str(loc), f"{HOST}:{rem}",
])
if kg.exists():
    subprocess.check_call([
        "scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25",
        str(kg), f"{HOST}:{rem_dir}/independent_keygenes.tsv",
    ])

cmd = (
    f"sed -i 's/\\r$//' {rem}; "
    f"cd {rem_dir}; "
    "nohup bash -lc '"
    "PY=/opt/bioinfo/envs/omicverse-core/bin/python; "
    "$PY hulu_fig467F_pipeline.py > logs/fig467F.log 2>&1"
    "' > logs/nohup_launch.log 2>&1 < /dev/null & "
    "echo STARTED:$!; "
    "sleep 1; "
    "head -n 5 logs/fig467F.log 2>/dev/null || true"
)
p = subprocess.run(
    ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", HOST, cmd],
    capture_output=True, timeout=90,
)
print((p.stdout or b"").decode("utf-8", "replace"))
print((p.stderr or b"").decode("utf-8", "replace"))
print("rc", p.returncode)
