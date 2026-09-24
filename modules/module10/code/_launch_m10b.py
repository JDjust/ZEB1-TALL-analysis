#!/usr/bin/env python3
import subprocess
from pathlib import Path
HOST = "hulu-workstation"
loc = Path(r"d:\_bioinformation\modules\module10\code\01c_scrna.py")
rs = Path(r"d:\_bioinformation\modules\module10\code\02c_analyze_r2.R")
rem = "/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code"
for p in (loc, rs):
    p.write_text(p.read_text(encoding="utf-8").replace("\r\n", "\n"), encoding="utf-8")
    subprocess.check_call(["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", str(p), f"{HOST}:{rem}/{p.name}"])
    print("scp", p.name, flush=True)
cmd = r"""
sed -i 's/\r$//' /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/01c_scrna.py /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/code/02c_analyze_r2.R
# patient-level composition note (do not treat cells as n)
python3 - <<'PY'
from pathlib import Path
p=Path('/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/tables/M5_r2_composition_note.txt')
p.write_text(
'chi2 on cell counts is NOT a valid test (cells are not independent).\n'
'Use patient-level Spearman of mapped stage probabilities vs ZEB1 (M5_r2_ZEB1_vs_mapped.tsv).\n'
'n_patients=10. Nominal map_entropy P=0.025 does not pass FDR.\n',
encoding='utf-8')
PY
cd /data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10
nohup bash -lc 'PY=/opt/bioinfo/envs/omicverse-core/bin/python; RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
$PY code/01c_scrna.py > logs/01c_scrna.log 2>&1
echo PREP:$?
tail -n 30 logs/01c_scrna.log
$RS code/02c_analyze_r2.R > logs/02c_analyze_r2.log 2>&1
echo FIG:$?
tail -n 20 logs/02c_analyze_r2.log
echo M10_DONE
' > logs/run_m10b.log 2>&1 < /dev/null &
echo STARTED:$!
"""
p = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=25", HOST, cmd],
                   capture_output=True, timeout=60)
print((p.stdout or b"").decode("utf-8", "replace"))
print("rc", p.returncode, flush=True)
