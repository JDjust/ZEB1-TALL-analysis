#!/usr/bin/env bash
set -uo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10
mkdir -p "$REMOTE"/{logs,figures,tables,processed,code}
cd "$REMOTE"
"$PY" code/01_prepare_module10.py > logs/01_prepare_module10.log 2>&1
echo PREP:$?
tail -n 40 logs/01_prepare_module10.log
"$RS" code/02_analyze_module10.R > logs/02_analyze_module10.log 2>&1
echo FIG:$?
tail -n 20 logs/02_analyze_module10.log
echo M10_DONE
