#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module7
cd "$REMOTE"
mkdir -p logs figures tables processed
"$PY" code/01_prepare_module7.py > logs/01_prepare_module7.log 2>&1
echo PREPARE_EXIT:$?
tail -n 30 logs/01_prepare_module7.log
"$RS" code/02_analyze_module7.R > logs/02_analyze_module7.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 40 logs/02_analyze_module7.log
echo FIG_N:$(ls figures | wc -l)
