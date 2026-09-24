#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module9
mkdir -p "$REMOTE/logs" "$REMOTE/figures" "$REMOTE/tables" "$REMOTE/processed"
cd "$REMOTE"
"$PY" code/01_prepare_module9.py > logs/01_prepare_module9.log 2>&1
echo PREPARE_EXIT:$?
tail -n 80 logs/01_prepare_module9.log
"$RS" code/02_analyze_module9.R > logs/02_analyze_module9.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 80 logs/02_analyze_module9.log
echo FIG_N:$(ls figures 2>/dev/null | wc -l)
echo TAB_N:$(ls tables 2>/dev/null | wc -l)
