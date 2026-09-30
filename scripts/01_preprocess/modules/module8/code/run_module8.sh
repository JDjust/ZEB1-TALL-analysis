#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module8
mkdir -p "$REMOTE/logs" "$REMOTE/figures" "$REMOTE/tables" "$REMOTE/processed"
cd "$REMOTE"
"$PY" code/01_prepare_module8.py > logs/01_prepare_module8.log 2>&1
echo PREPARE_EXIT:$?
tail -n 40 logs/01_prepare_module8.log
"$RS" code/02_analyze_module8.R > logs/02_analyze_module8.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 40 logs/02_analyze_module8.log
echo FIG_N:$(ls figures | wc -l)
