#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5
mkdir -p "$REMOTE/logs" "$REMOTE/figures" "$REMOTE/tables" "$REMOTE/processed"
cd "$REMOTE"
"$PY" code/00_inspect_env.py > logs/00_inspect_env.log 2>&1 || true
tail -n 40 logs/00_inspect_env.log
"$PY" code/01_prepare_module5.py > logs/01_prepare_module5.log 2>&1
echo PREPARE_EXIT:$?
tail -n 40 logs/01_prepare_module5.log
"$RS" code/02_analyze_module5.R > logs/02_analyze_module5.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 50 logs/02_analyze_module5.log
echo FIG_N:$(ls figures | wc -l)
