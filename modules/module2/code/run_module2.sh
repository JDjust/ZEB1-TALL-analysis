#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2
cd "$REMOTE"
mkdir -p logs figures tables processed
if [[ "${1:-}" != "r-only" ]]; then
  "$PY" code/01_prepare_module2.py > logs/01_prepare_module2.log 2>&1
  echo PREPARE_EXIT:$?
  tail -n 50 logs/01_prepare_module2.log
fi
rm -f figures/* tables/*
"$RS" code/02_analyze_module2.R > logs/02_analyze_module2.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 40 logs/02_analyze_module2.log
echo FIG_N:$(ls figures | wc -l)
echo TAB_N:$(ls tables | wc -l)
