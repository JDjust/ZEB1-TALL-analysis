#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3
cd "$REMOTE"
mkdir -p logs figures tables processed
if [[ "${1:-}" != "r-only" ]]; then
  "$PY" code/01d_fix_labels.py > logs/01d_fix_labels.log 2>&1
  echo FIX_EXIT:$?
  tail -n 20 logs/01d_fix_labels.log
fi
rm -f figures/* tables/*
"$RS" code/02_analyze_module3.R > logs/02_analyze_module3.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 40 logs/02_analyze_module3.log
echo FIG_N:$(ls figures | wc -l)
echo TAB_N:$(ls tables | wc -l)
