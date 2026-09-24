#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6
mkdir -p "$REMOTE/logs" "$REMOTE/figures" "$REMOTE/tables" "$REMOTE/processed"
cd "$REMOTE"
"$PY" code/01_prepare_module6.py > logs/01_prepare_module6.log 2>&1
echo PREPARE_EXIT:$?
tail -n 40 logs/01_prepare_module6.log
if [ -f code/02_analyze_module6.R ]; then
  "$RS" code/02_analyze_module6.R > logs/02_analyze_module6.log 2>&1
  echo ANALYZE_EXIT:$?
  tail -n 30 logs/02_analyze_module6.log
fi
echo FIG_N:$(ls figures 2>/dev/null | wc -l)
