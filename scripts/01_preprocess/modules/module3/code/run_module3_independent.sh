#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module3
cd "$REMOTE"
mkdir -p logs figures tables processed
"$PY" code/01b_prepare_independent.py > logs/01b_prepare_independent.log 2>&1
echo PREPARE_EXIT:$?
tail -n 40 logs/01b_prepare_independent.log
"$RS" code/02b_analyze_independent.R > logs/02b_analyze_independent.log 2>&1
echo ANALYZE_EXIT:$?
tail -n 40 logs/02b_analyze_independent.log
echo R2_FIG:$(ls figures/M3_r2_* 2>/dev/null | wc -l)
echo TAB:$(ls tables/M3_r2_* 2>/dev/null | wc -l)
