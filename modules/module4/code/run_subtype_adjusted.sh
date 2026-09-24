#!/usr/bin/env bash
set -euo pipefail
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
M4=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4
cd "$M4"
mkdir -p logs figures tables
"$RS" code/04_subtype_adjusted.R > logs/04_subtype_adjusted.log 2>&1
echo ADJ_EXIT:$?
tail -n 50 logs/04_subtype_adjusted.log
echo R2_FIG:$(ls figures/M4_r2_* 2>/dev/null | wc -l)
