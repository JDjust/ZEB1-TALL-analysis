#!/usr/bin/env bash
set -euo pipefail
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1
cd "$REMOTE"
mkdir -p logs figures
$RS code/04_enrich_figures_code2.R > logs/04_enrich_figures_code2.log 2>&1
echo ENRICH_EXIT:$?
tail -n 30 logs/04_enrich_figures_code2.log
