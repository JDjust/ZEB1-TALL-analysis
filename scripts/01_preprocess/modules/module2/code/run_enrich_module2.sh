#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2
cd "$REMOTE"
mkdir -p logs figures tables processed
"$PY" code/01b_fix_pharmacotype.py > logs/01b_fix_pharmacotype.log 2>&1
echo PHARM_EXIT:$?
tail -n 40 logs/01b_fix_pharmacotype.log
"$RS" code/03_enrich_module2.R > logs/03_enrich_module2.log 2>&1
echo ENRICH_EXIT:$?
tail -n 80 logs/03_enrich_module2.log
echo FIG_N:$(ls figures | wc -l)
# drop the three leftover violin/lollipop panels from 02_analyze
rm -f figures/M2_ZEB1_by_immunophenotype.pdf figures/M2_ZEB1_by_immunophenotype.png
rm -f figures/M2_ZEB1_effect_lollipop.pdf figures/M2_ZEB1_effect_lollipop.png
echo CLEANED_VIOLIN_LOLLIPOP
ls figures | sed -n '1,80p'
