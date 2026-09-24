#!/usr/bin/env bash
set -uo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
ROOT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis

echo "[$(date)] M1 r2 GSE107011"
mkdir -p "$ROOT/module1/logs" "$ROOT/module1/figures" "$ROOT/module1/tables"
cd "$ROOT/module1"
"$RS" code/06_gse107011.R > logs/06_gse107011.log 2>&1
echo M1_EXIT:$?
tail -n 40 logs/06_gse107011.log
echo M1_FIG:$(ls figures/M1_r2_* 2>/dev/null | wc -l)

echo "[$(date)] M8 r2 DepMap 26Q1"
mkdir -p "$ROOT/module8/logs" "$ROOT/module8/figures" "$ROOT/module8/tables" "$ROOT/module8/processed"
cd "$ROOT/module8"
"$PY" code/01b_prepare_depmap26q1.py > logs/01b_prepare_depmap26q1.log 2>&1
echo M8_PREP_EXIT:$?
tail -n 40 logs/01b_prepare_depmap26q1.log
"$RS" code/02_analyze_module8.R > logs/02_analyze_module8_26q1.log 2>&1
echo M8_FIG_EXIT:$?
tail -n 30 logs/02_analyze_module8_26q1.log

echo "[$(date)] M2 r2 thymus atlas"
mkdir -p "$ROOT/module2/logs" "$ROOT/module2/figures" "$ROOT/module2/tables" "$ROOT/module2/processed"
cd "$ROOT/module2"
"$PY" code/01b_prepare_thymus.py > logs/01b_prepare_thymus.log 2>&1
echo M2_PREP_EXIT:$?
tail -n 60 logs/01b_prepare_thymus.log
"$RS" code/02b_analyze_thymus.R > logs/02b_analyze_thymus.log 2>&1
echo M2_FIG_EXIT:$?
tail -n 40 logs/02b_analyze_thymus.log
echo M2_FIG:$(ls figures/M2_r2_* 2>/dev/null | wc -l)

echo "[$(date)] ALL_REMAINING_CORE_DONE"
