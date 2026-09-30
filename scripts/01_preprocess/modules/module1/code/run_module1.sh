#!/usr/bin/env bash
set -euo pipefail
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
PY=/opt/bioinfo/envs/omicverse-core/bin/python
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module1
mkdir -p "$REMOTE/processed" "$REMOTE/figures" "$REMOTE/tables" "$REMOTE/logs" "$REMOTE/code"
cd "$REMOTE"
echo "[$(date)] export clinical"
$RS code/01a_export_clinical.R > logs/01a_export_clinical.log 2>&1
echo "[$(date)] prepare"
$PY code/01_prepare_module1.py > logs/01_prepare_module1.log 2>&1
echo "[$(date)] analyze/plot"
$RS code/02_analyze_module1.R > logs/02_analyze_module1.log 2>&1
echo "[$(date)] DONE"
