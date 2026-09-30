#!/usr/bin/env bash
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
DATA=/data-b/liangfuhua/projects/TALL_dataset_download/data/DepMap_26Q1
PY=/opt/bioinfo/envs/omicverse-core/bin/python
"$PY" "$PROJECT/module8/code/depmap26q1_consistent.py" \
  --data "$DATA" \
  --genes "$PROJECT/module8/code/depmap26q1_genes.txt" \
  --output "$JOB_OUTPUT/depmap26q1_consistent"
"$PY" -m pip freeze > "$JOB_OUTPUT/depmap26q1_consistent/python_packages.txt"
