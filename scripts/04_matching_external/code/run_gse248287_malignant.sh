#!/usr/bin/env bash
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
INPUT="$PROJECT/module5/processed/gse248287_dx"
SCRIPT="$PROJECT/module5/code/gse248287_malignant_pseudobulk.py"
PY=/opt/bioinfo/envs/omicverse-core/bin/python
"$PY" "$SCRIPT" --input "$INPUT" --output "$JOB_OUTPUT/gse248287_malignant"
"$PY" -m pip freeze > "$JOB_OUTPUT/gse248287_malignant/python_packages.txt"
