#!/usr/bin/env bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
REMOTE=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4
mkdir -p "$REMOTE/code" "$REMOTE/logs" "$REMOTE/figures" "$REMOTE/tables" "$REMOTE/processed"
cd "$REMOTE"
"$PY" code/01_prepare_module4.py > logs/01_prepare_module4.log 2>&1
echo PREPARE_EXIT:$?
tail -n 40 logs/01_prepare_module4.log
