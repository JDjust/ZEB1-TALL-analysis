#!/usr/bin/env bash
set -u
PY=/opt/bioinfo/envs/singlecell-py/bin/python
echo "DECOUPLER"
"$PY" -c 'import decoupler; print(decoupler.__version__)' || true
echo "SKLEARN"
"$PY" -c 'import sklearn,scipy; print(sklearn.__version__, scipy.__version__)' || true
echo "260697"
find /data-b/liangfuhua/projects -maxdepth 5 -iname '*260697*' -o -iname '*GSE260697*' 2>/dev/null | head -30
echo "287751"
find /data-b/liangfuhua/projects/ZEB1_TALL_analysis -maxdepth 4 -iname '*287751*' 2>/dev/null | head -30
echo "THYMUS"
find /data-b/liangfuhua/projects/ZEB1_TALL_analysis -maxdepth 3 -type d -iname '*thymus*' 2>/dev/null | head -20
