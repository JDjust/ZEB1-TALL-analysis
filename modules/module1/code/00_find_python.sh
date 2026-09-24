#!/bin/bash
set -euo pipefail
for e in omicverse-core pharmaco-core-2026.08 singlecell-py workflow agent-harness r-bio; do
  echo "ENV $e"
  /opt/bioinfo/envs/$e/bin/python -c 'import pandas,numpy,scipy,openpyxl; print("ok", pandas.__version__)' 2>/dev/null || echo fail
done
echo '---sklearn matplotlib seaborn---'
/opt/bioinfo/envs/omicverse-core/bin/python -c 'import sklearn,matplotlib,seaborn; print("omicverse extra ok")' 2>/dev/null || echo omicverse extra fail
/opt/bioinfo/envs/pharmaco-core-2026.08/bin/python -c 'import sklearn,matplotlib,seaborn; print("pharmaco extra ok")' 2>/dev/null || echo pharmaco extra fail
/opt/bioinfo/envs/singlecell-py/bin/python -c 'import sklearn,matplotlib,seaborn,scanpy; print("sc extra ok")' 2>/dev/null || echo sc extra fail
