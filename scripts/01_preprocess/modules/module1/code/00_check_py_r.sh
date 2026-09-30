#!/bin/bash
set -euo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
echo "=== omicverse ==="
$PY -c 'import pandas; print("pandas", pandas.__version__)'
$PY -c 'import numpy; print("numpy", numpy.__version__)'
$PY -c 'import scipy; print("scipy", scipy.__version__)'
$PY -c 'import openpyxl; print("openpyxl ok")' || echo no_openpyxl
$PY -c 'import sklearn; print("sklearn ok")'
$PY -c 'import matplotlib; print("mpl", matplotlib.__version__)'
$PY -c 'import seaborn; print("sns ok")'
echo "=== r-bio R extra pkgs ==="
/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript -e 'pkgs=c("readxl","ggbeeswarm","survival","survminer","ComplexHeatmap","edgeR","DESeq2","sva"); for(p in pkgs) cat(p, as.character(requireNamespace(p, quietly=TRUE)), "\n")'
echo "=== user R lib ==="
/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript -e 'cat(.libPaths(), sep="\n")'
