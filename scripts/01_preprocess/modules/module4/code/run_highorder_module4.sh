#!/usr/bin/env bash
set -euo pipefail
GMT=/data-b/liangfuhua/projects/TALL_dataset_download/data/MSigDB
M4=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4
mkdir -p "$GMT" "$M4/logs" "$M4/figures" "$M4/tables"
cd "$GMT"
base=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs
for f in h.all.v2023.2.Hs.symbols.gmt \
         c2.cp.reactome.v2023.2.Hs.symbols.gmt \
         c5.go.bp.v2023.2.Hs.symbols.gmt \
         c3.tft.v2023.2.Hs.symbols.gmt
do
  if [ ! -s "$f" ]; then
    echo FETCH "$f"
    curl -fsSL -o "$f" "$base/$f" || echo FAIL "$f"
  fi
done
ls -lh "$GMT"
/opt/bioinfo/envs/omicverse-core/bin/python "$M4/code/rename_figures.py"
/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript "$M4/code/03_highorder_module4.R" > "$M4/logs/03_highorder_module4.log" 2>&1
echo HIGHORDER_EXIT:$?
tail -n 60 "$M4/logs/03_highorder_module4.log"
echo FIG_N:$(ls "$M4/figures" | wc -l)
