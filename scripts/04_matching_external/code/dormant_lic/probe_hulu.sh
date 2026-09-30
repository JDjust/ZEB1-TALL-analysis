#!/usr/bin/env bash
set -u
echo "HOST $(hostname)"
for p in \
  /opt/bioinfo/envs/singlecell-py/bin/python \
  /opt/bioinfo/envs/omicverse-core/bin/python \
  /opt/bioinfo/envs/infercnvpy/bin/python
do
  echo "PY $p"
  "$p" -c 'import scanpy,anndata; print("scanpy", scanpy.__version__)' || true
done
echo "RSEURAT"
/opt/bioinfo/envs/r-bio/bin/Rscript -e 'cat(as.character(packageVersion("Seurat")), "\n")' || true
echo "TAR"
tar -tvf /data-b/liangfuhua/projects/TALL_dataset_download/data/GSE260948/suppl/GSE260948_RAW.tar
echo "COLLECTRI"
find /data-b/liangfuhua/projects/ZEB1_TALL_analysis -iname '*collectri*' -o -iname '*dorothea*' 2>/dev/null | head -20
