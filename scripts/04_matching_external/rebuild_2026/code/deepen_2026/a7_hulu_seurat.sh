#!/bin/bash
set -euo pipefail
DIR=/data-b/liangfuhua/zeb1_a7
mkdir -p "$DIR"
cd "$DIR"
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253355/suppl/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz'
# full hematopoietic atlas only; never MSC subset / SRA
nohup wget -c --timeout=30 --tries=8 --retry-connrefused -O GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz "$URL" > seurat.log 2>&1 &
echo "pid $! started"
ls -l GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds.gz 2>/dev/null || true
