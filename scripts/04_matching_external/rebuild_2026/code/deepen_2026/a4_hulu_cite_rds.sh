#!/bin/bash
set -euo pipefail
DIR=/data-b/liangfuhua/zeb1_a4
mkdir -p "$DIR"
cd "$DIR"
# processed CITE-seq Seurat object only; never RAW.tar / SRA
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271304/suppl/GSE271304_HTSA_CITEseq_seurObj.rds.gz'
nohup wget -c --timeout=30 --tries=8 --retry-connrefused -O GSE271304_HTSA_CITEseq_seurObj.rds.gz "$URL" > cite_rds.log 2>&1 &
echo "pid $! started"
ls -l GSE271304_HTSA_CITEseq_seurObj.rds.gz 2>/dev/null || true
