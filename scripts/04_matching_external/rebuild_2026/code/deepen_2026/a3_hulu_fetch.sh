#!/bin/bash
set -euo pipefail
BASE=/data-b/liangfuhua/zeb1_a3
mkdir -p "$BASE/gse173432" "$BASE/gse165016" "$BASE/gse165207"
echo "HOST=$(hostname) USER=$(whoami)"
df -h /data-b | head
which wget curl Rscript python3 || true

list_suppl() {
  local acc="$1"
  local nnn="GSE${acc:3:3}nnn"
  local url="https://ftp.ncbi.nlm.nih.gov/geo/series/${nnn}/${acc}/suppl/"
  echo "==== $acc $url"
  curl -fsSL --max-time 40 "$url" | sed -n 's/.*href="\([^"]*\)".*/\1/p' | grep -v '^\?' || echo "LIST_FAIL $acc"
}

list_suppl GSE173432
list_suppl GSE165016
list_suppl GSE165207

# smallest first: TPM
TPM_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173432/suppl/GSE173432_CD34_EV_BCL11B_TPM.txt.gz"
echo "GET $TPM_URL"
curl -fL --max-time 120 -o "$BASE/gse173432/GSE173432_CD34_EV_BCL11B_TPM.txt.gz" "$TPM_URL"
ls -lh "$BASE/gse173432"
