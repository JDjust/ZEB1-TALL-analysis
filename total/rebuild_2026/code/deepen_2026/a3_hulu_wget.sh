#!/bin/bash
set -eu
cd /data-b/liangfuhua/zeb1_a3
mkdir -p gse165016 gse165207
fetch() {
  dest="$1"; url="$2"
  if [ -s "$dest" ] && gzip -t "$dest" 2>/dev/null; then
    echo "ok $dest"
    return
  fi
  echo "WGET $url"
  wget -q --timeout=60 --tries=3 -O "$dest.tmp" "$url" && mv "$dest.tmp" "$dest"
  gzip -t "$dest"
  echo "wrote $dest $(wc -c < "$dest")"
}
fetch gse165016/GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5265nnn/GSM5265325/suppl/GSM5265325_BCL11B_E22952_peaks.narrowPeak.gz
fetch gse165016/GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5265nnn/GSM5265326/suppl/GSM5265326_BCL11B_E14259_peaks.narrowPeak.gz
for gsm in GSM5028224 GSM5028225 GSM5028226 GSM5028227 GSM5028228 GSM5028229 GSM5028230 GSM5028231 GSM5028232; do
  echo "list $gsm"
  curl -fsSL --max-time 40 "https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/${gsm}/suppl/" | sed -n 's/.*href="\([^"]*bedpe[^"]*\)".*/\1/p'
done
echo HULU_LIST_DONE
