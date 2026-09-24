#!/bin/bash
set -euo pipefail
DATA=/data-b/liangfuhua/projects/TALL_dataset_download/data
echo "=== GSE13159 gsm characteristics sample ==="
grep -E "Sample_characteristics|Sample_title|Sample_source|leukemia|diagnosis|disease" "$DATA/GSE13159/GSE13159_gsm_brief.txt" | head -n 40
echo "=== unique characteristic keys ==="
grep "^!Sample_characteristics" "$DATA/GSE13159/GSE13159_gsm_brief.txt" | sed 's/:.*//' | sort | uniq -c | sort -nr | head -n 30
echo "=== series matrix first 60 lines ==="
zcat "$DATA/GSE13159/GSE13159_series_matrix.txt.gz" | sed -n '1,60p'
echo "=== TARGET clinical files ==="
find "$DATA/TARGET-ALL-P2/Clinical_Supplement" -type f | head -n 30
echo "=== TARGET biospecimen files ==="
find "$DATA/TARGET-ALL-P2/Biospecimen_Supplement" -type f | head -n 20
echo "=== TARGET RNA n files ==="
find "$DATA/TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts" -type f | wc -l
echo "=== TARGET RNA first file ==="
f=$(find "$DATA/TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts" -type f | head -n 1)
echo "FILE=$f"
head -n 8 "$f"
echo "=== PanALL rlog colnames/nrow ==="
python3 - <<'PY'
import gzip
p="/data-b/liangfuhua/projects/TALL_dataset_download/data/StJude_PanALL/panALL-1988S-rlog.subtype.txt.gz"
with gzip.open(p,"rt") as f:
    header=f.readline().rstrip("\n").split("\t")
    print("n_cols", len(header))
    print("first_cols", header[:15])
    n=0
    genes=[]
    for i,line in enumerate(f):
        n+=1
        g=line.split("\t",1)[0]
        if g.upper() in {"ZEB1","ZEB2","LMO2"} or "ZEB1" in g.upper() or "LMO2" in g.upper():
            genes.append(g)
        if i<2:
            print("row", i, "gene", g, "nfields", len(line.split("\t")))
print("n_rows", n)
print("hit_genes", genes[:20])
PY
echo "=== R packages ==="
Rscript -e 'pkgs=c("ggplot2","ggpubr","pheatmap","pROC","limma","metafor","GEOquery","affy","hgu133plus2.db"); for(p in pkgs) cat(p, ": ", requireNamespace(p, quietly=TRUE), "\n")'
echo "=== python ==="
python3 -c "import pandas,numpy,scipy,matplotlib,sklearn; print('pyok', pandas.__version__)"
echo "=== which R ==="
which Rscript; Rscript --version
