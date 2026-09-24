#!/usr/bin/env python3
import gzip
p="/data-b/liangfuhua/projects/TALL_dataset_download/data/GSE13159/GPL570.annot.gz"
with gzip.open(p,"rt",encoding="utf-8",errors="replace") as f:
    for i,line in enumerate(f):
        print(i, line[:300].rstrip())
        if i>40:
            break
