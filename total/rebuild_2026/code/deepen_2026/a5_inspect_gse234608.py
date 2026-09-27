from pathlib import Path
import gzip

p = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a5_gse234608\GSE234608_raw_counts.csv.gz")
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
    header = f.readline().rstrip("\n")
    cols = header.split(",")
    print("ncols", len(cols))
    print("col0", cols[0])
    print("samples", cols[1:])
    print("n_samples", len(cols) - 1)
    n = 0
    hits = []
    want = {"ZEB1", "ZEB2", "CD1A", "CD34", "LYL1", "BCL11B"}
    for line in f:
        n += 1
        g = line.split(",", 1)[0].strip().strip('"')
        g0 = g.split(".")[0]
        if g.upper() in want or g0 in {
            "ENSG00000148516", "ENSG00000169554", "ENSG00000158477",
            "ENSG00000174059", "ENSG00000104903", "ENSG00000127152",
        }:
            hits.append(g)
    print("nrows", n)
    print("hits", hits)
