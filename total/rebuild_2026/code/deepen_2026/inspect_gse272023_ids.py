from pathlib import Path
import pandas as pd

p = Path(r"D:\_bioinformation\ZEB1\total\data\validation\GSE272023_MatrixRNAseqVSTcountsFilted.txt.gz")
mat = pd.read_csv(p, sep="\t")
print("shape", mat.shape)
print("columns_head", list(mat.columns[:4]))
print("columns_tail", list(mat.columns[-6:]))
print("id_col", mat.columns[0], "example", mat.iloc[:3, 0].tolist())
need = ["ENSG00000148516", "ENSG00000169554", "ENSG00000158477",
        "ENSG00000174059", "ENSG00000104903", "ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"]
ids = mat.iloc[:, 0].astype(str)
for g in need:
    hit = ids[ids.str.contains(g, regex=False)]
    print(g, "n=", len(hit), "values=", hit.head(5).tolist())
