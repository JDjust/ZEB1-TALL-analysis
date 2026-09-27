"""Write AIEOP counts as samples x genes for clinTALL RNA-only inference."""
from pathlib import Path
import pandas as pd

src = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\clintall\counts_matrix.csv")
out = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\clintall\aieop_gene_expression.tsv")
cnt = pd.read_csv(src, sep=";")
if cnt.columns[0] != "Geneid":
    raise ValueError("unexpected AIEOP count schema")
mat = cnt.set_index("Geneid").T
mat.index.name = "sample_id"
if mat.shape[0] != 120:
    raise ValueError(f"expected 120 samples, got {mat.shape[0]}")
if not pd.api.types.is_integer_dtype(mat.to_numpy().ravel()):
    mat = mat.round().astype("int64")
mat.to_csv(out, sep="\t")
print("wrote", out, "shape", mat.shape)
