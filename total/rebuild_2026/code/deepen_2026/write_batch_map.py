from pathlib import Path
import pandas as pd

root = Path(r"D:\_bioinformation\ZEB1")
batch = pd.read_csv(root / "data" / "TALL_X01_batch.tsv", header=None, names=["batch"])
cols = pd.read_csv(root / "data" / "TALL_X01_counts.tsv", sep="\t", nrows=0).columns[1:].tolist()
if len(batch) != len(cols):
    raise ValueError("batch rows do not match count columns")
out = pd.DataFrame({"sample_id": cols, "batch": batch["batch"]})
dest = root / "total" / "rebuild_2026" / "data" / "deepen_2026" / "gate2_program" / "tall_x01_sample_batch.tsv"
dest.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(dest, sep="\t", index=False)
print(out.batch.value_counts().to_string())
print("wrote", dest)
