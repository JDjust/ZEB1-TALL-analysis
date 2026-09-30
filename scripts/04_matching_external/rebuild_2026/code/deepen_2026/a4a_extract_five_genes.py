"""Stream 5 genes from CSR without loading the 4.5GB index array."""
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import anndata as ad

P = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
h5 = P / "thymus_scrna_tcell_subset.h5ad"
adata = ad.read_h5ad(h5, backed="r")
obs = adata.obs
ghent = (obs["study"].astype(str) == "HTSA_Ghent").to_numpy()
sym = adata.var["feature_name"].astype(str).to_numpy()
genes = ["ZEB1", "ZEB2", "CD1A", "CD34", "LYL1"]
jmap = {g: int(np.where(sym == g)[0][0]) for g in genes}
print("jmap", jmap, flush=True)
meta = obs.loc[ghent, ["donor_id", "cell_type_level_4_explore"]].copy()
adata.file.close()

want = {j: g for g, j in jmap.items()}
out = {g: np.zeros(int(ghent.sum()), dtype=np.float32) for g in genes}
ii = np.flatnonzero(ghent)
pos = {int(i): k for k, i in enumerate(ii)}

with h5py.File(h5, "r") as f:
    print("X keys", list(f["X"].keys()) if "X" in f else list(f.keys()), flush=True)
    grp = f["X"]
    indptr = grp["indptr"][:]
    indices = grp["indices"]
    data = grp["data"]
    n = 0
    for i in ii:
        s, e = int(indptr[i]), int(indptr[i + 1])
        idx = indices[s:e]
        dat = data[s:e]
        k = pos[int(i)]
        for j, val in zip(idx, dat):
            g = want.get(int(j))
            if g is not None:
                out[g][k] = float(val)
        n += 1
        if n % 20000 == 0:
            print("rows", n, flush=True)

meta = meta.reset_index(drop=True)
for g in genes:
    meta[g] = out[g]
meta.to_csv(P / "a4a_ghent_five_genes.tsv", sep="\t", index=False)
print("wrote", len(meta), "nonzero", {g: int((out[g] > 0).sum()) for g in genes}, flush=True)
