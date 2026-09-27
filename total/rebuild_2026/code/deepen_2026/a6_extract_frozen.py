"""Stream frozen side50 symbols from HTSA_Ghent. Exact symbol match only."""
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import anndata as ad

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026")
out = P / "a6_yayon"
out.mkdir(parents=True, exist_ok=True)
h5 = P / "a4_yayon/thymus_scrna_tcell_subset.h5ad"
frozen = pd.read_csv(P / "gate2_program/frozen_ZEB_side50.tsv", sep="\t")
want_sym = set(frozen["symbol"].astype(str))

adata = ad.read_h5ad(h5, backed="r")
ghent = (adata.obs["study"].astype(str) == "HTSA_Ghent").to_numpy()
sym = adata.var["feature_name"].astype(str).to_numpy()
jmap = {}
for i, s in enumerate(sym):
    if s in want_sym:
        jmap[s] = i
print("recovered", len(jmap), "/", len(want_sym), flush=True)
rec = frozen.loc[frozen.symbol.isin(jmap)].copy()
print(rec.groupby("side").size().to_string(), flush=True)
pd.DataFrame({"symbol": sorted(want_sym - set(jmap)), "status": "absent_exact_symbol"}).to_csv(
    out / "a6_unrecovered_symbols.tsv", sep="\t", index=False
)
meta = adata.obs.loc[ghent, ["donor_id", "cell_type_level_4_explore"]].copy()
adata.file.close()

want = {j: s for s, j in jmap.items()}
ii = np.flatnonzero(ghent)
pos = {int(i): k for k, i in enumerate(ii)}
mat = {s: np.zeros(int(ghent.sum()), dtype=np.float32) for s in jmap}

with h5py.File(h5, "r") as f:
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
                mat[g][k] = float(val)
        n += 1
        if n % 20000 == 0:
            print("rows", n, flush=True)

meta = meta.reset_index(drop=True)
meta["donor_id"] = meta["donor_id"].astype(str)
meta["stage"] = meta["cell_type_level_4_explore"].astype(str)
rows = []
for (don, st), sub in meta.groupby(["donor_id", "stage"], observed=True):
    pos = sub.index.to_numpy(dtype=int)
    if pos.size < 10:
        continue
    rec = {"donor_id": don, "stage": st, "n_cells": int(pos.size)}
    for s, arr in mat.items():
        rec[s] = float(arr[pos].mean())
    rows.append(rec)
pb = pd.DataFrame(rows)
pb.to_csv(out / "a6_donor_stage_pseudobulk.tsv", sep="\t", index=False)
print("pseudobulks", pb.shape, flush=True)
