"""A4-C eligibility only. No MAGIC. No fill-in. Does not rescue A/B."""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

P = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a4_yayon")
adata = ad.read_h5ad(P / "merged_thymus_visium_pediatric.h5ad", backed="r")
genes = ["ZEB1", "ZEB2", "CD34", "LYL1", "CD1A"]
var = adata.var
# symbol column
sym_col = None
for c in ["feature_name", "gene_symbols", "symbol", "name"]:
    if c in var.columns:
        sym_col = c
        break
print("var_cols_head", list(var.columns)[:20], flush=True)
if sym_col:
    idx = {str(s): i for i, s in enumerate(var[sym_col].astype(str))}
else:
    idx = {str(s): i for i, s in enumerate(var.index.astype(str))}
print("found", {g: g in idx for g in genes}, flush=True)

tcols = [c for c in adata.obs.columns if str(c).startswith("T_")]
tfrac = adata.obs[tcols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
tot = pd.to_numeric(adata.obs["tot_cell_abundance"], errors="coerce")
t_enr = (tfrac / tot.replace(0, np.nan)) >= 0.5
print("T_enriched_spots", int(t_enr.sum()), "/", len(t_enr), flush=True)

X = adata.X
rows = []
for g in genes:
    if g not in idx:
        rows.append({"gene": g, "eligible": False, "reason": "absent"})
        continue
    j = idx[g]
    # backed sparse: get one column
    col = np.asarray(X[:, j].todense()).ravel() if hasattr(X[:, j], "todense") else np.asarray(X[:, j]).ravel()
    n_ok = 0
    fracs = []
    for lib, m in adata.obs.groupby("library_id", observed=True).groups.items():
        ii = adata.obs.index.get_indexer(m)
        te = t_enr.to_numpy()[ii]
        if te.sum() < 10:
            continue
        nz = (col[ii][te] > 0).mean()
        fracs.append((str(lib), float(nz)))
        if nz >= 0.10:
            n_ok += 1
    elig = n_ok >= 4
    rows.append({
        "gene": g,
        "n_sections_ge10pct": n_ok,
        "n_sections_tested": len(fracs),
        "eligible": elig,
        "reason": "ok" if elig else "sparse_or_below_gate",
    })
    print(g, "sections>=10%", n_ok, "/", len(fracs), "eligible", elig, flush=True)
pd.DataFrame(rows).to_csv(P / "a4c_detection_eligibility.tsv", sep="\t", index=False)
adata.file.close()
