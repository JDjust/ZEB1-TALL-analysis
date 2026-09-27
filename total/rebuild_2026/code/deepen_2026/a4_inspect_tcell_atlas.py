"""Inspect CZI T-cell atlas for author stages and paediatric donors. Backed read only."""
from pathlib import Path
import anndata as ad

p = Path(r"D:\_bioinformation\ZEB1\total\rebuild_2026\data\deepen_2026\a4_yayon\thymus_scrna_tcell_subset.h5ad")
print("bytes", p.stat().st_size, flush=True)
adata = ad.read_h5ad(p, backed="r")
print("n_obs", adata.n_obs, "n_vars", adata.n_vars, flush=True)
print("obs_cols", list(adata.obs.columns), flush=True)
obs = adata.obs
for c in obs.columns:
    if any(k in c.lower() for k in ("author", "cell_type", "donor", "assay", "stage", "anno", "cite", "dataset", "study", "age", "development")):
        nunq = int(obs[c].nunique(dropna=True))
        print(f"COL {c}\t{obs[c].dtype}\t{nunq}", flush=True)
        vc = obs[c].astype(str).value_counts().head(40)
        print(vc.to_string(), flush=True)
        print("---", flush=True)
adata.file.close()
