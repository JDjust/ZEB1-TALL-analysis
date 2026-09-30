"""Inspect paediatric Visium obs for CMA and author T-state deconvolution. No gene imputation."""
from pathlib import Path
import anndata as ad

p = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon\merged_thymus_visium_pediatric.h5ad")
out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
adata = ad.read_h5ad(p, backed="r")
print("n_obs", adata.n_obs, "n_vars", adata.n_vars, flush=True)
print("obs_cols", list(adata.obs.columns), flush=True)
print("obsm", list(adata.obsm.keys()) if adata.obsm is not None else None, flush=True)
print("uns_keys", list(adata.uns.keys())[:40], flush=True)
print("layers", list(adata.layers.keys()) if adata.layers is not None else None, flush=True)
obs = adata.obs
obs.head(3).to_csv(out / "visium_pediatric_obs_head.tsv", sep="\t")
# write column summary
rows = []
for c in obs.columns:
    s = obs[c]
    nunq = int(s.nunique(dropna=True))
    rows.append(f"{c}\t{s.dtype}\t{nunq}\t{s.isna().mean():.3f}")
(out / "visium_pediatric_obs_schema.tsv").write_text("column\tdtype\tnunique\tfrac_na\n" + "\n".join(rows), encoding="utf-8")
print("\n".join(rows), flush=True)
for c in obs.columns:
    cl = c.lower()
    if any(k in cl for k in ("cma", "axis", "cortex", "medulla", "cmj", "organ", "l2", "cell2loc", "abundance", "donor", "sample")):
        print("VALUE", c, obs[c].dropna().astype(str).value_counts().head(8).to_dict(), flush=True)
adata.file.close()
