from pathlib import Path
import anndata as ad

out = Path(r"D:/_bioinformation/ZEB1/data/analysis\total\rebuild_2026\data\deepen_2026\a4_yayon")
adata = ad.read_h5ad(out / "thymus_scrna_tcell_subset.h5ad", backed="r")
obs = adata.obs
keep = [
    "donor_id", "donor_age", "age_group", "study", "assay", "author_cell_type",
    "cell_state", "cell_type_level_0", "cell_type_level_1", "cell_type_level_2",
    "cell_type_level_3", "cell_type_level_4_explore", "cell_type", "sample",
    "enrichment", "is_primary_data",
]
keep = [c for c in keep if c in obs.columns]
sub = obs[keep].copy()
sub.to_csv(out / "tcell_atlas_obs_sample.tsv", sep="\t", index=False)
for c in keep:
    vc = sub[c].astype(str).value_counts()
    (out / f"tcell_atlas_{c}_counts.tsv").write_text(
        vc.rename_axis(c).reset_index(name="n").to_csv(sep="\t", index=False),
        encoding="utf-8",
    )
    print(c, int(vc.shape[0]), "top", vc.index[0], int(vc.iloc[0]), flush=True)
adata.file.close()
