"""Correct GSE206710 display IDs using the original sample-stage table.

This changes labels only. It asserts row-wise stage/value identity before
editing the two frozen display tables used in Figure 3.
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
source = pd.read_csv(ROOT / "modules/module2/tables/M2_r2_GSE206710_sample_means.tsv", sep="\t")
source = source[source.stage.isin(["DN_early", "DP_immature", "DP_mature"])].reset_index(drop=True)
assert len(source) == 12
for rel in [
    "total/data/validation/zeb_developmental_residual/normal_stage_units_scores.tsv",
    "total/rebuild_2026/data/source_data_rebuilt/F3/F3A_normal_stage_units.tsv",
]:
    path = ROOT / rel
    table = pd.read_csv(path, sep="\t")
    idx = table.dataset.eq("GSE206710")
    assert idx.sum() == 12
    assert table.loc[idx, "stage"].tolist() == source.stage.tolist()
    assert table.loc[idx, "unit"].str.split("|").str[0].tolist() == source.donor.tolist()
    table.loc[idx, "unit"] = source["sample"].str.cat(source.stage, sep="|").to_numpy()
    assert table.loc[idx, "unit"].is_unique
    table.to_csv(path, sep="\t", index=False, float_format="%.16g")
    print(path)
