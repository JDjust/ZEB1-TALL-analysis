#!/usr/bin/env python3
"""Rename module figures to Mx_x.y_* so they sort by checklist ID."""
from __future__ import annotations

from pathlib import Path

ROOTS = [
    Path(r"D:\_bioinformation"),
    Path("/data-b/liangfuhua/projects/ZEB1_TALL_analysis"),
]

# extra/unnumbered -> checklist ID. Duplicates of an already-numbered file are deleted.
M2 = {
    "M2_PCA_density": "M2_2.10_PCA_density",
    "M2_ZEB1_stages_vs_TALL": "M2_2.2_stages_vs_TALL",
    "M2_arrow_CD34_to_SP": "M2_2.6_arrow_CD34_to_SP",
    "M2_bar_beeswarm_immunophenotype": "M2_2.12_bar_beeswarm_IP",
    "M2_bidirectional_TALL_vs_stages": "M2_2.16_bidirectional",
    "M2_bump_marker_ranks": "M2_2.2_bump_marker_ranks",
    "M2_circular_heatmap_markers": "M2_2.6_circular_heatmap",
    "M2_complexheatmap_stages": "M2_2.6_complexheatmap",
    "M2_corr_tile_TALL": "M2_2.12_corr_tile",
    "M2_dualaxis_LMO2_ZEB1": "M2_2.4_dualaxis_LMO2_ZEB1",
    "M2_dumbbell_CD34_vs_SP": "M2_2.6_dumbbell_CD34_vs_SP",
    "M2_forest_OLS": "M2_2.16_forest_OLS",
    "M2_heatmap_line_combo": "M2_2.6_heatmap_line_combo",
    "M2_radar_stage_markers": "M2_2.6_radar",
    "M2_radial_nearest_stage": "M2_2.11_radial",
    "M2_scaled_area_trajectories": "M2_2.2_scaled_area",
    "M2_scatter_LMO2_ZEB2_maturity": "M2_2.13_2.14_facet",
    "M2_scatter_ZEB1_vs_LMO2": "M2_2.15_ZEB1_vs_LMO2",
    "M2_scatter_equation_ZEB1_maturity": "M2_2.12_scatter_equation",
    "M2_signed_stage_spearman": "M2_2.8_signed_spearman",
    "M2_stacked_nearest_by_IP": "M2_2.11_stacked",
    "M2_volcano_TALL_vs_normal": "M2_2.8_volcano",
    "M2_waffle_nearest_stage": "M2_2.11_waffle",
}
M3 = {
    "M3_HOXA_vs_rest": "M3_3.9_HOXA_vs_rest",
    "M3_LMO2_LYL1_vs_rest": "M3_3.8_LMO2_LYL1_vs_rest",
    "M3_TAL_vs_rest": "M3_3.10_TAL_vs_rest",
    "M3_mut_FBXW7": "M3_3.11_mut_FBXW7",
    "M3_mut_NOTCH1": "M3_3.11_mut_NOTCH1",
    "M3_mut_PHF6": "M3_3.11_mut_PHF6",
    "M3_mut_PTEN": "M3_3.12_mut_PTEN",
    "M3_mut_RAS_JAK": "M3_3.13_mut_RAS_JAK",
    "M3_pharmacotype_ETP": "M3_3.7_pharmacotype_ETP",
    "M3_pharmacotype_subtype": "M3_3.2_pharmacotype_subtype",
    "M3_subtype_gene_heatmap": "M3_3.2_subtype_heatmap",
}
M4 = {
    "M4_fgsea_h_all_v2023_2_Hs_symbols": "M4_4.6_hallmark_gsea_nes",
}
M7 = {
    "M7_scatter_Dasatinib": "M7_7.7_scatter_Dasatinib",
    "M7_scatter_Mercaptopurine": "M7_7.7_scatter_Mercaptopurine",
    "M7_scatter_Nelarabine": "M7_7.7_scatter_Nelarabine",
    "M7_scatter_Panobinostat": "M7_7.7_scatter_Panobinostat",
    "M7_scatter_Ruxolitinib": "M7_7.7_scatter_Ruxolitinib",
    "M7_scatter_Thioguanine": "M7_7.7_scatter_Thioguanine",
    "M7_scatter_Venetoclax": "M7_7.7_scatter_Venetoclax",
    "M7_scatter_Vorinostat": "M7_7.7_scatter_Vorinostat",
}
MAPS = {"module2": M2, "module3": M3, "module4": M4, "module7": M7}


def apply(root: Path) -> None:
    if not root.exists():
        return
    for mod, mp in MAPS.items():
        fig = root / mod / "figures"
        if not fig.exists():
            continue
        for src_stem, dst_stem in mp.items():
            for ext in (".pdf", ".png"):
                src = fig / f"{src_stem}{ext}"
                dst = fig / f"{dst_stem}{ext}"
                if not src.exists():
                    continue
                if dst.exists() and src.resolve() != dst.resolve():
                    src.unlink()
                    print(f"drop duplicate {src}")
                    continue
                src.rename(dst)
                print(f"{src.name} -> {dst.name}")


def main():
    for r in ROOTS:
        apply(r)


if __name__ == "__main__":
    main()
