"""
Collect the exact source tables used by each main figure into total/data/,
renamed as Source Data for the manuscript. This makes the figures fully
reproducible: every panel maps to one file here.
"""
import os
import shutil
import _style as S

MAP = {
    # Figure 1
    "SourceData_Fig1A_disease_ranking.tsv": ("module1", "M1_1.6_ZEB1_disease_ranking.tsv"),
    "SourceData_Fig1B_TALL_vs_comparators.tsv": ("module1", "M1_1.2_1.5_pairwise_vs_TALL.tsv"),
    "SourceData_Fig1C_cohort_effects.tsv": ("module1", "M1_1.13_cohort_effect_sizes.tsv"),
    "SourceData_Fig1C_meta.tsv": ("module1", "M1_1.14_random_effects_meta.tsv"),
    "SourceData_Fig1D_ZEB_axis.tsv": ("module1", "M1_1.7_ZEB1_ZEB2_LMO2_effects.tsv"),
    "SourceData_Fig1E_ROC.tsv": ("module1", "M1_1.9_multiROC_three_genes.tsv"),
    "SourceData_Fig1F_immune_atlas.tsv": ("module1", "M1_r2_1.1_celltype_ZEB1_rank.tsv"),
    # Figure 2
    "SourceData_Fig2A_thymus_trajectory.tsv": ("module2", "M2_stage_keygene_rho.tsv"),
    "SourceData_Fig2B_scRNA_stage.tsv": ("module2", "M2_r2_GSE195812_stage_gene_means.tsv"),
    "SourceData_Fig2C_OLS.tsv": ("module2", "M2_2.16_2.17_regressions.tsv"),
    "SourceData_Fig2D_replication.tsv": ("module2", "M2_2.18_replication.tsv"),
    "SourceData_Fig2E_nearest_stage.tsv": ("module2", "M2_2.11_nearest_stage.tsv"),
    "SourceData_Fig2F_scRNA_stage_coupling.tsv": ("module2", "M2_r2_GSE195812_gene_vs_stage.tsv"),
    # Figure 3
    "SourceData_Fig3A_subtype_counts.tsv": ("module3", "M3_3.1_subtype_counts.tsv"),
    "SourceData_Fig3BC_subtype_effects.tsv": ("module3", "M3_3.2_3.4_subtype_effects.tsv"),
    "SourceData_Fig3D_effect_matrix.tsv": ("module3", "M3_r2_3.5_effect_matrix.tsv"),
    "SourceData_Fig3E_meta.tsv": ("module3", "M3_r2_3.5_meta.tsv"),
    "SourceData_Fig3F_ETP_model.tsv": ("module3", "M3_pharmacotype_ETP_age_model.tsv"),
    # Figure 4
    "SourceData_Fig4A_keygene_rho.tsv": ("module4", "M4_keygene_stats.tsv"),
    "SourceData_Fig4B_keygene_adjusted.tsv": ("module4", "M4_r2_keygenes_adjusted.tsv"),
    "SourceData_Fig4C_genomewide_adjusted.tsv": ("module4", "M4_r2_subtype_adjusted_ZEB1.tsv"),
    "SourceData_Fig4D_signatures_TARGET.tsv": ("module4", "M4_4.9_4.19_signature_vs_ZEB1.tsv"),
    "SourceData_Fig4D_signatures_Pharmacotype.tsv": ("module4", "M4_4.32_pharmacotype_signature_vs_ZEB1.tsv"),
    "SourceData_Fig4E_hallmark_gsea.tsv": ("module4", "M4_4.6_hallmark_fgsea.tsv"),
    # Figure 5
    "SourceData_Fig5A_composition.tsv": ("module5", "M5_5.13_composition.tsv"),
    "SourceData_Fig5B_discovery_assoc.tsv": ("module5", "M5_5.41_patient_ZEB1_associations.tsv"),
    "SourceData_Fig5DE_stage_mapping.tsv": ("module5", "M5_r2_map_patient.tsv"),
    # Figure 6
    "SourceData_Fig6B_methylation.tsv": ("module6", "M6_6.17_ZEB1_promoter_beta_summary.tsv"),
    "SourceData_Fig6C_occupancy.tsv": ("module6", "M6_r2_GSE154675_promoter_signal.tsv"),
    "SourceData_Fig6D_TLX1_KD.tsv": ("module9", "M9_keygene_effects.tsv"),
    "SourceData_Fig6E_Lmo2_TG.tsv": ("module9", "M9_gse188225_keygenes.tsv"),
    # Figure 7
    "SourceData_Fig7A_drug_spearman.tsv": ("module7", "M7_7.7_7.8_drug_spearman.tsv"),
    "SourceData_Fig7B_MRD.tsv": ("module7", "M7_7.18_ZEB1_MRD.tsv"),
}

n = 0
for dest, (mod, name) in MAP.items():
    src = S.tbl(mod, name)
    if os.path.exists(src):
        shutil.copyfile(src, os.path.join(S.DATA_DIR, dest))
        n += 1
    else:
        print("MISSING:", src)
print("copied %d/%d source-data tables to total/data/" % (n, len(MAP)))
locked = os.path.join(S.DATA_DIR, "depmap26q1_job1702")
for dest, source in [
    ("SourceData_Fig7C_ZEB1_dependency.tsv", "tall_zeb1_dependency.tsv"),
    ("SourceData_Fig7DE_drug_target_GE.tsv", "target_gene_effect_summary.tsv"),
]:
    src = os.path.join(locked, source)
    if not os.path.exists(src):
        raise FileNotFoundError("DepMap 26Q1 locked output missing: " + src)
    shutil.copyfile(src, os.path.join(S.DATA_DIR, dest))
print("Supplementary Source Data are written by compute_supp_analyses.py (FigS2/S4/S6/S10–S12).")
