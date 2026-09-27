# Editorial rebuild: locked panel source map

This map implements `D:/_bioinformation/ZEB1/revison.txt`. It is the single working specification for the seven new main figures. All paths below are relative to `D:/_bioinformation/ZEB1/total/rebuild_2026/data/`. `SD/` means `source_data_rebuilt/`; `D/` means `deepen_2026/`. The scientific boundaries in `D/FINAL_EVIDENCE_AUDIT.md` and the A6/A7/A8 locks govern wording. The old figure numbers are not evidence of final placement.

## Global figure contract

- Final-size, two-column plates: 180 mm wide; heights F1 190, F2 205, F3 210, F4 195, F5 200, F6 185, F7 215 mm. Export PDF with editable text; rasterize only dense points or tissue images. Do not rescale fonts or geometry in a finalizer.
- Negative residual and ZEB2-side use muted blue; positive residual and ZEB1-side use burnt orange. BCL11B, ETP-like and TLX3 use one fixed highlight scheme throughout; all other subtypes are grey unless the panel needs their full categorical palette. Heatmaps use a zero-centred blue–white–orange scale, with no pink fill.
- Each plate has one visual centre and 3–5 supporting panels. Panel headers are short labels; the legend explains the question, data source, unit, n, statistical method, uncertainty and boundary for every letter.
- No panel may show a number solely copied from this plan. The cited TSV or lock must contain it. A caption must distinguish within-cohort residuals from external cohort scores and distinguish donor, patient, case and gene denominators.

## Main figures

| Panel | Exact source | Graphic and required labels | Boundary / QA |
|---|---|---|---|
| F1A | `SD/F1/F1D_stage_tracks.tsv`; `D/a4_yayon/a4a_stage_balance.tsv` | Clean thymic entry → cortex/DP → selection → SP/medulla schematic anchored to measured stage names | Diagram is anatomical orientation, not inferred cell lineage or pseudotime |
| F1B | `SD/F1/F1D_stage_tracks.tsv`; `SD/F1/F1_cross_atlas_window_summary.tsv` | Four aligned stage trajectories, one shared ZEB balance scale; mark source and biological unit | GSE142522, GSE195812, GSE206710, Park/HTA are separate references, not 20 independent donors |
| F1C | `D/a4_yayon/a4a_stage_balance.tsv`; `D/a4_yayon/a4a_donor_stage_pseudobulk.tsv`; `D/a4_yayon/a4a_gate_lock.tsv` | Eight prespecified early/DP/SP poles drawn in biological order from 33 eligible author states; stage-group medians early −0.37, DP +0.52, SP −0.20 | The source table has 33 eligible states and is sorted by balance, not stage order; never connect its raw rows as a developmental trajectory. ETP has one donor. Full 33-state audit belongs in S1. |
| F1D | `D/a4_yayon/a4a_donor_level_balance.tsv`; `D/a4_yayon/a4a_gate_lock.tsv` | Five donor paired early–DP slope lines, 4/5 same direction | Wilcoxon P=0.063 is a small sensitivity label, not significance marker |
| F1E | `D/a4_yayon/a4b_tstate_cma_by_donor.tsv`; `D/a4_yayon/a4b_tstate_cma_localization.tsv` | Cortex–medulla CMA distributions by author state; six spatial donors | Cell2location/CMA localization, not direct ZEB2 Visium detection; see `a4c_detection_eligibility.tsv` |
| F1F | `D/a4_yayon/a4ab_stage_balance_and_cma.tsv` | Integrated position-versus-stage-balance topology, with early/DP/SP tracks | Visual centre; show spatial association, not directional migration |
| F2A | `SD/F2/F2B_subtype_counts.tsv`; `SD/F2/F2C_developmental_medians.tsv` | Compact 17-subtype strip with n and median developmental proxy; 1,309 total | No pie/Sankey; subtype order fixed through D/E |
| F2B | `SD/F2/F2C-D_patients.tsv`; `SD/F3/F3B_within_cohort_curve.tsv` | All-patient B versus D scatter and cohort-internal natural spline; highlight BCL11B, ETP-like, TLX3 | B=z(ZEB1)−z(ZEB2); D=z(CD1A)−[z(CD34)+z(LYL1)]/2; D is expression proxy |
| F2C | `SD/F2/F2F_model_variance.tsv`; `SD/revision3_robustness/model_comparison.tsv` | Raw, adjusted and CV R² for immunophenotype, D, subtype, D+subtype; emphasize incremental 13.9 percentage points in raw R² | Label model metric beside each bar; do not imply nested models from raw R² alone |
| F2D | `SD/F3/F3B-C_patient_level.tsv`; `SD/F3/F3D-F_subtype_statistics.tsv` | 17-subtype residual raincloud/beeswarm in fixed order; individual patients visible | Show n for small groups; residual is from common leukemia-cohort spline |
| F2E | `SD/F3/F3D-F_subtype_statistics.tsv`; `SD/revision2_statistics/subtype_effects_bootstrap.tsv` | Bootstrap median forest, all 17 rows; BCL11B −1.76, ETP-like −0.23, TLX3 +0.64 highlighted | Median residual and joint-model subtype coefficient are different estimands |
| F2F | `D/aieop120_official/aieop120_residual_with_official_labels.tsv`; `D/aieop120_official/aieop120_subtype_residual_medians.tsv`; `D/aieop120_official/aieop120_validation_summary.tsv` | External architecture strip: BCL11B 4, ETP-like 24, TLX3 25 and six classes with n≥5 | 120 samples, 16 official classes; moderate external support, not a 17-class rank replication |
| F3A | `D/gate2_program/freeze_summary.tsv`; `D/gate2_program/frozen_ZEB_side50.tsv` | Frozen selection rule as compact quantitative protocol: gene model, FDR<.01, |βR|≥.25, top 50 partial R² each side | Freeze precedes validation; source models and exact gene set in table S5 |
| F3B | `D/gate2_program/all_genes_residual_effects.tsv`; `D/gate2_program/frozen_ZEB_side50.tsv` | Genome-wide βR versus partial R², with 100 frozen genes emphasized and sparse gene labels | Avoid DEG-style P-value volcano |
| F3C | `D/a7_lineage/a7_0_patient_scores.tsv` joined by `sample_id` to `SD/F3/F3B-C_patient_level.tsv` | Two frozen-score-versus-residual continuous plots, with subtype density rug | Existing locked 1,309-patient scores; join cardinality must remain one-to-one, with no reselection |
| F3D | `D/a6_yayon/a6_donor_stage_scores.tsv`; `D/a6_yayon/a6_1_donor_program.tsv` | Early→DP→SP normal trajectories for both frozen scores, donor points and stage summaries | Normal program scores are donor×stage pseudobulk, not leukemia patient values |
| F3E | `D/a6_yayon/a6_2_gene_effects.tsv`; `D/a6_yayon/a6_2_lock.tsv` | **Hero:** 90-gene Δnormal versus βR quadrant; 13 concordant, 32 discordant, 45 weak, Spearman ρ=−.50 | Neutral threshold follows A6 lock; rho descriptive; no mechanistic verb |
| F3F | `D/a6_yayon/a6_3_cameraPR.tsv`; `D/a6_yayon/a6_3_lock.tsv` | Two frozen ZEB2 subsets in matched BCL11B–ETP contrast: 38 discordant/neutral P=3.93×10⁻¹⁶; 10 early-concordant P=.019 | CAMERA gene-set evidence; do not call either subset causal |
| F4A | `SD/revision3_robustness/bcl11b_etp_matched_pairs.tsv`; `SD/revision3_robustness/bcl11b_etp_matched_summary.tsv` | 18 pairs: paired developmental-coordinate dumbbell, small absolute distance | Matching is internal, not external replication |
| F4B | Same pair files as F4A | Paired residual plot, BCL11B−ETP-like median contrast −1.92 with bootstrap interval from summary | Explain paired statistic and whether CI is bootstrap |
| F4C | `D/gate2_sensitivity/heldout_matched_cameraPR.tsv`; `D/gate2_sensitivity/heldout_vs_original_overlap.tsv` | Leave-BCL11B/ETP-out frozen-program CAMERA effect in held-out pairs | Held-out within-cohort phenotype bridge, not independent validation |
| F4D | `SD/revision4_targeted/matched_BCL11B_ETP_gene_results.tsv`; `SD/revision4_targeted/matched_BCL11B_ETP_hallmark_camera.tsv` | Compact top effect genes plus Hallmark context; pair-blocked model | IFN/JAK–STAT means transcriptional context, not a proven pathway mechanism |
| F4E | `SD/F4/F4C_GSE162280_cases.tsv`; `SD/F4/F4E_lesion_group_ranges.tsv` | **Hero:** 12-case ZEB1/ZEB2 CPM dumbbell, partner/phenotype annotation; 12/12 ratio<0; 7/7 non-ZEB2 partners | Mixed AML/MPAL/ETP case series; no ordinary T-ALL comparator; cannot infer sufficiency |
| F4F | `D/a8_adult/a8_2_merged.tsv`; `D/a8_adult/A82_CLOSED.md` | Adult high-confidence BCL11B n=7 residual strip: median −2.34, all seven negative, age 31–75, 0/7 ETP-like | ALLCatchR2 context, not second adult generalization analysis |
| F5A | `D/FINAL_EVIDENCE_AUDIT.md`; `D/a65_external/a65_summary.tsv`; `D/a7b_malignant/a7b_summary.tsv`; `D/a8_adult/a8_1_summary.tsv` | Evidence strip: independent bulk → ETP/MPAL bulk → blast-enriched → author malignant → adult | Identify independent units; do not claim zero patient overlap with ALLCatchR2 training |
| F5B | `D/gate3_validation/gse146901_program_scores.tsv`; `D/gate3_validation/gse146901_program_tests.tsv` | Frozen ZEB1-side/ZEB2-side score and balance in GSE146901 using fixed visual template | Healthy peripheral T comparator limitation in legend |
| F5C | `D/a5_gse234608/a5_patient_frozen_scores.tsv`; `D/a5_gse234608/a5_pregate_tests.tsv` | Same frozen-score template in GSE234608 | Keep ETP/MPAL specimen definitions exact |
| F5D | `D/a7b_malignant/a7b05_sample_scores.tsv`; `D/a7b_malignant/a7b05_summary.tsv` | Same template in blast-enriched GSE243914; ETP 6 versus non-ETP 18 | Blast-enriched, not malignant-pure |
| F5E | `D/a65_external/a65_GSE146901_genes.tsv`; `D/a65_external/a65_GSE234608_genes.tsv`; `D/a7b_malignant/a7b05_genes.tsv`; `D/a7b_malignant/a7b_genes.tsv`; `D/a8_adult/a8_genes.tsv` | **Hero:** expected-direction fraction and gene-level rho for five cohorts (87/.60, 80/.64, 91.8/.77, 77.7/.62, 98.7/.86) | Denominator of recovered/tested genes must appear per cohort; directional structure, not each gene independently replicated |
| F5F | `D/a7b_malignant/a7b_patient_scores.tsv`; `D/a7b_malignant/a7b_bootstrap.tsv` | Author-defined malignant 15 patient pseudobulks, two score-versus-B scatter plots: Z1 ρ+.73, Z2 ρ−.21 | Z2 CI includes zero; show it, no invented CNV/TCR labels |
| F6A | `D/a7_lineage/scores/a7a_unit_scores_min20.tsv`; `D/a7_lineage/scores/a7a_summary_min20.tsv` | Normal marrow four broad-lineage paired-score landscape: HSPC, CLP, T, myeloid | Donor×lineage pseudobulk; use locked author label map |
| F6B | `D/a7_lineage/scores/a7a_donor_min20.tsv`; `D/a7_lineage/A7_A_LOCK.md` | Donor-paired ZEB2-side differences versus T: myeloid 12/12, HSPC 11/11, CLP 0/8 | Sample units and paired n visible; do not infer origin |
| F6C | `D/a7_lineage/scores/a7a_genes_min20.tsv`; `D/a7_lineage/A7_A_LOCK.md` | Gene Δmyeloid−T versus −βR: 75% expected, ρ=.60 | ZEB1-side HSPC caveat in legend/S7 |
| F6D | `D/a7_lineage/a7_0_frozen_genes_specimen_adjusted.tsv`; `D/a7_lineage/a7_0_overall.tsv` | Primary versus specimen-adjusted β, identity line; 100/100 same sign, ρ=.9995 | BM 928 vs blood 381; adjustment alone does not prove malignant intrinsic |
| F6E | F6A–D plus F5D/F, `D/a7c_decompose/A7C_CLOSED.md` | Evidence-boundary ladder with quantitative anchors: lineage affinity → specimen audit → blast-enriched → author malignant | This is a synthesis of measured results; never say contamination is completely excluded |
| F7A | `D/a8_adult/a8_patient_scores.tsv`; `D/a8_adult/a8_1_summary.tsv` | Adult cohort n=79 balance B versus D and its **own** ns(D,3) curve | Refit B/D/R within GSE280250; seven cases <18 in author Adult series |
| F7B | Same adult patient/summary files | ZEB1-side score versus adult R; ρ+.60, bootstrap CI .43–.73 | Patient is the independent unit |
| F7C | Same adult patient/summary files | ZEB2-side score versus adult R; ρ−.66, bootstrap CI −.76 to −.52 | Do not reuse pediatric residual values |
| F7D | `D/a8_adult/a8_genes.tsv`; `D/a8_adult/a8_recovered.tsv` | **Hero:** βadult versus βdiscovery; 86/100 recovered, 78/86 tested, 98.7% expected, ρ=.86 | Label recovered and tested denominators separately |
| F7E | `D/a8_adult/a8_1_summary.tsv`; `D/a8_adult/A81_LOCK.md` | Compact age sensitivity, all 79 / age≥18 72 / age≥21 70 | All 79 is primary; adult author label includes seven <18 |
| F7F | `D/a8_adult/a8_2_merged.tsv`; `D/a8_adult/a8_2_subtype_map.tsv`; `D/a8_adult/A82_CLOSED.md` | Only prespecified high-confidence BCL11B, ETP-like, TAL1 DP-like, TLX3 immature context | Only three unique mapped n≥5 classes contribute ordering; do not headline ρ=1.0 |
| F7G | Main F1–F7 locked results plus `D/FINAL_EVIDENCE_AUDIT.md` | Evidence relationship schematic: normal scaffold → within-cohort residual → non-equivalent side programs → external/adult support | Explicit boundaries: no cell of origin, no sole admixture explanation, no acute BCL11B-OE induction, no replicated chromatin mechanism |

## Supplementary renumbering and one-page scope

| Final figure | One-page question | Principal source families |
|---|---|---|
| S1 | Normal harmonization, biological units, proxy marker eligibility, spatial detection limit | `SD/Supplement_S1/`, `D/a4_yayon/` |
| S2 | Developmental-proxy/residual model robustness | `SD/revision3_robustness/`, `SD/revision4_targeted/spline_df_2_3_4.tsv`, `D/gate2_sensitivity/` |
| S3 | AIEOP official classifier and external architecture audit | `D/aieop120_official/`, `D/gate1_classifier/` |
| S4 | Frozen program selection and held-out internal robustness | `D/gate2_program/`, `D/gate2_sensitivity/` |
| S5 | BCL11B matched transcriptome and complete 12-case lesion audit | `SD/revision4_targeted/`, `SD/F4/`, `SD/S11/` |
| S6 | Independent bulk and blast-enriched frozen-program audit | `D/gate3_validation/`, `D/a5_gse234608/`, `D/a65_external/`, `D/a7b_malignant/a7b05_*` |
| S7 | Normal marrow lineage, specimen sensitivity and 45+3 decomposition | `D/a7_lineage/`, `D/a7c_decompose/` |
| S8 | Author malignant-cell pseudobulk and Lim denominators/Day28 null | `D/a7b_malignant/`, `SD/S7/`, `SD/Supplement_S3/` |
| S9 | Adult age and ALLCatchR2 subtype context | `D/a8_adult/` |
| S10 | Acute BCL11B-OE and chromatin boundaries | `D/a3_gse165209/`, `SD/S10/`; clinical attenuation stays in S10 table |

The historical bulk meta-analysis loses its standalone figure. Old chromatin S5 and S9 merge into S10. Clinical attenuation and complete raw subtype expression distributions move to tables/source data. All S1–S10 are **single pages**; no `page1/page2` exports.

## Supplementary workbook (one file, ten numbered tables)

| Table | Required contents and source |
|---|---|
| S1 | Dataset accession, population, unit, role and overlap audit; `SD/supplementary_tables/Table_S2_dataset_units.tsv` plus audit locks |
| S2 | Normal stage harmonization, donor/stage units, spatial state mapping and detection eligibility; `SD/F1/`, `D/a4_yayon/` |
| S3 | 17-subtype residual, adjusted effects, CI, HC3 and FDR; `SD/F3/F3D-F_subtype_statistics.tsv` |
| S4 | AIEOP official calls, 33/33 S9 fusion check and subtype QC; `D/aieop120_official/` |
| S5 | Frozen 100 genes, selection stats, held-out derivation and A6 classification; `D/gate2_program/`, `D/gate2_sensitivity/`, `D/a6_yayon/` |
| S6 | 18 matched pairs, pair-blocked transcriptome context and 12 lesion cases; `SD/revision3_robustness/`, `SD/revision4_targeted/`, `SD/F4/` |
| S7 | External cohort score and constituent-gene validation with recovered/tested denominators; F5 source family |
| S8 | Marrow donor×lineage, specimen adjustment and 45+3 decomposition; `D/a7_lineage/`, `D/a7c_decompose/` |
| S9 | GSE248287 author malignant patient pseudobulk, adult A8-1, ALLCatchR2 context; `D/a7b_malignant/`, `D/a8_adult/` |
| S10 | Acute OE, chromatin positive/negative context, exploratory clinical model warnings and attenuation; `D/a3_gse165209/`, `SD/S10/`, `SD/S6/` |

## Items to resolve before declaring submission-ready

1. Verify F3C one-to-one join of the existing 1,309-patient frozen side-score table to the residual table; do not select new genes.
2. For F5E and F7D, print recovered and model-tested gene counts separately. The A8 lock gives 86 recovered but only 78/86 limma-tested; the displayed 98.7% applies to the tested set.
3. Verify the spatial panel uses CMA/cell2location state localization; direct ZEB2 Visium detection is not evaluable.
4. Resolve all four journal-specific abstract, word, display-item and Methods limits against current official instructions before exporting versions. A master with long Methods is not a submission file for BJH/JCMM.
5. Render all 17 plates at final physical size and check text/axes/panel alignment; detailed letter-by-letter legends and Word-embedded main figures are required in all journal versions.
