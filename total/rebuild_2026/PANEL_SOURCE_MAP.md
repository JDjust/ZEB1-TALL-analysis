# Final editorial figure source map

This map describes the **current exported plates**, not older figure drafts. Frozen analysis tables live under `total/rebuild_2026/data/`; `D/` below means `deepen_2026/`, and `SD/` means `source_data_rebuilt/`. Figures live in `submission/article/00_Figure_Master/{main,supplementary}/`. All drawing scripts live together in `total/rebuild_2026/code/` and render directly at 180 mm width. No scientific thresholds or gene memberships were chosen during drawing.

| Plate | Literal panel content | Frozen inputs / drawing script |
|---|---|---|
| F1A | 6-resource × 10-window coverage matrix | `SD/F1/`, `D/a4_yayon/`; `editorial_upgrade_1.py` |
| F1B | Four source-specific normal trajectories | `SD/F1/F1D_stage_tracks.tsv`; `editorial_upgrade_1.py` |
| F1C,D | All 33 Yayon states, five-donor paired inset | `D/a4_yayon/a4a_*`; `editorial_upgrade_1.py` |
| F1E,F | All 24 spatially mapped T states; 19 exact RNA/CMA joins | `D/a4_yayon/a4b_*`, `a4ab_*`; `editorial_upgrade_1.py` |
| F2A | 17-subtype integrated measured-median matrix | `SD/F2/`, `SD/F3/`, `D/a7_lineage/a7_0_patient_scores.tsv`; `editorial_upgrade_2.py` |
| F2B,C | 1,309 patient $B$ versus $D$; model variance lollipop | `SD/F2/`, `SD/F3/F3B_within_cohort_curve.tsv`; `editorial_upgrade_2.py` |
| F2D,E | Residual half-violins with patients; 17 bootstrap median intervals | `SD/F3/F3B-C_patient_level.tsv`, `F3D-F_subtype_statistics.tsv`; `editorial_upgrade_2.py` |
| F2F | Six comparable discovery–AIEOP subtype slopes; small-class anchor | `D/aieop120_official/`; `editorial_upgrade_2.py` |
| F3A,B | partial-$R^2$ rank curves with small count flow; genome-wide effects | `D/gate2_program/`; `editorial_upgrade_347.py` |
| F3C | Frozen 100 genes × ten discovery residual deciles | `D/gate2_program/`, `D/a7_lineage/a7_0_patient_scores.tsv`; `editorial_upgrade_347.py` |
| F3D,E | Five-donor normal inset; 90-gene normal–leukemia effect map | `D/a6_yayon/`; `editorial_upgrade_347.py` |
| F3F | All 48 evaluable matched ZEB2-side gene effects | `D/a6_yayon/`, `SD/revision4_targeted/`; `editorial_upgrade_347.py` |
| F4A,B | Fixed 18-pair developmental and residual comparisons | `SD/revision3_robustness/bcl11b_etp_matched_*`; `editorial_upgrade_347.py` |
| F4C | Two genuine Venn diagrams: 39/50 and 29/50 intersections | `D/gate2_program/`, `D/gate2_sensitivity/`; `editorial_upgrade_347.py` |
| F4D | Pair-blocked gene volcano and Hallmark CAMERA bubble plot | `SD/revision4_targeted/`; `editorial_upgrade_347.py` |
| F4E,F | Twelve lesion-case dumbbells; seven adult BCL11B patient matrix | `SD/F4/F4C_GSE162280_cases.tsv`, `D/a8_adult/`; `editorial_upgrade_347.py` |
| F5A,B | Five-cohort validation matrix; three-cohort score effect matrix | `D/a65_external/`, `D/a7b_malignant/`, `D/a8_adult/`; `editorial_upgrade_5.py` |
| F5C,D,E | Gene-preservation bubble forest; 15 malignant patients; 94-gene matrix | same frozen families; `editorial_upgrade_5.py` |
| F6A,B,C | Donor × lineage matrices; paired differences; 93-gene affinity | `D/a7_lineage/scores/`; `editorial_upgrade_6.py` |
| F6D,E | BM/PB coefficient audit; measured 48-gene × four-cohort decomposition | `D/a7_lineage/`, `D/a7c_decompose/`; `editorial_upgrade_6.py` |
| F7A–D | Adult $B$/$D$ spline, both score relations with age inset, gene effects | `D/a8_adult/`; `editorial_upgrade_347.py` |
| F7E,F | All adult HC cluster patients; 100-gene × eight-context effect matrix | `D/a8_adult/`, `D/a6_yayon/`, `D/a7_lineage/`, external cohort effects; `editorial_upgrade_347.py` |

| Supplement | Literal panel inventory | Drawing script |
|---|---|---|
| S1 | A coverage; B 33-state aligned matrix; C marker-direction matrix; D 16-library T-enriched Visium detection matrix | `editorial_upgrade_supplement_rest.py` |
| S2 | A cross-fit; B 17×5 rank heatmap; C model $R^2$; D fixed spline df; E batch coefficient audit; F normal range | `editorial_upgrade_supplement_rest.py` |
| S3 | A official calls; B 33-case fusion→official-class alluvial; C slopegraph; D assignment probabilities | `editorial_upgrade_supplement.py` |
| S4 | A genome effects; B binary original/held-out gene membership; C CAMERA; D PC overadjustment hexbin | `editorial_upgrade_supplement_rest.py` |
| S5 | A 18-pair bipartite match; B residual pairs; C Hallmark bubbles; D complete 12-case matrix | `editorial_upgrade_supplement_rest.py` |
| S6 | A six-vector correlation/shared-gene matrix; B–D three constituent-effect comparisons | `editorial_upgrade_supplement.py` |
| S7 | A donor×lineage score matrix; B 93-gene expression; C donor differences; D cell gate | `editorial_upgrade_supplement_rest.py` |
| S8 | A 15-patient matrix; B separate Z1/Z2 scatter; C cell counts; D 41 state pairs; E 12 Day-28 pairs | `editorial_upgrade_supplement_rest.py` |
| S9 | A 79-case cluster→immature-call→residual-bin alluvial; B cluster counts; C age; D immature-call residuals | `editorial_upgrade_supplement.py` |
| S10 | A three-donor acute OE; B six-track ZEB locus peak plots; C ChIP peak-hit/assay availability matrix; D six-endpoint attenuation forest | `editorial_upgrade_10.py` |

**Audit rules.** `S1D_Visium_library_detection.tsv` was recalculated on the original T-enriched spot mask and reproduces the frozen 2/16 ZEB2 gate. Figure 7F external contrast columns are sign-oriented and then scaled within each column for display; heatmap color magnitudes must not be compared across assays. S9 residual bins are visualization categories only. The archived HiChIP 3-versus-4 exact P=0.057, pooled-loop and scATAC boundaries remain in Table S10 and the Supplementary text; S10C does not recast their readouts as ChIP peak counts.

Build order: run `editorial_upgrade_1.py`, `_2.py`, `_347.py`, `_5.py`, `_6.py`, `editorial_upgrade_supplement.py`, `editorial_upgrade_supplement_rest.py`, `_10.py`; then `sync_editorial_visual_text.py`, `sync_editorial_master_refs.py`, `build_editorial_journal_packages.py`, `build_chinese_advisor_docs.py`. The plot builder scripts assert key source dimensions before export. Open every PNG at final size and check literal types, text, whitespace and clipping before submission.
