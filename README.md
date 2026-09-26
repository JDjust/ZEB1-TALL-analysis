# ZEB1–ZEB2 axis in T-cell acute lymphoblastic leukemia — analysis and figure code

Code accompanying the manuscript

> **Subtype-associated ZEB1–ZEB2 expression patterns in T-cell acute lymphoblastic leukemia after developmental-expression adjustment**
> Zhou Q\*, Ni Q\*, Liang F, Wang H, Liang X, Li Y, Xu H, Wu C, Li L (\*equal contribution: Zhou Q and Ni Q).

The four lean journal submission packages correspond to tag `v1.0.1-submission`. The numerical sensitivity analyses are unchanged from `v1.0-submission`; the later tag removes duplicate document formats from the local packaging workflow.

The study is a secondary analysis of public, de-identified datasets. This repository contains code and aggregate reader-facing results; no patient-level data are redistributed.

## Repository layout

| Path | Contents |
|---|---|
| `modules/module1` … `modules/module10` | Upstream analysis modules (normal thymus references, Pölönen cohort balance and residual models, independent bulk cohorts, Lim single-cell pseudobulk, HiChIP/scATAC, clinical models). |
| `total/code/` | Targeted validation and sensitivity analyses that produce the analysis tables used by the figures (developmental coordinate, residual models, GSE146901, GSE162280, Lim pseudobulk, chromatin summaries, clinical models). |
| `total/rebuild_2026/code/` | Figure and manuscript pipeline for the submitted version (R/ggplot2). |
| `supplementary_tables/` | Final reader-facing Tables S1–S9, aggregate Figure 7/Table S8 summaries, and revision-four source summaries. |

### Figure pipeline (`total/rebuild_2026/code/`)

| Script | Output |
|---|---|
| `theme_nature.R` | Shared journal theme, palette and export helpers. |
| `build_f1_nature.R` | Figure 1 (normal thymopoiesis). |
| `build_f2_ggplot.R` | Figure 2 (subtype ZEB configurations). |
| `build_f3_revision2.R` | Figure 3 (within-cohort developmental residuals, bootstrap intervals, joint model). |
| `build_f4_ggplot.R` | Figure 4 (BCL11B-rearranged leukemias). |
| `build_main_validation_relayout.R`, `build_f6_gateA.R` | Figures 5–6 (independent bulk RNA and Lim patient-level pseudobulk). |
| `revision3_reviewer_robustness.R`, `build_f7_revision3.R` | Leave-one-subtype-out residuals, alternative developmental scores, adjusted/cross-validated R², fixed and matched BCL11B contrasts, and current Figure 7. |
| `revision4_targeted_sensitivity.R` | Fixed spline df 2/3/4, donor-aware normal-proxy summaries and ridge Cox sensitivity for sparse subtype coefficients. |
| `revision4_matched_program.R` | Pair-blocked BCL11B versus ETP-like limma–voom and Hallmark competitive tests. |
| `build_revision4_supplementary_figures.R` | Final Supplementary Figures S10–S11 in print-size PDF/TIFF/PNG. |
| `build_four_journal_packages.py` | Local editorial assembly of npj Systems Biology and Applications, JCMM, British Journal of Haematology and Scientific Reports variants. |
| `build_s3_audit.R`, `build_s8_clinical_audit.R`, `build_supplement_consolidated.R` | Supplementary Figures S1–S8; the former HiChIP main plate is retained as S9. |
| `revision2_statistical_hardening.R` | 3,000 stratified patient bootstraps and the spline-plus-subtype joint model (seed 20260926). |
| `revision2_clinical_sensitivity.R` | Age-, sex- and white-cell-count-adjusted clinical sensitivity (Supplementary Table S4). |
| `finalize_journal_figures.R` | Runs current figure builders at their native print dimensions (vector PDF, 600-dpi TIFF, PNG); no second geometry or font scaling. |
| `build_reader_tables.py`, `render_supplement_tables_pdf.py` | Reader-facing Supplementary Tables S1–S8 from full-precision derived TSVs; `build_four_journal_packages.py` adds S9. |
| `build_submission.py`, `build_vector_preview.py` | Local manuscript and PDF assembly from the author-maintained Markdown, tables and figures. |

## Reproducing the figures

Requirements: R ≥ 4.4 with `ggplot2`, `patchwork`, `ragg`, `readr`, `ggalluvial`, `ggridges`, `ggrepel`; Python ≥ 3.10 with `pandas`, `python-docx`, `openpyxl`; `pandoc` ≥ 3.

```bash
cd total/rebuild_2026/code
Rscript revision3_reviewer_robustness.R     # requires original 1,309-patient files
Rscript revision4_targeted_sensitivity.R    # spline/proxy and ridge Cox sensitivity
Rscript revision4_matched_program.R         # exploratory matched expression programs
Rscript build_revision4_supplementary_figures.R
Rscript finalize_journal_figures.R          # all current main and supplementary figures
Rscript finalize_journal_figures.R F3 S2    # selected figures
python build_reader_tables.py               # requires derived TSVs
```

The scripts read analysis tables under `total/data/validation/`, `total/rebuild_2026/data/source_data_rebuilt/`, `modules/module*/tables/`, and the original Pölönen count object under `data/polonen_syn54032669/`. Patient-level inputs and patient-level derived tables are not redistributed here. Obtain the source studies under their access terms and regenerate tables with the module and `total/code/` scripts. The current Figure 7 builder and robustness analysis derive the repository root from their script location. Document assembly also requires author-maintained manuscript Markdown and a local Pandoc/Chrome setup; it is not necessary to rerun the numerical analyses.

The `supplementary_tables/` directory provides the submitted S1–S9 Excel workbook, five aggregate S8 TSVs, and revision-four full-precision model and gene-set summaries. Individual-patient residuals and matched-pair identifiers remain outside this public repository. The ridge Cox bootstrap intervals condition on the selected penalty; they do not establish independent clinical prediction.

## Source datasets

| Resource | Accession / source | Role |
|---|---|---|
| Pölönen et al., *Nature* 2024 | St. Jude Cloud / as described in the original publication | Primary diagnostic cohort (1,309 patients, 17 molecular subtypes) |
| Normal human thymus | GSE142522, GSE195812, GSE206710, Park et al. *Science* 2020 (Human Thymus Atlas) | Developmental reference |
| Yang et al., *Nat Commun* 2021 | GSE146901 | Independent ETP/non-ETP bulk RNA; pooled loop calls |
| Di Giacomo et al., *Blood* 2021 | GSE162280 | BCL11B-rearranged case series |
| Xu et al., *Nat Cancer* 2024 | As described in the original publication | Malignant-cell pseudobulk states |
| Gambi et al., *Mol Cell* 2025 | As described in the original publication | HiChIP and scATAC |

## Statistical notes

- Balance is z(ZEB1) − z(ZEB2) on TMM-normalized log2 CPM. The developmental coordinate is z(CD1A) − [z(CD34) + z(LYL1)]/2 and excludes ZEB1, ZEB2 and LMO2.
- Developmental residuals come from `balance ~ ns(coordinate, df = 3)` fitted within the diagnostic cohort; subtype effects from the same spline plus 17 subtypes with sum-to-zero contrasts, HC3 robust tests and Benjamini–Hochberg correction.
- Figure 7 reports same-cohort leave-one-subtype-out fitting, three leave-one-marker-out scores, an expanded score, repeated subtype-stratified five-fold model comparison, and fixed BCL11B comparisons. These checks do not constitute independent replication of the complete 17-subtype ranking. HiChIP is supplementary and exploratory.
- Fully adjusted EFS/OS Cox fits emitted sparse-subtype convergence warnings. Their estimates, P values and FDR remain in Supplementary Tables S1/S4 as audit outputs; the manuscript does not use them to infer survival associations. Supplementary Figure S8 marks these outputs separately from the binary endpoints.
- Single-cell comparisons are made on patient-level pseudobulks, not on cells.

## Contact

Chao Wu (chaowutjmuch@163.com); Limei Li (lilimei116@126.com).
