# ZEB1–ZEB2 axis in T-cell acute lymphoblastic leukemia — analysis and figure code

Code accompanying the manuscript

> **Molecular subtypes of T-cell acute lymphoblastic leukemia differentially reconfigure a developmentally patterned ZEB1–ZEB2 axis**
> Zhou Q\*, Ni Q\*, Wang H, Liang X, Li Y, Xu H, Wu C, Li L (\*equal contribution).

The study is a secondary analysis of public, de-identified datasets. This repository contains code only; no patient-level data are redistributed.

## Repository layout

| Path | Contents |
|---|---|
| `modules/module1` … `modules/module10` | Upstream analysis modules (normal thymus references, Pölönen cohort balance and residual models, independent bulk cohorts, Lim single-cell pseudobulk, HiChIP/scATAC, clinical models). |
| `total/code/` | Targeted validation and sensitivity analyses that produce the analysis tables used by the figures (developmental coordinate, residual models, GSE146901, GSE162280, Lim pseudobulk, chromatin summaries, clinical models). |
| `total/rebuild_2026/code/` | Figure and manuscript pipeline for the submitted version (R/ggplot2). |

### Figure pipeline (`total/rebuild_2026/code/`)

| Script | Output |
|---|---|
| `theme_nature.R` | Shared journal theme, palette and export helpers. |
| `build_f1_nature.R` | Figure 1 (normal thymopoiesis). |
| `build_f2_ggplot.R` | Figure 2 (subtype ZEB configurations). |
| `build_f3_revision2.R` | Figure 3 (within-cohort developmental residuals, bootstrap intervals, joint model). |
| `build_f4_ggplot.R` | Figure 4 (BCL11B-rearranged leukemias). |
| `build_f5_ggplot.R`, `build_f6_ggplot.R`, `build_s8_ggplot.R`, `build_s9_ggplot.R`, `build_s12_ggplot.R` + `build_main_validation_relayout.R` | Figures 5–7 (independent bulk RNA, Lim single-cell states, HiChIP). |
| `build_s1_ggplot.R` … `build_s12_ggplot.R` + `build_supplement_consolidated.R` | Supplementary Figures S1–S7. |
| `revision2_statistical_hardening.R` | 3,000 stratified patient bootstraps and the spline-plus-subtype joint model (seed 20260926). |
| `revision2_clinical_sensitivity.R` | Age-, sex- and white-cell-count-adjusted clinical sensitivity (Supplementary Table S4). |
| `finalize_journal_figures.R` | Runs every figure builder and exports the 180-mm journal versions (PDF with editable text, 600-dpi TIFF, PNG). |
| `build_submission.py` | Assembles manuscript, Supplementary Information and supplementary-table files. |

## Reproducing the figures

Requirements: R ≥ 4.4 with `ggplot2`, `patchwork`, `ragg`, `readr`, `ggalluvial`, `ggridges`, `ggrepel`; Python ≥ 3.10 with `pandas`, `python-docx`, `openpyxl`; `pandoc` ≥ 3.

```bash
cd total/rebuild_2026/code
Rscript finalize_journal_figures.R          # all main and supplementary figures
Rscript finalize_journal_figures.R F3 S2    # selected figures
python build_submission.py                  # documents and tables
```

The figure scripts read the analysis tables under `total/data/validation/` and `modules/module*/tables/`. These tables are derived from the public datasets below and are not redistributed here; regenerate them with the module and `total/code/` scripts after obtaining the source data under each resource's access terms. Paths are resolved relative to the repository root.

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
- Single-cell comparisons are made on patient-level pseudobulks, not on cells.

## Contact

Chao Wu (chaowutjmuch@163.com); Limei Li (lilimei116@126.com).
