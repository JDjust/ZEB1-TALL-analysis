# Developmentally contextualized ZEB1–ZEB2 states in T-ALL

Code and aggregate reader-facing outputs for the article by Zhou Q, Ni Q, Liang F, Wang H, Liang X, Li Y, Xu H, Wu C and Li L. Zhou and Ni contributed equally; Wu and Li are corresponding authors. Fuhua Liang is the third author.

The current visual-redesign snapshot is `v2.0.2-visual-redesign`. It contains the seven main and ten supplementary figure builders used for the current Word/PDF package. The numerical analyses come from the frozen A6/A7/A8 evidence locks; this release changes their presentation and manuscript organization, not the underlying scientific questions.

## Scientific scope

The paper tests whether molecular T-ALL subtypes have ZEB1–ZEB2-associated expression states after conditioning on a three-gene developmental-expression **proxy**. The primary shared residual is fitted within 1,309 pediatric/young-adult diagnostic patients. Frozen 50-gene programs on each residual side are compared with normal thymus, developmentally matched BCL11B and ETP-like patients, independent leukemia cohorts, normal marrow, malignant-cell pseudobulks and an adult T-ALL cohort.

The analysis does **not** infer cell of origin, direct ZEB regulation, a replicated chromatin mechanism or clinical utility. No external dataset fully reproduces the 17-subtype primary ranking. Acute BCL11B overexpression opposes the observed lesion-associated state; sparse-subtype Cox estimates are not used for independent survival inference.

## Current figure and table map

| Item | Role |
|---|---|
| Figure 1 | Normal thymus stage and spatial scaffold |
| Figure 2 | 1,309-patient subtype landscape and within-cohort residual |
| Figure 3 | Frozen programs versus normal developmental change |
| Figure 4 | Developmentally matched BCL11B pole and lesion-defined convergence |
| Figure 5 | Frozen program preservation across independent leukemia contexts |
| Figure 6 | Normal marrow affinity and passive-admixture audit |
| Figure 7 | Adult T-ALL directional architecture, complete cluster context and 100-gene integration |
| Figures S1–S10 | Eligibility, sensitivity, case, assay and negative-boundary audits; one page each |
| Tables S1–S10 | Dataset manifest, normal/stage units, subtype statistics, classifier audit, frozen genes, matched cases, external programs, marrow/specimen, malignant/adult and perturbation/chromatin/clinical audits |

The exact panel-to-source map is [`total/rebuild_2026/PANEL_SOURCE_MAP.md`](total/rebuild_2026/PANEL_SOURCE_MAP.md). The workbook is [`supplementary_tables/Supplementary_Tables_S1-S10.xlsx`](supplementary_tables/Supplementary_Tables_S1-S10.xlsx). Each numbered table may span multiple worksheet tabs; full-precision TSVs are regenerated locally from the original sources and are not all redistributed here.

## Code layout

| Path | Content |
|---|---|
| `modules/module1` … `modules/module10` | Legacy upstream dataset preparation and analyses retained for provenance |
| `total/code/` | Primary cohort, external and sensitivity analysis scripts |
| `total/rebuild_2026/code/` | Revision-two to revision-four models, A6/A7/A8 targeted analyses, final plate builders and submission assembly |
| `supplementary_tables/` | Aggregate reader-facing workbook and non-identifying result summaries |

The current plate and package builders are:

```text
total/rebuild_2026/code/editorial_upgrade_1.py
total/rebuild_2026/code/editorial_upgrade_2.py
total/rebuild_2026/code/editorial_upgrade_347.py
total/rebuild_2026/code/editorial_upgrade_5.py
total/rebuild_2026/code/editorial_upgrade_6.py
total/rebuild_2026/code/editorial_upgrade_supplement.py
total/rebuild_2026/code/editorial_upgrade_supplement_rest.py
total/rebuild_2026/code/editorial_upgrade_10.py
total/rebuild_2026/code/editorial_rebuild_tables.py
total/rebuild_2026/code/build_editorial_journal_packages.py
total/rebuild_2026/code/build_chinese_advisor_docs.py
```

The Python plate builders use Matplotlib at the final 180 mm physical width and preserve editable PDF text. They read frozen derived TSVs under `total/rebuild_2026/data/`, which are excluded from this repository when they contain patient-level or source-controlled data. The packaging script also reads author-maintained manuscript text and local figure files under `submission/`; it cannot run from this public repository alone. The manifest and script paths document provenance, while users must obtain the original public data under each source study's terms to regenerate patient-level intermediates.

## Reproduction sequence

1. Obtain the source datasets and follow the upstream module and `total/code/` scripts to generate the patient-, donor- and gene-level derived TSVs locally.
2. Run the frozen robustness/targeted scripts under `total/rebuild_2026/code/` (including `revision3_reviewer_robustness.R`, `revision4_targeted_sensitivity.R` and `revision4_matched_program.R`) and the A6/A7/A8 scripts documented in the panel map.
3. Generate the final reader-facing workbook with `editorial_rebuild_tables.py`, then render the plates with the corresponding `editorial_upgrade_*.py` scripts in the panel map order.
4. Use `build_editorial_journal_packages.py` with the author-maintained manuscript source to assemble journal-specific Word/PDF files.

These instructions describe the execution order and source dependencies. A clean-room end-to-end rerun from source archives has **not** been independently certified; readers should inspect eligibility manifests, script warnings and locked audit notes before interpreting a regenerated result.

## Public source accessions and roles

| Resource | Role |
|---|---|
| GSE142522, GSE195812, GSE206710, Park/HTA, Yayon human thymus atlas | Normal developmental stage and spatial context |
| Pölönen et al. diagnostic T-ALL | Primary 1,309-patient, 17-subtype discovery cohort |
| AIEOP120 official-classifier series | Moderate external pediatric architecture audit |
| GSE146901, GSE234608 | Independent bulk expression context |
| GSE243914 | Blast-enriched T-ALL context |
| GSE248287 | Source-author malignant-cell pseudobulk context |
| GSE280250 | Author adult T-ALL series, 79 diagnostic cases |
| GSE253355 | Normal marrow donor-by-lineage context |
| GSE162280 | Twelve BCL11B lesion-defined cases |
| GSE165209 | Acute BCL11B-overexpression boundary |
| Lim single-cell resource, HiChIP, pooled loops, scATAC | Supporting patient-state and chromatin-boundary analyses |

The manuscript and Table S1 specify source papers, biological units, eligibility and whether each dataset is independent, overlapping or contextual. The GSE280250 adult cohort was absent from ALLCatchR2's listed development cohorts; complete patient non-overlap has not been proven.

## Contact

Chao Wu: chaowutjmuch@163.com. Limei Li: lilimei116@126.com.
