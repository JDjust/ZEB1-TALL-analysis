# Current ZEB1/T-ALL manuscript code candidate

Manuscript: Stage-dependent ZEB1 expression across developmental programs and molecular subtypes in T-cell acute lymphoblastic leukemia.

This directory matches local manuscript revision20 (2026-10-12), based on scientific CURRENT generation196. It is separate from the repository's historical release. Original results and data outside this directory remain historical and must not be mixed with the revised manuscript.

- analysis/run_presubmission_analysis.py implements the subtype-adjusted/within-subtype, donor, fixed-member and cross-cohort analyses. Install analysis/requirements.txt. The entry command is python -X utf8 analysis/run_presubmission_analysis.py --stage all after obtaining the expected processed inputs through lawful source access.
- reproduction_code/render.py, render_layout1.py and render_template.py render the revised six main and six supporting figures from saved source/result tables. Use the separate requirements.rendering.txt environment; its pandas version intentionally differs from the statistical replay environment. The remaining scripts export editable figures/documents or perform bounded QA.
- method_records contains the saved scientific method-script snapshots supporting earlier analyses. These retain their original environment/path requirements. Source-qualified clinical/SU/MC/D results are not automatically refitted by the new-analysis entry command.

This code directory does not contain participant/donor tables, expression matrices, manuscripts, author/ethics forms, credentials, internal logs or the complete local research ZIP. Missing processed inputs mean a clean checkout cannot reproduce all figures/models alone. Local processed-input replay was previously verified; that is distinct from independent raw-data reproduction and from public minimum-data completeness.

See INPUTS_AND_ACCESS.md, FILE_MANIFEST.json and LICENSE_SCOPE_PENDING.md. This is a code-only publication for source inspection. New analysis scripts have no newly assigned reuse license; the repository's limited historical MIT grant is not expanded by this publication. Public minimum observation inputs remain incomplete, so code visibility alone does not establish full reproducibility.
