# ZEB1 / T-ALL analysis code

Private author-review code repository. Public release and a code license have not yet been finalized. No DOI has been assigned.

## Contents

- `modules/module*/code/`: existing upstream analysis and plotting scripts recovered from the project modules.
- `total/code/`: manuscript figure generation, targeted sensitivity analyses, evidence audits and document assembly.
- `CODE_MANIFEST.json`: SHA-256 hashes of the uploaded source files.

## Reuse existing results

From `total/`, inspect the workflow with `python code/build_package.py --list-steps`.
The default `python code/build_package.py` redraws figures and rebuilds documents using existing analysis exports. It does not require rerunning upstream analysis. The `--refresh-derived` option is a separate explicit analysis refresh and should not be used for routine figure edits.

## Required inputs and environment

This is a code deposit, not a self-contained raw-data reproduction package. Data tables, author metadata, manuscript files, third-party PDFs and large matrices are not included. Restore the existing local `total/data/` and sibling `modules/module*/tables/` exports before building; `_style.py` also supports the bundled `total/data/module_tables/` fallback. Missing inputs will prevent execution.

Python figure/document scripts use packages including matplotlib, numpy, pandas, scipy and python-docx; individual analysis scripts have additional requirements. R and Slurm scripts retain their original package and server-path requirements. Paths such as local Windows locations and Hulu work directories must be adapted to another system. Word pagination checks require Windows and Microsoft Word. This deposit does not claim clean-machine execution or invent an environment lockfile.

## Interpretation and access

Source and dataset identifiers are retained in the scripts. Access and redistribution terms must be checked with each original data resource. Private GitHub access is not public data/code availability. The repository does not include patient-level data, raw sequencing, DepMap matrices, credentials or author contact metadata.

No new statistical analyses were run to create this deposit. Existing results and their limitations remain in the local author-review package.
