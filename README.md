# ZEB1 T-ALL: local-data reproducible analysis and figures

This release replaces multiple presentation versions with one current figure workflow, records dataset/QC provenance and preserves the latest seven main figures and thirteen supplementary figures. It follows the supplied repository template while retaining real analysis module names.

## Run the verified figure workflow

Python requirements are recorded in requirements.txt. R entries use base R to call the original Python renderer; they preserve the current figure data and layout.

```powershell
$env:ZEB1_DATA_ROOT = 'D:/_bioinformation/ZEB1/data/analysis'
$env:ZEB1_PYTHON = 'D:/Python/python.exe'
$env:ZEB1_RSCRIPT = 'D:/R/R-4.6.0/bin/Rscript.exe'
./run_all.ps1 -Figures 2
# One panel:
& $env:ZEB1_RSCRIPT --vanilla scripts/06_figures/Figure/Figure2/A/panelA.R
```

Complete figures are under results/figures; panel entries under scripts/06_figures/Figure and SFigure. The full existing datasets must already be present at ZEB1_DATA_ROOT. Missing inputs cause an error; no remote download is attempted. Use `python src/figure_engine/render_figure.py 2 --out <directory>` to isolate all exports.

## Structure and evidence

- data/: accession index, fixed 100-gene definition, empty sample metadata schema; raw and processed matrices remain local.
- scripts/00_download–07_tables/: real retained module code, stage notes and 92 panel entries. analysis_code_index.tsv maps original to release paths. Upstream preprocessing is retained, not claimed to have been rerun or made portable.
- src/: shared figure engine and R bridge.
- results/: current editable PDF/SVG and PNG figures; large objects and individual-level tables excluded.
- docs/: QC boundaries, observed environment, verification receipts and applied skill methodology.
- manuscript/: pointer to locally retained Full_Master and two current journal packages.

The completed checks cover 20 frozen-data renders, 92 panel deliverables and 64 input hashes. They do not certify original-read QC, Harmony convergence, malignant annotation, independent cohort provenance or a clean-machine environment restore. See docs/data_verification.md and docs/environment_limits.md. Historical local/SSH paths in upstream code require dataset-specific review before rerunning.

The full local Article documents and local data are preserved separately. The existing GitHub visibility is unchanged. License scope is limited; see LICENSE. No manuscript DOI or unverified accession version is invented.
