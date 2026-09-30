# Local working project

The manuscript/data workstation now keeps exactly five directories at `D:/_bioinformation/ZEB1`: `Article`, `Image`, `data`, `script`, and `GitHub`. This implements the requested removal of other items from the daily project root.

- `Article/PROJECT_README.md` is the local project guide; complete Full_Master and both current journal packages remain under Article.
- `data/documentation` holds verification, version selection, applied skills and recovery records. `data/metadata` holds dataset/code/QC indexes and file hashes.
- `script/environment` holds the observed environment files; `script/run_all.ps1` and `script/run_all.sh` run the R entry. Panel entry paths and Image paths remain unchanged.
- `script/skills/zeb1-reproducibility-framework` holds the inspected local skill source; a copy is installed in the actual Codex home. Nuwa and Nature figure skills are already installed.
- Earlier JSON receipts preserve the actual command/path at the time of their execution. The current local acceptance record is `data/documentation/COMPLETION_AUDIT.json`.

This GitHub repository keeps its own root README, environment files, scripts and `.agents` under the supplied publication template. These belong inside the local `GitHub` subdirectory and do not add items to the local ZEB1 root. The data root is still `D:/_bioinformation/ZEB1/data/analysis`; no Hulu download is required.

The current checks validate delivery and frozen-input plotting. Original-read processing, unresolved Harmony convergence and dataset permissions are separate from folder cleanup; see data_verification.md.
