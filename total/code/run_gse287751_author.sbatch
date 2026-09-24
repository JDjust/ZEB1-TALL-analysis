#!/usr/bin/env bash
#SBATCH --job-name=zeb1_author_hto
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/logs/author_hto_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module10/logs/author_hto_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
/opt/bioinfo/envs/omicverse-core/bin/python "$PROJECT/module10/code/gse287751_author_pseudobulk.py" \
  --raw "$PROJECT/module10/processed/gse287751_scrna" \
  --metadata "$PROJECT/module10/processed/gse287751_scrna/author_cell_metadata.tsv" \
  --output "$PROJECT/module10/tables/author_condition_pseudobulk"
