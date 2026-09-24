#!/usr/bin/env bash
#SBATCH --job-name=zeb1_thymus_reference
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/logs/reference_validation_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module2/logs/reference_validation_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
/opt/bioinfo/envs/omicverse-core/bin/python "$PROJECT/module2/code/thymus_reference_aggregate.py" \
 --reference /data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F/Park_HTA/Human_T_cells.h5ad \
 --query "$PROJECT/module5/processed/GSE227122_annotated.h5ad" \
 --output "$PROJECT/module2/tables/reference_validation_${SLURM_JOB_ID}"
