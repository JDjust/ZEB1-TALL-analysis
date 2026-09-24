#!/usr/bin/env bash
#SBATCH --job-name=zeb1_label_stability
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=00:45:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs/harmony_sensitivity_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs/harmony_sensitivity_%j.err
set -euo pipefail
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
/opt/bioinfo/envs/omicverse-core/bin/python "$PROJECT/module5/code/harmony_label_sensitivity.py" \
 --input "$PROJECT/module5/processed/GSE227122_annotated.h5ad" \
 --output "$PROJECT/module5/tables/harmony_label_sensitivity_${SLURM_JOB_ID}"
