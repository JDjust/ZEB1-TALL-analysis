#!/usr/bin/env bash
#SBATCH --job-name=zeb1-qc-limma
#SBATCH --partition=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/logs/target_qc_limma_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/logs/target_qc_limma_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4
/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript "$PROJECT/code/target_qc_limma_sensitivity.R"
