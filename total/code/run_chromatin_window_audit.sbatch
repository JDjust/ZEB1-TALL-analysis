#!/usr/bin/env bash
#SBATCH --job-name=zeb1_window_audit
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6/logs/window_audit_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module6/logs/window_audit_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
/opt/bioinfo/envs/omicverse-core/bin/python "$PROJECT/module6/code/chromatin_window_audit.py" \
 --project "$PROJECT" --output "$PROJECT/module6/tables/window_audit_${SLURM_JOB_ID}"
