#!/usr/bin/env bash
#SBATCH --job-name=zeb1_patient_qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs/patient_qc_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs/patient_qc_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
/opt/bioinfo/envs/omicverse-core/bin/python "$PROJECT/module5/code/scrna_patient_qc.py" \
  --project "$PROJECT" --output "$PROJECT/module5/tables/patient_qc_validation"
