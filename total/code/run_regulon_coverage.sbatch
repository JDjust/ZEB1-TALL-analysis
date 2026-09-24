#!/usr/bin/env bash
#SBATCH --job-name=zeb1_regulon_audit
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module9/logs/regulon_audit_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module9/logs/regulon_audit_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
/opt/bioinfo/envs/omicverse-core/bin/python "$PROJECT/module9/code/regulon_coverage_audit.py" \
 --project "$PROJECT" \
 --network /data-b/liangfuhua/projects/TALL_dataset_download/data/external_fig467F/CollecTRI_omnipath.tsv \
 --output "$PROJECT/module9/tables/regulon_audit_${SLURM_JOB_ID}"
