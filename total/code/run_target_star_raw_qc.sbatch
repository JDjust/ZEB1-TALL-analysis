#!/usr/bin/env bash
#SBATCH --job-name=zeb1-target-qc
#SBATCH --partition=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/logs/target_star_raw_qc_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4/logs/target_star_raw_qc_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module4
RAW=/data-b/liangfuhua/projects/TALL_dataset_download/data/TARGET-ALL-P2/Gene_Expression_Quantification/STAR_-_Counts
python3 "$PROJECT/code/target_star_raw_qc.py" \
  --metadata "$PROJECT/processed/target_samples.tsv" \
  --raw-dir "$RAW" \
  --output "$PROJECT/tables/target_star_raw_qc"
