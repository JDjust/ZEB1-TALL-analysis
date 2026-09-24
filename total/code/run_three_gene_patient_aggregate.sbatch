#!/usr/bin/env bash
#SBATCH --job-name=zeb1_3gene_patients
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs/three_gene_patient_%j.out
#SBATCH --error=/data-b/liangfuhua/projects/ZEB1_TALL_analysis/module5/logs/three_gene_patient_%j.err
set -euo pipefail
PROJECT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
PY=/opt/bioinfo/envs/omicverse-core/bin/python
mkdir -p "$PROJECT/module5/tables/three_gene_patient_validation"
"$PY" "$PROJECT/module5/code/three_gene_patient_aggregate.py" \
  --gse227122 "$PROJECT/module5/processed/GSE227122_annotated.h5ad" \
  --gse248287 "$PROJECT/module5/processed/gse248287_dx" \
  --output "$PROJECT/module5/tables/three_gene_patient_validation/patient_raw_umi_four_genes.tsv"
