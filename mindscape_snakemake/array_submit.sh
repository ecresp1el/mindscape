#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-config/config.yaml}"

# Dynamically extract the sample names
if command -v yq &>/dev/null; then
    SAMPLES=( $(yq e '.samples[]' "$CONFIG_FILE") )
else
    SAMPLES=( $(python3 -c "import yaml; print('\n'.join(yaml.safe_load(open('$CONFIG_FILE'))['samples']))") )
fi

NUM_SAMPLES=${#SAMPLES[@]}
echo "🔁 Submitting array job for ${NUM_SAMPLES} samples."

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=seurat_per_sample
#SBATCH --account=parent0
#SBATCH --array=0-$((NUM_SAMPLES-1))
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=12:00:00
#SBATCH --partition=largemem
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/nfs/turbo/umms-parent/MindscapeProjects/logs/array_%A_%a.out
#SBATCH --error=/nfs/turbo/umms-parent/MindscapeProjects/logs/array_%A_%a.err

SAMPLES=( ${SAMPLES[*]} )
sample=\${SAMPLES[\$SLURM_ARRAY_TASK_ID]}

echo "🎯 Task \$SLURM_ARRAY_TASK_ID → sample: \$sample"

snakemake \
  --configfile "$CONFIG_FILE" \
  --config sample="\$sample" \
  /nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/\$sample.rds \
  --use-conda \
  --cores 1
EOF