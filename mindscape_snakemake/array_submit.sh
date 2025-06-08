#!/usr/bin/env bash
set -euo pipefail
set -x  # Enable debug printing

CONFIG_FILE="${1:-config/config.yaml}"

# Load samples
SAMPLES=( $(python3 -c "import yaml; print('\n'.join(yaml.safe_load(open('$CONFIG_FILE'))['samples']))") )
NUM_SAMPLES=${#SAMPLES[@]}
echo "$(date) | 🧬 Samples loaded: ${SAMPLES[*]}"
echo "🔁 Submitting array job for $NUM_SAMPLES samples..."

SAMPLE_STR="${SAMPLES[*]}"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=seurat_per_sample
#SBATCH --account=parent0
#SBATCH --array=0-$((NUM_SAMPLES-1))%2
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=06:00:00
#SBATCH --partition=largemem
#SBATCH --output=/nfs/turbo/umms-parent/MindscapeProjects/logs/array_%A_%a.out
#SBATCH --error=/nfs/turbo/umms-parent/MindscapeProjects/logs/array_%A_%a.err

set -euo pipefail
set -x  # Enable debug printing
samples=($SAMPLE_STR)
sample=\${samples[\$SLURM_ARRAY_TASK_ID]}

echo "$(date) | 🔧 SLURM_ARRAY_TASK_ID = \$SLURM_ARRAY_TASK_ID"
echo "$(date) | 🎯 Running sample: \$sample"

if [[ -z "\$sample" ]]; then
  echo "$(date) | ❌ ERROR: sample is empty!"
  exit 1
fi

# Activate environment
source ~/.bashrc
echo "$(date) | 🔄 Activating conda environment: mindscape-env"
conda activate mindscape-env || { echo "$(date) | ❌ Conda activate failed"; exit 1; }

echo "$(date) | 🔍 Checking for Snakemake..."
which snakemake || { echo "$(date) | ❌ Snakemake not found"; exit 1; }

echo "$(date) | ✅ Conda activated"
which snakemake
snakemake --version

echo "$(date) | 🚀 Running Snakemake for sample: \$sample"
snakemake \
  --configfile "$CONFIG_FILE" \
  --config sample="\$sample" \
  /nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/\$sample.rds \
  --use-conda \
  --forceall \
  --cores 1 \
  --printshellcmds \
  --reason \
  --verbose || {
    echo "$(date) | ❌ Snakemake execution failed for \$sample"
    exit 1
  }
echo "$(date) | ✅ Snakemake finished successfully for \$sample"
EOF