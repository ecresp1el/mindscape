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

echo "$(date) | 🧾 Submitting with job array range 0-$((NUM_SAMPLES-1))"
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=seurat_per_sample
#SBATCH --account=parent0
#SBATCH --array=0-$((NUM_SAMPLES-1))%6
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

echo "$(date) | 🧪 CONFIG_FILE = $CONFIG_FILE"
echo "$(date) | 🧪 Sample string: \${samples[*]}"
echo "$(date) | 🧪 Selected sample: \$sample"

echo "$(date) | 🔧 SLURM_ARRAY_TASK_ID = \$SLURM_ARRAY_TASK_ID"

echo "$(date) | 🎯 Running sample: \$sample"
echo "$(date) | 📂 Present working directory: \$PWD"
echo "$(date) | 📜 Listing contents of current dir:"
ls -lh

echo "👤 Whoami: \$(whoami)"
echo "🏠 HOME: \$HOME"
echo "🐍 Snakemake at: \$(which snakemake)"

if [[ -z "\$sample" ]]; then
  echo "$(date) | ❌ ERROR: sample is empty!"
  exit 1
fi

# Activate environment
export BASHRCSOURCED=1
source ~/.bashrc
echo "$(date) | 🔄 Activating conda environment: mindscape-env"
conda activate mindscape-env || { echo "$(date) | ❌ Conda activate failed"; exit 1; }

echo "$(date) | 🔍 Checking for Snakemake..."
which snakemake || { echo "$(date) | ❌ Snakemake not found"; exit 1; }

echo "$(date) | ✅ Conda activated"
which snakemake
snakemake --version

echo "$(date) | 🧪 Snakemake command:"
echo snakemake --configfile "$CONFIG_FILE" --config sample="\$sample" --use-conda --forceall --cores 1 --printshellcmds --verbose /nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/\$sample.rds

echo "$(date) | 🚀 Running Snakemake for sample: \$sample"
export PYTHONUNBUFFERED=1
snakemake \
  /nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/\$sample.rds \
  --configfile "$CONFIG_FILE" \
  --config sample="\$sample" \
  --use-conda \
  --forceall \
  --cores 1 \
  --printshellcmds \
  --verbose || {
    echo "$(date) | ❌ Snakemake execution failed for \$sample"
    exit 1
  }
echo "$(date) | ✅ Snakemake finished successfully for \$sample"
EOF