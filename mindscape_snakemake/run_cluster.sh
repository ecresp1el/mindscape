#!/bin/bash
# Usage: bash run_cluster.sh [optional_config.yaml] [additional snakemake flags]

set -euo pipefail

if [[ "${1:-}" == *.yaml ]]; then
    CONFIG_FILE="$1"
    shift
else
    CONFIG_FILE="config/config.yaml"
fi

echo "📁 Using config file: $CONFIG_FILE"
mkdir -p /nfs/turbo/umms-parent/MindscapeProjects/logs

# Run Snakemake using the SLURM profile
snakemake \
  --profile profiles/greatlakes \
  --configfile "$CONFIG_FILE" \
  "$@"