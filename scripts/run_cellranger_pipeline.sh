#!/bin/bash

###############################################################################
# Run the MindScape Snakemake Cell Ranger workflow (UMich Great Lakes)
#
# This script:
# - Loads the Snakemake module
# - Runs the Snakefile located at mindscape/workflows/cellranger.smk
# - Uses config/config.yaml for input and output paths
#
# HOW TO USE:
#   bash scripts/run_cellranger_pipeline.sh
###############################################################################

echo "Loading Snakemake module..."
module load snakemake/7.32.4

echo "Running MindScape Cell Ranger workflow..."
snakemake -j 1 \
  --snakefile mindscape/workflows/cellranger.smk \
  --configfile config/config.yaml

echo "Workflow complete."
