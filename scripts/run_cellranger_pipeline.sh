#!/bin/bash

echo "Loading Bioinformatics module..."
module purge
module load Bioinformatics cellranger
module load snakemake

# Define test folder location
TEST_DIR="/nfs/turbo/umms-parent/Manny_test"
OUTPUT_DIR="$TEST_DIR/10496-MW-reanalysis"
CONFIG_DEST="$TEST_DIR/multi_config.csv"

# Copy the original multi config file to your test folder
echo "Copying config.csv to $CONFIG_DEST"
mkdir -p "$TEST_DIR"
cp -f "/nfs/turbo/umms-parent/Carmen_Miranda_scRNAseq /90 Day results/10x_analysis_10496-MW/Sample_10496-MW-Pool01/config.csv" "$CONFIG_DEST"

# Run Snakemake with everything inside the test folder
echo "Running MindScape Cell Ranger workflow with cellranger multi..."
snakemake -j 1 \
  --snakefile mindscape/workflows/cellranger.smk \
  "$OUTPUT_DIR/multi/multiplexing_analysis"

echo "Workflow complete."
