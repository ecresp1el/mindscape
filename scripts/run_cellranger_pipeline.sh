#!/bin/bash

###############################################################################
# Run the MindScape Snakemake pipeline (UMich Great Lakes, module-based)
# - Loads modules
# - Copies and modifies the Cell Ranger multi config
# - Uses a custom reference and Flex probe set file
# - Runs everything in /nfs/turbo/umms-parent/Manny_test
###############################################################################

echo "Loading Bioinformatics module..."
module purge
module load Bioinformatics cellranger
module load snakemake

# Define test folder and reference locations
TEST_DIR="/nfs/turbo/umms-parent/Manny_test"
OUTPUT_DIR="$TEST_DIR/10496-MW-reanalysis"
CONFIG_DEST="$TEST_DIR/multi_config.csv"
SNAKEFILE="/home/elcrespo/Desktop/githubprojects/MindScape/mindscape/workflows/cellranger.smk"

# ✅ Manually downloaded probe-set file for 10x Flex (already placed in this path)
PROBE_FILE="/nfs/turbo/umms-parent/Manny_probe_set/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"

# Create test folder if it doesn't exist
mkdir -p "$TEST_DIR"

# 📄 Copy original config from the core into your reproducible test folder
echo "Copying config.csv to $CONFIG_DEST"
cp -f "/nfs/turbo/umms-parent/Carmen_Miranda_scRNAseq /90 Day results/10x_analysis_10496-MW/Sample_10496-MW-Pool01/config.csv" "$CONFIG_DEST"

# 🧬 Insert 'create-bam,true' after [gene-expression] to support Flex chemistry
echo "Injecting 'create-bam,true' directly after [gene-expression] section..."
awk '
BEGIN { inserted = 0 }
{
  print
  if ($0 ~ /^\[gene-expression\]/) {
    getline nextline
    print "create-bam,true"
    print nextline
    inserted = 1
  }
}
' "$CONFIG_DEST" > "$CONFIG_DEST.tmp" && mv "$CONFIG_DEST.tmp" "$CONFIG_DEST"

# 🧬 Replace reference path to use your reproducible local reference
echo "Replacing reference path with custom reproducible one..."
sed -i 's|^reference,.*|reference,/nfs/turbo/umms-parent/Manny_human_ref/refdata-gex-GRCh38-2020-A|' "$CONFIG_DEST"

# 🧬 Patch the probe-set path to use the manually downloaded Flex file
echo "Injecting probe-set path for Flex chemistry..."
sed -i "s|^probe-set,.*|probe-set,$PROBE_FILE|" "$CONFIG_DEST"

# 🚀 Run Snakemake (from inside TEST_DIR) to launch cellranger multi
echo "Running MindScape Cell Ranger workflow with cellranger multi..."
snakemake -j 1 \
  --snakefile "$SNAKEFILE" \
  --directory "$TEST_DIR" \
  10496-MW-reanalysis

echo "Workflow complete."
