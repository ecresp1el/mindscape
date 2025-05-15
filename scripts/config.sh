#!/usr/bin/env bash
#
# mindscape/scripts/config.sh
#

# Where to run the test
TEST_DIR="/nfs/turbo/umms-parent/Manny_test"

# Original multi‐config CSV
TURBO_CONFIG_SOURCE="/nfs/turbo/umms-parent/Carmen_Miranda_scRNAseq /90 Day results/10x_analysis_10496-MW/Sample_10496-MW-Pool01/config.csv"

# Flex probe set (on Turbo storage)
PROBE_PATH="/nfs/turbo/umms-parent/Manny_probe_set/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"

# Where your Snakemake workflow lives (relative)
SNAKEFILE="mindscape/workflows/cellranger.smk"

# Auto‐detect HPC resources
CORES=${SLURM_CPUS_ON_NODE:-$(nproc)}
MEMORY=${SLURM_MEM_PER_NODE:-$(free -g | awk 'NR==2{print $2}')}  # GB

# **Only** the base folder containing all your refdata‐gex builds
TURBO_REF_BASE="/nfs/turbo/umms-parent/Manny_human_ref"

# Output ID for Cell Ranger multi
OUTPUT_ID="10496-MW-reanalysis"