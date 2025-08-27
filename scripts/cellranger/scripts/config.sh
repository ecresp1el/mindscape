#!/usr/bin/env bash
#
# mindscape/scripts/config.sh
#

# Where to run the test
TEST_DIR="/nfs/turbo/umms-parent/Manny_test"

# Original multi‐config CSV
TURBO_CONFIG_SOURCE="/nfs/turbo/umms-parent/Accessible_multi-config_csvs/Carmen_Miranda_scRNAseq /30 Day results/fastqs_9583-MW/multi-config.csv"

# Flex probe set (on Turbo storage)
PROBE_PATH="/nfs/turbo/umms-parent/10X_Human_Refs/2020-A/Probe_set/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"

# Where your Snakemake workflow lives (relative)
SNAKEFILE="mindscape/workflows/cellranger.smk"

# Auto‐detect HPC resources
CORES=${SLURM_CPUS_ON_NODE:-$(nproc)}
MEMORY=${SLURM_MEM_PER_NODE:-$(free -g | awk 'NR==2{print $2}')}  # GB

# The base folder containing all your refdata‐gex builds
TURBO_REF_BASE="/nfs/turbo/umms-parent/10X_Human_Refs"

# Reference subpath
REF_SUBPATH="2020-A/Ref_genome/refdata-gex-GRCh38-2020-A"

# Output ID for Cell Ranger multi
OUTPUT_ID="10496-MW-reanalysis"