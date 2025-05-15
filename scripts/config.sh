#!/usr/bin/env bash
#
# mindscape/scripts/config.sh
# — All paths & settings in one place —
#

# Where to run the test
TEST_DIR="/nfs/turbo/umms-parent/Manny_test"

# The CellRanger “multi” config source you copy & patch
TURBO_CONFIG_SOURCE="/nfs/turbo/umms-parent/Carmen_Miranda_scRNAseq /90 Day results/10x_analysis_10496-MW/Sample_10496-MW-Pool01/config.csv"

# Flex probe set (on Turbo storage)
PROBE_PATH="/nfs/turbo/umms-parent/Manny_probe_set/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"

# Where your Snakemake workflow lives (relative to repo root)
SNAKEFILE="mindscape/workflows/cellranger.smk"

# Auto‐detect HPC resources
CORES=${SLURM_CPUS_ON_NODE:-$(nproc)}
MEMORY=${SLURM_MEM_PER_NODE:-$(free -g | awk 'NR==2{print $2}')}

# Reference genome (on Turbo storage)
TURBO_REF_BASE="/nfs/turbo/umms-parent/Manny_human_ref"
REF_SUBPATH="GRCh38-GENCODEv35_build/refdata-gex-GRCh38-2020-A"

# Output ID for Cell Ranger multi
OUTPUT_ID="10496-MW-reanalysis"