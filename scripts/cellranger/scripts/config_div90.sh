#!/usr/bin/env bash
#
# scripts/cellranger/scripts/config_div90.sh
# Parameterized configuration for DIV90 Cell Ranger pipeline.
#
# Usage:
#   - Edit values below OR export env vars before running the wrapper.
#   - Required: TEST_DIR, TURBO_CONFIG_SOURCE, PROBE_PATH, and a reference
#     (either REF_GENOME or TURBO_REF_BASE+REF_SUBPATH).
#
# Quick examples:
#   Local dry-run (plan only):
#     export TURBO_CONFIG_SOURCE="/nfs/turbo/.../config.csv"
#     DRY_RUN=1 bash scripts/cellranger/scripts/create_project_cellranger_div90.sh
#
#   Slurm submission (actual run):
#     TIME=48:00:00 CPUS=32 MEM=128G ACCOUNT=parent0 \
#     scripts/cellranger/div90/drive_div90.sh slurm
#
# Notes:
#   - Avoid hardcoding absolute paths here when possible. Prefer exporting
#     paths per environment or mounting shared locations.
#   - All variables here can be overridden via environment variables.

# Where to run the test (required)
# - Default safely under /nfs/turbo/umms-parent to avoid clobbering others.
# - Creates a user-scoped, timestamped workdir unless overridden (RUNSTAMP controls suffix).
RUNSTAMP=${RUNSTAMP:-$(date +%Y%m%d_%H%M%S)}
TEST_DIR=${TEST_DIR:-/nfs/turbo/umms-parent/${USER:-unknown}/mindscape_div90/${RUNSTAMP}} # Test directory per run ans per user

# Source multi-config CSV for DIV90 (required)
# Default: prior DIV90 config used in cellranger pipeline notes
# Note: path contains spaces; keep quotes when overriding.
TURBO_CONFIG_SOURCE=${TURBO_CONFIG_SOURCE:-"/nfs/turbo/umms-parent/Accessible_multi-config_csvs/Carmen_Miranda_scRNAseq /90 Day results/fastqs_10496-MW/multi-config.csv"} # For running cellranger

# FASTQ path handling (optional)
# - Some historical multi_config.csv files reference FASTQs under restricted or
#   cluster-inaccessible prefixes (e.g., /nfs/turbo/agc-data/...).
# - Set one of the following to normalize the fastqs column in [libraries]:
#   1) FASTQS_DIR: force the fastqs column for all rows to this directory
#   2) FASTQS_REPLACE_FROM + FASTQS_REPLACE_TO: prefix replacement within fastqs
#
# Examples:
#   export FASTQS_DIR="/nfs/turbo/umms-parent/Accessible_fastqs/10496-MW/fastqs_10496-MW"
#   export FASTQS_REPLACE_FROM="/nfs/turbo/agc-data/processing"
#   export FASTQS_REPLACE_TO="/nfs/turbo/umms-parent/Accessible_fastqs"
FASTQS_DIR=${FASTQS_DIR:-}
FASTQS_REPLACE_FROM=${FASTQS_REPLACE_FROM:-}
FASTQS_REPLACE_TO=${FASTQS_REPLACE_TO:-}

# Accessibility enforcement (DIV90 default)
# - ENFORCE_ACCESSIBLE=1 ensures curated accessible CSV and FASTQ prefixes are used.
# - ALLOWED_FASTQS_PREFIX restricts FASTQ directories to a known accessible root.
ENFORCE_ACCESSIBLE=${ENFORCE_ACCESSIBLE:-1}
ALLOWED_FASTQS_PREFIX=${ALLOWED_FASTQS_PREFIX:-/nfs/turbo/umms-parent}

# Flex probe set (required)
# Default: previously used 10X Human Refs (2020-A) path on turbo
# You can override via env if your cluster uses a different location.
PROBE_PATH=${PROBE_PATH:-/nfs/turbo/umms-parent/10X_Human_Refs/2020-A/Probe_set/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv}

# Cell Ranger reference
# Option A: Set REF_GENOME directly
# Default: previously used 10X Human Refs (2020-A) genome on turbo
REF_GENOME=${REF_GENOME:-/nfs/turbo/umms-parent/10X_Human_Refs/2020-A/Ref_genome/refdata-gex-GRCh38-2020-A}
# Option B: Or set these and the wrapper derives REF_GENOME="$TURBO_REF_BASE/$REF_SUBPATH"
TURBO_REF_BASE=${TURBO_REF_BASE:-/nfs/turbo/umms-parent/10X_Human_Refs}
REF_SUBPATH=${REF_SUBPATH:-2020-A/Ref_genome/refdata-gex-GRCh38-2020-A}

# Snakemake snakefile path inside the cellranger folder
# You can override with an absolute path, or a path relative 
# to scripts/cellranger/ (e.g., "cellranger.smk"). Ensure relative paths are readable with a path base.
SNAKEFILE=${SNAKEFILE:-cellranger.smk}

# Output ID for Cell Ranger multi target (Snakemake rule name)
OUTPUT_ID=${OUTPUT_ID:-div90-reanalysis}

# Resources (auto-detected; can be overridden)
# - CORES controls Snakemake parallelism for local runs.
# - For Slurm, CPUS/MEM are set on submission; Snakemake itself still uses CORES inside the job.
CORES=${CORES:-${SLURM_CPUS_ON_NODE:-$(command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu 2>/dev/null || echo 1)}}
MEMORY_GB=${MEMORY_GB:-${SLURM_MEM_PER_NODE:-}}
