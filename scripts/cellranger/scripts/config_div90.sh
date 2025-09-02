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
# Notes:
#   - Avoid hardcoding absolute paths here when possible. Prefer exporting
#     paths per environment or mounting shared locations.
#   - All variables here can be overridden via environment variables.

# Where to run the test (required)
# Default safely under /nfs/turbo/umms-parent to avoid clobbering others.
# Creates a user-scoped, timestamped workdir unless overridden.
RUNSTAMP=${RUNSTAMP:-$(date +%Y%m%d_%H%M%S)}
TEST_DIR=${TEST_DIR:-/nfs/turbo/umms-parent/${USER:-unknown}/mindscape_div90/${RUNSTAMP}}

# Source multi-config CSV for DIV90 (required)
# Example: /nfs/.../90 Day results/fastqs_10496-MW/multi-config.csv
TURBO_CONFIG_SOURCE=${TURBO_CONFIG_SOURCE:-}

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
# to scripts/cellranger/ (e.g., "cellranger.smk").
SNAKEFILE=${SNAKEFILE:-cellranger.smk}

# Output ID for Cell Ranger multi target (Snakemake rule name)
OUTPUT_ID=${OUTPUT_ID:-div90-reanalysis}

# Resources (auto-detected; can be overridden)
CORES=${CORES:-${SLURM_CPUS_ON_NODE:-$(command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu 2>/dev/null || echo 1)}}
MEMORY_GB=${MEMORY_GB:-${SLURM_MEM_PER_NODE:-}}
