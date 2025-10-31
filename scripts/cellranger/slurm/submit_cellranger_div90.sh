#!/usr/bin/env bash
set -euo pipefail
# Ensures accurate execution of shell script

# Portable Slurm submit wrapper that avoids hardcoded #SBATCH lines.
# Configure via env vars or flags; passes options on sbatch CLI.

usage() {
  cat <<EOF
Usage: ACCOUNT=... TIME=... CPUS=... MEM=... ./scripts/cellranger/slurm/submit_cellranger_div90.sh

Env vars (optional):
  ACCOUNT     Slurm account (e.g., parent0)
  PARTITION   Slurm partition (e.g., standard)
  QOS         Slurm QoS
  TIME        Walltime (e.g., 24:00:00)
  CPUS        CPUs per task (e.g., 32)
  MEM         Memory (e.g., 128G)
  MAIL_USER   Email address for notifications
  JOB_NAME    Job name (default: cellranger_multi_div90)
  LOG_DIR     Logs directory (default: $TEST_DIR/logs or ./logs)

Required runtime config for the job body comes from:
  scripts/cellranger/scripts/config_div90.sh (edit or export env vars)

Examples:
  TEST_DIR=/scratch/$USER/mindscape \
  TURBO_CONFIG_SOURCE=/path/to/DIV90/multi-config.csv \
  PROBE_PATH=/path/to/Probe_set.csv \
  REF_GENOME=/path/to/refdata-gex-GRCh38-2020-A \
  TIME=24:00:00 CPUS=32 MEM=128G ACCOUNT=parent0 \
  ./scripts/cellranger/slurm/submit_cellranger_div90.sh
EOF
}

# Exit upon help request
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage; exit 0
fi

# Resolve job script path relative to this submit script and ensure its presence
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
JOB_SCRIPT="$SUBMIT_DIR/job_cellranger_div90.sh"
[[ -f "$JOB_SCRIPT" ]] || { echo "❌ Missing job script: $JOB_SCRIPT" >&2; exit 1; }

JOB_NAME=${JOB_NAME:-cellranger_multi_div90}

# Prefer logs under TEST_DIR if defined, otherwise local ./logs
# This keeps job logs alongside the run directory for easier tracking.
if [[ -n "${TEST_DIR:-}" ]]; then
  LOG_DIR=${LOG_DIR:-"$TEST_DIR/logs"}
else
  LOG_DIR=${LOG_DIR:-"./logs"}
fi
mkdir -p "$LOG_DIR"

# sbatch argument list (no hardcoded #SBATCH in job body)
ARGS=(
  "--job-name=$JOB_NAME"
  "--output=$LOG_DIR/%x_%j.out"
  "--error=$LOG_DIR/%x_%j.err"
  "--mail-type=END,FAIL"
)

[[ -n "${ACCOUNT:-}"   ]] && ARGS+=("--account=$ACCOUNT")
[[ -n "${PARTITION:-}" ]] && ARGS+=("--partition=$PARTITION")
[[ -n "${QOS:-}"       ]] && ARGS+=("--qos=$QOS")
[[ -n "${TIME:-}"      ]] && ARGS+=("--time=$TIME")
[[ -n "${CPUS:-}"      ]] && ARGS+=("--cpus-per-task=$CPUS")
[[ -n "${MEM:-}"       ]] && ARGS+=("--mem=$MEM")
[[ -n "${MAIL_USER:-}" ]] && ARGS+=("--mail-user=$MAIL_USER")

# Ensure environment variables propagate to the job (portable cluster default)
# This allows config variables (e.g., TEST_DIR, PROBE_PATH) to be visible to the job body.
ARGS+=("--export=ALL")

echo "Submitting with: sbatch ${ARGS[*]} $JOB_SCRIPT"
sbatch "${ARGS[@]}" "$JOB_SCRIPT"
