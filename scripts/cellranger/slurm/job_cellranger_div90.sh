#!/usr/bin/env bash
set -euo pipefail

# Minimal, portable Slurm job body without #SBATCH directives.
# Submit with options via sbatch command line (see submit_cellranger_div90.sh).

echo "⏳ [$(date)] Starting DIV90 Cell Ranger job (PID $$)"

# Resolve wrapper path robustly in Slurm (script is spooled to /var/spool/...)
resolve_wrapper() {
  local candidates=()
  # 1) Explicit override
  if [[ -n "${WRAPPER_PATH:-}" ]]; then
    candidates+=("$WRAPPER_PATH")
  fi
  # 2) Based on submit dir (where sbatch was invoked)
  if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    candidates+=("$SLURM_SUBMIT_DIR/../scripts/create_project_cellranger_div90.sh")
  fi
  # 3) Based on current working directory (often equals SLURM_SUBMIT_DIR)
  candidates+=("$PWD/../scripts/create_project_cellranger_div90.sh")
  # 4) Fallback: attempt relative to the job script location (may fail in spool)
  local job_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  candidates+=("$job_dir/../scripts/create_project_cellranger_div90.sh")

  for c in "${candidates[@]}"; do
    if [[ -f "$c" ]]; then
      echo "$c"; return 0
    fi
  done
  return 1
}

WRAPPER="$(resolve_wrapper)" || { echo "❌ Wrapper not found via SLURM_SUBMIT_DIR/PWD/BASH_SOURCE fallbacks" >&2; exit 1; }

# Optional: requeue on pre-timeout signal if available
_requeue_on_timeout() {
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "⚠️ Caught pre-timeout signal. Requeuing job ID $SLURM_JOB_ID..."
    scontrol requeue "$SLURM_JOB_ID" || true
  fi
  exit 99
}
trap _requeue_on_timeout SIGUSR1 || true

# Execute wrapper (no dependency on submit dir)
bash "$WRAPPER"

echo "✅ [$(date)] DIV90 Cell Ranger job finished"
