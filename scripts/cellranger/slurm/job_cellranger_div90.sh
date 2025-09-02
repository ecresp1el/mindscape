#!/usr/bin/env bash
set -euo pipefail

# Minimal, portable Slurm job body without #SBATCH directives.
# Submit with options via sbatch command line (see submit_cellranger_div90.sh).

echo "⏳ [$(date)] Starting DIV90 Cell Ranger job (PID $$)"

# Resolve wrapper path relative to this job script (self-contained)
JOB_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WRAPPER="$JOB_DIR/../scripts/create_project_cellranger_div90.sh"
if [[ ! -f "$WRAPPER" ]]; then
  echo "❌ Wrapper not found next to job script: $WRAPPER" >&2
  exit 1
fi

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
