#!/usr/bin/env bash
#
#SBATCH --job-name=cellranger_multi_jeyoon_day_30
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/cellranger_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/cellranger_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --signal=B:USR1@300             # 🔁 Requeue 5 minutes before timeout
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=oltij@umich.edu

# ─────────────────────────────────────────────────────────────────────────────
# cellranger_submission_jeyoon_day_30_slurm.sh
# Slurm wrapper for scripts/cellranger/scripts/create_project_cellranger_jeyoon_day_30.sh
# Adds safe requeueing support if timeout is near
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# 🛑 Requeue trap
_requeue_on_timeout() {
  echo "⚠️ Caught pre-timeout signal. Requeuing job ID $SLURM_JOB_ID..."
  scontrol requeue "$SLURM_JOB_ID"
  exit 99
}
trap _requeue_on_timeout SIGUSR1

echo "⏳ [$(date)] Starting Slurm job $SLURM_JOB_NAME (ID $SLURM_JOB_ID)"

# 🧭 Ensure working directory is the submission dir
cd "$SLURM_SUBMIT_DIR"

# 🧼 Load environment modules (if needed)
module purge

# 📍 Relative path to your wrapper script
WRAPPER="scripts/cellranger/scripts/create_project_cellranger_jeyoon_day_30.sh"

# 🧪 Sanity check
echo "🔍 SLURM_SUBMIT_DIR = $SLURM_SUBMIT_DIR"
echo "🔍 Looking for wrapper at: $WRAPPER"
if [[ ! -f "$WRAPPER" ]]; then
  echo "❌ ERROR: Wrapper script not found at $WRAPPER"
  exit 1
fi

# ✅ Make sure it's executable
chmod +x "$WRAPPER"

# 🚀 Launch pipeline
echo "🚀 Launching cellranger multi workflow..."
bash "$WRAPPER"

echo "✅ [$(date)] Slurm job $SLURM_JOB_NAME (ID $SLURM_JOB_ID) completed"