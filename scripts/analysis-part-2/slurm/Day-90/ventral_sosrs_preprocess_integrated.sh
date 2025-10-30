#!/bin/bash
#SBATCH --job-name=ventral_sosrs_preprocess_integrated_day90
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/preprocess_integrated_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/preprocess_integrated_%j.err
#SBATCH --partition=largemem
#SBATCH --time=4:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

source ~/.bashrc
conda activate mindscape-env

# Verify Seurat v5 and basic deps (quick check in R)
echo "🔍 Verifying Seurat v5 is available"
Rscript - <<'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("Seurat", quietly=TRUE) || packageVersion("Seurat") < "5.0.0") {
  stop("Seurat v5 not found in environment 'mindscape-env'. Please install it.")
}
message("Seurat OK: ", packageVersion("Seurat"))
EOF
if [ $? -ne 0 ]; then
  echo "❌ Seurat v5 or dependencies missing. Aborting."
  exit 1
fi

# Set variables
BASE_DIR="/nfs/turbo/umms-parent/Manny_test"
export MINDSCAPE_INPUT_DIR="$BASE_DIR/ventral_sosrs_output__only_normalized"
export MINDSCAPE_OUTPUT_DIR="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon"

SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/Day-90/mindscape_preprocess_integrated_jeyoon.R"

echo "🔁 Running preprocess/integrate script"
echo "📄 Script: $SCRIPT_PATH"
echo "📁 Output dir: $MINDSCAPE_OUTPUT_DIR"

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Rscript exited with status $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ Preprocessing + integration job finished successfully"