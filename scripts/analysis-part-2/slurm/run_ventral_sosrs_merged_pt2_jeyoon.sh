#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape Part 2 - Integration and Full Analysis
# ==============================================================================

#SBATCH --job-name=ventral_sosrs_integration_cell_cycle_reg_jeyoon
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/integration_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/integration_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment and set paths
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Verify Seurat v5 is available in this environment
# ------------------------------------------------------------------------------
echo "🔍 Verifying Seurat v5 is available"
Rscript - <<'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suc <- requireNamespace("Seurat", quietly=TRUE) && packageVersion("Seurat") >= "5.0.0"
if (!suc) {
  stop("Seurat v5 not found. Please install Seurat v5 in the 'mindscape-env' conda environment before running this job.")
}
# Also check a dependency like 'spam'
if (!requireNamespace("spam", quietly=TRUE)) {
  stop("Dependency 'spam' not found. Please install it in the conda env (e.g., conda install -c conda-forge r-spam).")
}
# Print versions
library(Seurat)
message("Seurat version: ", packageVersion("Seurat"))
EOF

if [ $? -ne 0 ]; then
  echo "❌ Seurat v5 or dependencies missing in mindscape-env. Aborting."
  exit 1
fi

BASE_DIR="/nfs/turbo/umms-parent/Manny_test"
OUTPUT_DIR="$BASE_DIR/ventral_sosrs_output"
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/mindscape_integrate_samples_cell_cycle_reg_jeyoon.R"

echo "🔁 Running integration script"
echo "📄 Script: $SCRIPT_PATH"
echo "📁 Output dir: $OUTPUT_DIR"

export MINDSCAPE_INPUT_DIR="$BASE_DIR/ventral_sosrs_output__only_normalized"
export MINDSCAPE_OUTPUT_DIR="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon"

# ------------------------------------------------------------------------------
# Run the R integration script
# ------------------------------------------------------------------------------
Rscript "$SCRIPT_PATH"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Rscript exited with status $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ Integration and analysis complete"