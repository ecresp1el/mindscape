#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape Part 2 - Integration and Full Analysis
# ==============================================================================

#SBATCH --job-name=ventral_sosrs_integration_cell_cycle_reg_both
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/integration_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/integration_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment and set paths
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

BASE_DIR="/nfs/turbo/umms-parent/Manny_test"
OUTPUT_DIR="$BASE_DIR/ventral_sosrs_output"
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/mindscape/workflows/ventral_sosrs_workflow_2025_merged/mindscape_integrate_samples_cell_cycle_reg_both.R"

echo "🔁 Running integration script"
echo "📄 Script: $SCRIPT_PATH"
echo "📁 Output dir: $OUTPUT_DIR"

export MINDSCAPE_INPUT_DIR="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output__only_normalized"
export MINDSCAPE_OUTPUT_DIR="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_both"

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
