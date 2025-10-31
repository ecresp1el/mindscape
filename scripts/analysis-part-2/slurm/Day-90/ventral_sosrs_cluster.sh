#!/bin/bash
#SBATCH --job-name=ventral_sosrs_cluster_day30
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/cluster_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/cluster_%j.err
#SBATCH --partition=standard
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

source ~/.bashrc
conda activate mindscape-env

# Set Variables
PREPROCESSED_RDS="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon/preprocessed_integrated_jeyoon_test_3.rds"
export MINDSCAPE_PREPROCESSED_RDS="$PREPROCESSED_RDS"
export MINDSCAPE_OUTPUT_DIR="/nfs/turbo/umms-parent/Manny_test/Day-90/ventral_sosrs_output_clustered_test_3"

SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/mindscape_cluster_jeyoon.R"

if [ ! -f "$PREPROCESSED_RDS" ]; then
  echo "❌ Preprocessed RDS not found at: $PREPROCESSED_RDS"
  exit 2
fi

echo "🔁 Running clustering script"
echo "📄 Script: $SCRIPT_PATH"
echo "📁 Output dir: $MINDSCAPE_OUTPUT_DIR"

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Rscript exited with status $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ Clustering job finished successfully"
