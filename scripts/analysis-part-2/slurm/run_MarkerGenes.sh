#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - DEGs (All Clusters)
# ==============================================================================

#SBATCH --job-name=degs_all_clusters_tryRDS
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_all_clusters_tryRDS_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_all_clusters_tryRDS_%j.err
#SBATCH --partition=standard
#SBATCH --time=04:00:00          # Increase if FindMarkers runs longer
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4        # Increase if you want faster DE analysis
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Set input/output paths
# ------------------------------------------------------------------------------
BASE_DIR="/nfs/turbo/umms-parent/Manny_test"

# Input: integrated Seurat object
export DEG_INPUT="$BASE_DIR/ventral_sosrs_output_clustered_test_1/clustered_jeyoon_test_1.rds"

# Output: directory for DEG results
export DEG_OUTPUT="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs_test_1/DEGs_all_clusters"
mkdir -p "$DEG_OUTPUT"

# Path to DEG R script
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/MarkerGenes_AllClusters.R"

# ------------------------------------------------------------------------------
# Run the DEG analysis script
# ------------------------------------------------------------------------------
echo "📊 Running DEG script with:"
echo "   INPUT      = $DEG_INPUT"
echo "   OUTPUT DIR = $DEG_OUTPUT"
echo "   SCRIPT     = $SCRIPT_PATH"
echo "------------------------------------------------------------------"

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ DEG script failed with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ DEG analysis complete"