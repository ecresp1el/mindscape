#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - Day 30 DEGs (All Clusters)
# ==============================================================================

#SBATCH --job-name=degs_all_clusters_day30
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_all_clusters_day30_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_all_clusters_day30_%j.err
#SBATCH --partition=standard
#SBATCH --time=24:00:00        # May increase if runtime is longer
#SBATCH --mem=180G
#SBATCH --cpus-per-task=4      # Increase if you want faster FindMarkers
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Set input/output paths (Day 30 version)
# ------------------------------------------------------------------------------
BASE_DIR="/nfs/turbo/umms-parent/Manny_test"

# Input: integrated Seurat object from Day 30 analysis
export DEG_INPUT="$BASE_DIR/ventral_sosrs_output_clustered_day_30/clustered_jeyoon_day_30_2.rds"

# Output: directory for DEG results
export DEG_OUTPUT="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_day_30/DEGs_all_clusters_res_0.2"
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

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ DEG script failed with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ DEG analysis complete (Day 30)"
