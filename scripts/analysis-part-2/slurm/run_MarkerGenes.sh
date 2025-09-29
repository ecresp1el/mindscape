#!/bin/bash

#SBATCH --job-name=degs_all_clusters
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_all_clusters_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_all_clusters_%j.err
#SBATCH --partition=standard
#SBATCH --time=04:00:00   # All clusters may take longer
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4  # Can increase to speed up FindMarkers
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ──────────────────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate mindscape-env

# Paths passed via environment variables
export DEG_INPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon/integrated_analysis_tryRDS.rds"

# OUTPUT should be a directory for multiple cluster CSVs
export DEG_OUTPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_tryRDS/DEGs_all_clusters"
mkdir -p "$DEG_OUTPUT"

SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/MarkerGenes_AllClusters.R"

echo "📊 Running DEG script with:"
echo "   INPUT  = $DEG_INPUT"
echo "   OUTPUT DIR = $DEG_OUTPUT"

Rscript "$SCRIPT_PATH"