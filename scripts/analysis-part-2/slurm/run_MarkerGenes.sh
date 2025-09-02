#!/bin/bash

#SBATCH --job-name=degs_cluster7
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_cluster7_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_cluster7_%j.err
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ──────────────────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate mindscape-env

# Paths passed via environment variables
export DEG_INPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_day_30/integrated_analysis_jeyoon_day_30.h5Seurat"
export DEG_OUTPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_day_30/cluster7_DEGs.csv"

SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/MarkerGenes.R"

echo "📊 Running DEG script with:"
echo "   INPUT  = $DEG_INPUT"
echo "   OUTPUT = $DEG_OUTPUT"

Rscript "$SCRIPT_PATH"