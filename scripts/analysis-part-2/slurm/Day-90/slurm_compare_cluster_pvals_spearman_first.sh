#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - Cluster Matching (Spearman-first Hungarian)
# ==============================================================================
#SBATCH --job-name=compare_clust_spear_hung_prop
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/compare_clust_spear_hung_prop_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/compare_clust_spear_hung_prop_%j.err
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Paths and script
# ------------------------------------------------------------------------------
BASE_DIR="/nfs/turbo/umms-parent/Manny_test"
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/Day-90/CompareClusterPvals_Hungarian_SpearmanFirst.R"

export REFERENCE_DEGS="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs/DEGs_all_clusters/DEGs_all_clusters.csv"
export NEW_DEGS="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs_test_1/DEGs_all_clusters/DEGs_all_clusters.csv"

export OUT_DIR="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs_test_1/match_comparison_spearman_first_prop"
mkdir -p "$OUT_DIR"

# Optional tuning
export MIN_SHARED=50
export SPEARMAN_THRESH=0.5
export MEAN_DIFF_THRESH=0.1
export MAD_THRESH=0.2
export MIN_SHARED_PROP=0.10   # new: minimum proportion of shared genes (0.0–1.0)

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
echo "📊 Running Spearman-first Hungarian cluster matching (with prop_shared filter)"
echo "   REF_DEGS        = $REFERENCE_DEGS"
echo "   NEW_DEGS        = $NEW_DEGS"
echo "   OUT_DIR         = $OUT_DIR"
echo "   MIN_SHARED      = $MIN_SHARED"
echo "   SPEARMAN_THRESH = $SPEARMAN_THRESH"
echo "   MEAN_DIFF_THRESH= $MEAN_DIFF_THRESH"
echo "   MAD_THRESH      = $MAD_THRESH"
echo "   MIN_SHARED_PROP = $MIN_SHARED_PROP"
echo "------------------------------------------------------------------"

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Script failed with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ Completed Spearman-first Hungarian mapping (prop filter active). Results in: $OUT_DIR"