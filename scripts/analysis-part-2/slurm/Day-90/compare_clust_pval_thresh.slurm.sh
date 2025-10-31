#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - Cluster Matching (Direct P-Value Similarity)
# ==============================================================================
#SBATCH --job-name=compare_clust_pval_sim
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/compare_clust_pval_sim_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/compare_clust_pval_sim_%j.err
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
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/Day-90/CompareClusterPvals_DirectThreshold_Summary.R"

export REFERENCE_DEGS="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs/DEGs_all_clusters/DEGs_all_clusters.csv"
export NEW_DEGS="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs_test_1/DEGs_all_clusters/DEGs_all_clusters.csv"

export OUT_DIR="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_DEGs_test_1/match_comparison_pval_similarity"
mkdir -p "$OUT_DIR"

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
# Maximum allowed |Δp| per gene for a match
export PVALUE_DIFF_THRESH=0.05   

# Minimum shared genes required between clusters
export MIN_SHARED=50             

# Require proportion of similar p-values (set to 1.0 = 100% for now)
export PROP_SIMILAR_REQUIRED=1.0 

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
echo "📊 Running direct p-value similarity cluster matching"
echo "   REF_DEGS              = $REFERENCE_DEGS"
echo "   NEW_DEGS              = $NEW_DEGS"
echo "   OUT_DIR               = $OUT_DIR"
echo "   PVALUE_DIFF_THRESH    = $PVALUE_DIFF_THRESH"
echo "   MIN_SHARED            = $MIN_SHARED"
echo "   PROP_SIMILAR_REQUIRED = $PROP_SIMILAR_REQUIRED"
echo "------------------------------------------------------------------"

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?

# ------------------------------------------------------------------------------
# Post-run check
# ------------------------------------------------------------------------------
if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Script failed with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ Completed direct p-value similarity mapping."
echo "📁 Results available in: $OUT_DIR"
