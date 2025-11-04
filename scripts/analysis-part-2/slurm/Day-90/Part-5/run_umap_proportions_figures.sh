#!/bin/bash

# ==============================================================================
# INPUTS / OUTPUTS
# ==============================================================================
# INPUTS:
#   - DEG_INPUT:
#       Path to clustered Seurat RDS file for Day 90 (with UMAP and cluster metadata).
#   - MINDSCAPE_OUTPUT_DIR:
#       Directory where the generated UMAP and proportion plots will be saved.
#   - Active environment:
#       mindscape-env (must include Seurat ≥5.0 and patchwork)
#
# OUTPUTS:
#   - Verified Seurat object with cluster names appended:
#       clustered_day90_with_cluster_names_2.rds
#   - Combined UMAP + cluster proportion plots:
#       day90_res_0.5_3_combined.{png,svg,pdf}
#   - Logs:
#       stdout → umap_props_<jobid>.out
#       stderr → umap_props_<jobid>.err
# ==============================================================================


#SBATCH --job-name=umap_props
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/umap_props_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/umap_props_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Verify Seurat v5 is available
# ------------------------------------------------------------------------------
echo "🔍 Verifying Seurat v5 is available"
Rscript - <<'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suc <- requireNamespace("Seurat", quietly=TRUE) && packageVersion("Seurat") >= "5.0.0"
if (!suc) stop("Seurat v5 not found in mindscape-env.")
if (!requireNamespace("patchwork", quietly=TRUE)) stop("patchwork not found in mindscape-env.")
EOF

if [ $? -ne 0 ]; then
  echo "❌ Required R packages missing in mindscape-env. Aborting."
  exit 1
fi

# ------------------------------------------------------------------------------
# Set input/output paths
# ------------------------------------------------------------------------------
export DEG_INPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_clustered/clustered_jeyoon.rds"
export MINDSCAPE_OUTPUT_DIR="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output"

mkdir -p "$MINDSCAPE_OUTPUT_DIR"

# ------------------------------------------------------------------------------
# Run the plotting script
# ------------------------------------------------------------------------------
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/Day-90/Part-5/umap_proportions_figures.R"

echo "🔁 Running plotting script"
echo "📄 Script: $SCRIPT_PATH"
echo "📂 Input RDS: $DEG_INPUT"
echo "📁 Output dir: $MINDSCAPE_OUTPUT_DIR"

Rscript "$SCRIPT_PATH"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Rscript exited with status $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ UMAP and cluster proportion plots complete"