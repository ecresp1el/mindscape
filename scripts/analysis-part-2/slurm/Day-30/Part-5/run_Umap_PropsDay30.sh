#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - Day 30 UMAP + Cluster Proportions
# ==============================================================================

#SBATCH --job-name=umap_props_day30
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/umap_props_day30_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/umap_props_day30_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Verify Seurat v5 and patchwork are available
# ------------------------------------------------------------------------------
echo "🔍 Verifying Seurat v5 and patchwork availability"
Rscript - <<'EOF'
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
ok <- requireNamespace("Seurat", quietly=TRUE) && packageVersion("Seurat") >= "5.0.0"
if (!ok) stop("Seurat v5 not found in mindscape-env.")
if (!requireNamespace("patchwork", quietly=TRUE)) stop("patchwork not found in mindscape-env.")
EOF

if [ $? -ne 0 ]; then
  echo "❌ Required R packages missing in mindscape-env. Aborting."
  exit 1
fi

# ------------------------------------------------------------------------------
# Set input/output paths (Day 30)
# ------------------------------------------------------------------------------
BASE_DIR="/nfs/turbo/umms-parent/Manny_test"

export DEG_INPUT="$BASE_DIR/ventral_sosrs_output_clustered_day_30/clustered_jeyoon_day_30_2.rds"
export MINDSCAPE_OUTPUT_DIR="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_day_30/umap_props_output_res_0.2"

mkdir -p "$MINDSCAPE_OUTPUT_DIR"

# ------------------------------------------------------------------------------
# Run the plotting script
# ------------------------------------------------------------------------------
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/umap_proportions_figuresDIV30.R"

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

echo "✅ UMAP and cluster proportion plots complete (Day 30)"
