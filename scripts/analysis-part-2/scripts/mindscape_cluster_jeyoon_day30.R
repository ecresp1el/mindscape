#!/usr/bin/env Rscript

# mindscape_cluster_jeyoon_day_30.R
# -------------------------------------------------------------------------------
# PURPOSE:
# Load preprocessed + integrated RDS, run clustering + UMAP,
# and save results + plots.
# This is step 2 of 5 steps in a pipeline to generate a UMAP and cell-type proportions figure for a day 30 timepoint. This step focuses on clustering preprocessed data. 
# -------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting clustering + UMAP analysis\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Retrieve environment variables
# ------------------------------------------------------------------------------
input_rds <- Sys.getenv("MINDSCAPE_PREPROCESSED_RDS")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")

if (input_rds == "" || output_dir == "") {
  stop("❌ MINDSCAPE_PREPROCESSED_RDS or MINDSCAPE_OUTPUT_DIR not set.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load checkpoint object
# ------------------------------------------------------------------------------
cat("📥 Loading checkpoint RDS:", input_rds, "\n")
obj <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# Clustering + UMAP
# ------------------------------------------------------------------------------
cat("🔍 Finding neighbors + clusters...\n")
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.1)

cat("🗺 Running UMAP...\n")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:10)

# ------------------------------------------------------------------------------
# Helper function for RDS verification
# ------------------------------------------------------------------------------
verify_rds_roundtrip <- function(obj, save_path) {
  reloaded <- readRDS(save_path)
  if (identical(obj, reloaded)) {
    cat("✅ Verification passed: reloaded identical to saved\n")
  } else {
    stop("❌ Verification failed: objects differ")
  }
}

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, "clustered_jeyoon_day_30_3.rds")
cat("💾 Saving clustered object:", save_path, "\n")
saveRDS(obj, save_path)

verify_rds_roundtrip(obj, save_path)

write.csv(as.data.frame(Idents(obj)),
          file = file.path(output_dir, "cluster_ids_jeyoon_day_30_3.csv"))
cat("✅ Cluster IDs saved\n")

png(file.path(output_dir, "umap_jeyoon_day_30_3.png"), width = 800, height = 600)
DimPlot(obj, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

cat("⏱ Completed in ", round(difftime(Sys.time(), script_start_time, units="mins"), 2), " mins\n")
