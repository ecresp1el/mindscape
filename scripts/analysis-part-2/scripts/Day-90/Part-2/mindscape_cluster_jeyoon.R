#!/usr/bin/env Rscript
# mindscape_cluster_jeyoon.R
# -------------------------------------------------------------------------------
# PURPOSE:
# Load preprocessed + integrated RDS, run clustering + UMAP,
# save results + plots. Includes jhelper verify functions.
# This is step 2 of 5 steps in a pipeline to generate a UMAP and cell-type proportions figure for a day 90 timepoint. This step focuses on clustering preprocessed data. 
# -------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# INPUTS / OUTPUTS
# ------------------------------------------------------------------------------
# INPUTS:
#   - Environment variable MINDSCAPE_PREPROCESSED_RDS:
#       Path to the preprocessed + integrated Seurat object (.rds)
#       generated from the previous pipeline step
#   - Environment variable MINDSCAPE_OUTPUT_DIR:
#       Directory where clustering outputs will be written
#
# OUTPUTS:
#   - Clustered + UMAP-embedded Seurat object:
#       clustered_jeyoon_test_3.rds
#   - Cluster identity table:
#       cluster_ids_jeyoon_test_3.csv
#   - UMAP visualization:
#       umap_jeyoon_test_3.png
#   - Verification log messages confirming object integrity
# ------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting clustering + UMAP analysis\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# jhelper functions
# ---------------------------
jhelper_check_seurat <- function(min_version = "5.0.0") {
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (packageVersion("Seurat") < package_version(min_version)) stop("Seurat version too old")
  message("Seurat version: ", packageVersion("Seurat"))
  invisible(TRUE)
}

jhelper_assert_file_exists <- function(path) {
  if (!file.exists(path)) stop(paste0("Required file not found: ", path))
  invisible(TRUE)
}

jhelper_verify_rds_roundtrip <- function(obj, save_path) {
  if (!file.exists(save_path)) stop("Saved file not found for verification.")
  reloaded <- readRDS(save_path)
  if (identical(obj, reloaded)) {
    cat("✅ Verification passed: reloaded identical to saved\n")
  } else {
    stop("❌ Verification failed: objects differ after reload")
  }
  invisible(TRUE)
}

# ---------------------------
# Env & inputs
# ---------------------------
input_rds <- Sys.getenv("MINDSCAPE_PREPROCESSED_RDS")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")
if (input_rds == "" || output_dir == "") {
  stop("❌ MINDSCAPE_PREPROCESSED_RDS or MINDSCAPE_OUTPUT_DIR not set.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

jhelper_check_seurat()
jhelper_assert_file_exists(input_rds)

cat("📥 Loading checkpoint RDS:", input_rds, "\n")
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("Loaded object is not a Seurat object.")

# ---------------------------
# Clustering + UMAP
# ---------------------------
cat("🔍 Finding neighbors + clusters...\n")
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)

cat("🗺 Running UMAP...\n")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:10)

# ---------------------------
# Save outputs + verify
# ---------------------------
save_path <- file.path(output_dir, "clustered_jeyoon_test_3.rds")
cat("💾 Saving clustered object:", save_path, "\n")
saveRDS(obj, save_path)

jhelper_verify_rds_roundtrip(obj, save_path)

cat("💾 Writing cluster IDs CSV\n")
write.csv(as.data.frame(Idents(obj)), file = file.path(output_dir, "cluster_ids_jeyoon_test_3.csv"))

png(file.path(output_dir, "umap_jeyoon_test_3.png"), width = 800, height = 600)
DimPlot(obj, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

cat("⏱ Completed in ", round(difftime(Sys.time(), script_start_time, units = "mins"), 2), " mins\n")
