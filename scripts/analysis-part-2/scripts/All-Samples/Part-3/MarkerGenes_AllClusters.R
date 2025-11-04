#!/usr/bin/env Rscript

# ==============================================================================
# MarkerGenes_AllClusters.R
#
# PURPOSE:
#   Load normalized Seurat objects (.rds) and run DEG analysis
#   for all clusters vs all others. Inputs and outputs are
#   provided via environment variables (set in SLURM).
# This is step 3 of 5 steps in a pipeline to generate a UMAP and cell-type proportions figure for day 30 and day 90 timepoints. This step focuses on finding all of the marker genes for all clusters.
# ==============================================================================

# ------------------------------------------------------------------------------
# INPUTS / OUTPUTS
# ------------------------------------------------------------------------------
# INPUTS:
#   - Environment variable DEG_INPUT:
#       Path to integrated + clustered Seurat object (.rds)
#       (output from the previous pipeline step)
#   - Environment variable DEG_OUTPUT:
#       Directory where DEG result files will be saved
#
# OUTPUTS:
#   - Per-cluster DEG tables:
#       DEGs_cluster_<cluster_id>.csv
#   - Combined DEG table (all clusters):
#       DEGs_all_clusters.csv
#   - Verification messages confirming file integrity
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Helper: robust CSV verification
# ------------------------------------------------------------------------------
verify_csv_integrity <- function(df, csv_path) {
  reloaded <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  # Ensure same column order
  reloaded <- reloaded[, colnames(df), drop = FALSE]
  
  # Convert numeric to character for columns that are character in df
  char_cols <- names(df)[sapply(df, is.character)]
  for (col in char_cols) {
    reloaded[[col]] <- as.character(reloaded[[col]])
  }
  
  # Compare approximately for numerics, exact for characters
  ok <- all.equal(df, reloaded, tolerance = 1e-8, check.attributes = FALSE)
  
  if (!isTRUE(ok)) {
    stop("❌ Verification failed for CSV: ", csv_path, "\nDetails: ", ok)
  }
  
  message("✅ Verification passed: ", basename(csv_path))
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Retrieve paths from environment variables
# ------------------------------------------------------------------------------
input_file  <- Sys.getenv("DEG_INPUT")   # .rds Seurat object
output_path <- Sys.getenv("DEG_OUTPUT")  # output folder

if (input_file == "" || output_path == "") {
  stop("❌ DEG_INPUT or DEG_OUTPUT not set.")
}
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load Seurat object (.rds, v5 format)
# ------------------------------------------------------------------------------
cat("📥 Loading Seurat object from:", input_file, "\n")
obj <- readRDS(input_file)
if (!inherits(obj, "Seurat")) {
  stop("❌ Input file is not a Seurat object: ", input_file)
}
DefaultAssay(obj) <- "RNA"
cat("✅ Seurat object loaded successfully\n")

# Join layers before DEG
cat("🔗 Joining RNA assay layers...\n")
obj <- JoinLayers(obj)
cat("✅ Layers joined successfully\n")

# ------------------------------------------------------------------------------
# DEG analysis
# ------------------------------------------------------------------------------
start_time <- Sys.time()
cat("📊 DEG analysis for all clusters started at:", format(start_time), "\n")

# Check for clusters
if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("❌ 'seurat_clusters' metadata not found.")
}
Idents(obj) <- obj$seurat_clusters
cat("✅ Identity set to 'seurat_clusters'\n")

# Sort clusters
clusters <- sort(unique(Idents(obj)))
cat("🔹 Found clusters:", paste(clusters, collapse = ", "), "\n")

# Find markers
all_markers <- list()
for (clust in clusters) {
  cat("🧪 Finding markers for cluster", clust, "vs all others...\n")
  markers <- tryCatch(
    FindMarkers(
      object          = obj,
      ident.1         = clust,
      ident.2         = NULL,
      only.pos        = TRUE,
      min.pct         = 0.01,
      logfc.threshold = 0.1,
      test.use        = "wilcox"
    ),
    error = function(e) {
      warning(paste("⚠️ Failed to calculate DEGs for cluster", clust, ":", conditionMessage(e)))
      return(NULL)
    }
  )
  
  if (!is.null(markers) && nrow(markers) > 0) {
    markers$gene <- rownames(markers)
    markers$cluster <- clust  # keep original type, verification handles character/numeric differences
    all_markers[[as.character(clust)]] <- markers
    cluster_file <- file.path(output_path, paste0("DEGs_cluster_", clust, ".csv"))
    write.csv(markers, file = cluster_file, row.names = FALSE)
    verify_csv_integrity(markers, cluster_file)
    cat("💾 Saved", nrow(markers), "markers for cluster", clust, "to", cluster_file, "\n")
  } else {
    cat("⚠️ No DEGs found for cluster", clust, "\n")
  }
}

if (length(all_markers) > 0) {
  combined_markers <- bind_rows(all_markers)
  combined_file <- file.path(output_path, "DEGs_all_clusters.csv")
  write.csv(combined_markers, combined_file, row.names = FALSE)
  verify_csv_integrity(combined_markers, combined_file)
  cat("💾 Saved combined DEG table with", nrow(combined_markers), "rows to", combined_file, "\n")
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
cat("⏱️ Elapsed time:", round(elapsed, 2), "minutes\n")

# ------------------------------------------------------------------------------
# Print any captured warnings explicitly
# ------------------------------------------------------------------------------
w <- warnings()
if (length(w) > 0) {
  cat("⚠️ Warnings generated during DEG analysis:\n")
  print(w)
}
