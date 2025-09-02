#!/usr/bin/env Rscript

# ==============================================================================
# MarkerGenes.R
#
# PURPOSE:
#   Run DEG analysis for cluster 7 vs all others using Seurat object.
#   Inputs and outputs are provided via environment variables (set in SLURM).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

start_time <- Sys.time()
cat("📊 DEG analysis for cluster 7 started at:", format(start_time), "\n")

# ------------------------------------------------------------------------------
# Get paths from environment variables
# ------------------------------------------------------------------------------
input_path  <- Sys.getenv("DEG_INPUT")
output_path <- Sys.getenv("DEG_OUTPUT")

if (input_path == "" || output_path == "") {
  stop("❌ Environment variables DEG_INPUT or DEG_OUTPUT not set.")
}

cat("📥 Loading Seurat object from:", input_path, "\n")
if (!file.exists(input_path)) {
  stop(paste("❌ Input file does not exist:", input_path))
}
merged <- LoadH5Seurat(input_path)
DefaultAssay(merged) <- "RNA"

print(Assays(merged))
print(slotNames(merged[["RNA"]]))
# ------------------------------------------------------------------------------
# Ensure clustering info is present
# ------------------------------------------------------------------------------
if (!"seurat_clusters" %in% colnames(merged@meta.data)) {
  stop("❌ 'seurat_clusters' metadata not found.")
}
Idents(merged) <- merged$seurat_clusters
cat("✅ Identity set to 'seurat_clusters'\n")

# ------------------------------------------------------------------------------
# Find DEGs for cluster 7
# ------------------------------------------------------------------------------
cat("🧪 Finding markers for cluster 7 vs. all other clusters...\n")
cluster7.markers <- FindMarkers(
  object          = merged,
  ident.1         = 7,
  ident.2         = NULL,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------
if (nrow(cluster7.markers) == 0) {
  warning("⚠️ No DEGs found for cluster 7.")
} else {
  cluster7.markers$gene <- rownames(cluster7.markers)
  cat("💾 Saving DEGs to:", output_path, "\n")
  write.csv(cluster7.markers, file = output_path, row.names = FALSE)
  cat("✅ Saved", nrow(cluster7.markers), "markers\n")
}

# ------------------------------------------------------------------------------
# Done
# ------------------------------------------------------------------------------
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
cat("⏱️ Elapsed time:", round(elapsed, 2), "minutes\n")