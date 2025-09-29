#!/usr/bin/env Rscript

# mindscape_integrate_and_analyze.R
# -------------------------------------------------------------------------------
# PURPOSE:
# Load normalized Seurat objects (.rds), integrate them, rebuild RNA layers,
# and run standard Seurat dimensionality reduction and clustering.
# -------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting data integration and analysis\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# Retrieve environment variables set by SLURM
# ------------------------------------------------------------------------------
input_base <- Sys.getenv("MINDSCAPE_INPUT_DIR")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")

if (input_base == "" || output_dir == "") {
  stop("❌ MINDSCAPE_INPUT_DIR or MINDSCAPE_OUTPUT_DIR is not set.")
}

cat(paste0("📂 Input directory: ", input_base, "\n"))
cat(paste0("📁 Output directory: ", output_dir, "\n"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load all .rds normalized objects
# ------------------------------------------------------------------------------
rds_files <- list.files(input_base, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(rds_files) == 0) stop("❌ No .rds files found in input directory.")

seurat_list <- list()
for (f in rds_files) {
  cat("📥 Loading file:", f, "\n")
  obj <- readRDS(f)
  if (!inherits(obj, "Seurat")) {
    stop(paste0("❌ File is not a Seurat object: ", f))
  }
  seurat_list[[basename(f)]] <- obj
}

# ------------------------------------------------------------------------------
# Downstream Seurat analysis
# ------------------------------------------------------------------------------
cat("🛠 Merging Seurat Objects...\n")
merged <- merge(x = seurat_list[[1]], y = seurat_list[-1])

cat("🧬 Normalize again...\n")
merged <- NormalizeData(merged)

cat("🔄 Scaling data...\n")
merged <- ScaleData(merged, vars.to.regress = "CC.Difference", features = rownames(merged))

cat("✨ Finding variable features...\n")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

cat("🎯 Running PCA...\n")
merged <- RunPCA(merged, features = VariableFeatures(merged))

cat("🤝 Integrating Data...\n")
merged <- IntegrateLayers(object = merged, method = CCAIntegration,
                          orig.reduction = "pca",
                          new.reduction = "integrated.cca",
                          verbose = FALSE)

cat("🤝 Finding neighbors and clustering...\n")
merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:10)
merged <- FindClusters(merged, resolution = 0.5)

cat("🗺 Running UMAP...\n")
merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:10)

cat(paste0("🔢 Number of cells in integrated object: ", ncol(merged), "\n"))
cat("✅ Integration and analysis complete\n")

# ------------------------------------------------------------------------------
# Helper: Verify save/load integrity (strict identical check)
# ------------------------------------------------------------------------------
verify_rds_roundtrip <- function(obj, save_path) {
  reload_start <- Sys.time()
  reloaded <- readRDS(save_path)
  reload_end <- Sys.time()
  
  if (identical(obj, reloaded)) {
    cat("✅ Verification passed: reloaded object is identical to saved object\n")
  } else {
    stop("❌ Verification failed: reloaded object differs from saved object")
  }
  
  elapsed_reload <- difftime(reload_end, reload_start, units = "secs")
  cat(paste0("⏱️ Verification completed in ", round(elapsed_reload, 2), " seconds\n"))
  invisible(reloaded)
}

# ------------------------------------------------------------------------------
# Save integrated object (native Seurat v5 .rds format)
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, "integrated_analysis_tryRDS.rds")
cat(paste0("💾 Saving integrated object to: ", save_path, "\n"))
saveRDS(merged, file = save_path)

# Run strict verification
verify_rds_roundtrip(merged, save_path)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write.csv(as.data.frame(Idents(merged)),
          file = file.path(output_dir, "integrated_cluster_ids_tryRDS.csv"))
cat("✅ Cluster IDs saved\n")

png(file.path(output_dir, "integrated_umap_tryRDS.png"),
    width = 800, height = 600)
DimPlot(merged, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

script_end_time <- Sys.time()
elapsed <- difftime(script_end_time, script_start_time, units = "mins")
cat(paste0("⏱️ Completed in ", round(elapsed, 2), " minutes\n"))