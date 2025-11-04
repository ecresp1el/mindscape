#!/usr/bin/env Rscript

# mindscape_integrate_and_analyze.R

# -------------------------------------------------------------------------------
# PURPOSE:
# Load normalized Seurat objects (.h5Seurat), integrate them, and run
# standard Seurat dimensionality reduction and clustering.
# -------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting data integration and analysis\n")

library(Seurat)
library(SeuratDisk)
library(dplyr)

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
# Load all .h5Seurat normalized objects
# ------------------------------------------------------------------------------
h5_files <- list.files(input_base, pattern = "\\.h5Seurat$", full.names = TRUE, recursive = TRUE)
if (length(h5_files) == 0) stop("❌ No .h5Seurat files found in input directory.")

seurat_list <- list()
for (f in h5_files) {
  cat(paste0("📥 Loading file: ", f, "\n"))
  obj <- LoadH5Seurat(f)
  seurat_list[[basename(f)]] <- obj
}

# ------------------------------------------------------------------------------
# Preprocessing for each object
# ------------------------------------------------------------------------------
cat("🧬 Normalize again...\n")
seurat_list <- lapply(seurat_list, NormalizeData)

cat("🔄 Scaling data...\n")
seurat_list <- lapply(seurat_list, function(obj) {
  ScaleData(obj, vars.to.regress = "CC.Difference", features = rownames(obj))
})

cat("✨ Finding variable features...\n")
seurat_list <- lapply(seurat_list, function(obj) {
  FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
})

cat("🎯 Running PCA...\n")
seurat_list <- lapply(seurat_list, function(obj) {
  RunPCA(obj, features = VariableFeatures(obj))
})

# ------------------------------------------------------------------------------
# Integration
# ------------------------------------------------------------------------------
cat("🤝 Integrating Data...\n")
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = 2000, reduction = "cca", dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Define all.genes after integration
all.genes <- rownames(integrated)

cat("🔄 Scaling integrated data...\n")
integrated <- ScaleData(integrated, features = all.genes)

cat("🎯 Running PCA on integrated object...\n")
integrated <- RunPCA(integrated, reduction.name = "pca", assay = "integrated")

cat("🔍 Finding neighbors and clustering...\n")
integrated <- FindNeighbors(integrated, dims = 1:10)
integrated <- FindClusters(integrated, resolution = 0.5)

cat("🗺 Running UMAP...\n")
integrated <- RunUMAP(integrated, dims = 1:10)

cat(paste0("🔢 Number of cells in integrated object: ", ncol(integrated), "\n"))
cat("✅ Integration and analysis complete\n")

# ------------------------------------------------------------------------------
# Save integrated object
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, "integrated_analysis_2.h5Seurat")
cat(paste0("💾 Saving integrated object to: ", save_path, "\n"))
SaveH5Seurat(integrated, filename = save_path, overwrite = TRUE)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write.csv(as.data.frame(Idents(integrated)), file = file.path(output_dir, "integrated_cluster_ids_2.csv"))
cat("✅ Cluster IDs saved\n")

png(file.path(output_dir, "integrated_umap_cell_cycle_reg_jeyoon_2.png"), width = 800, height = 600)
DimPlot(integrated, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

script_end_time <- Sys.time()
elapsed <- difftime(script_end_time, script_start_time, units = "mins")
cat(paste0("⏱️ Completed in ", round(elapsed, 2), " minutes\n"))