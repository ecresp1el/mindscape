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
  cat("📥 Loading file:", f, "\n")
  obj <- LoadH5Seurat(f)
  # Rebuild the RNA assay as a proper v5 Assay
  counts_mat <- GetAssayData(obj, layer = "counts")
  data_mat   <- GetAssayData(obj, layer = "data")
  new_assay <- CreateAssay5Object(counts = counts_mat)
  new_assay <- SetAssayData(new_assay, layer = "data", new.data = data_mat)
  obj[["RNA"]] <- new_assay
  DefaultAssay(obj) <- "RNA"
  seurat_list[[basename(f)]] <- obj
}

# ------------------------------------------------------------------------------
# Downstream Seurat analysis
# ------------------------------------------------------------------------------
cat("🛠 Merging Seurat Objects...\n")
if (length(seurat_list) == 0) stop("❌ Input list is empty.")
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
merged <- IntegrateLayers(object = merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

cat("🤝 Finding neighbors and clustering...\n")
merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:10)
merged <- FindClusters(merged, resolution = 0.5)

cat("🗺 Running UMAP...\n")
merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:10)

cat(paste0("🔢 Number of cells in integrated object: ", ncol(merged), "\n"))
cat("✅ Integration and analysis complete\n")

# ------------------------------------------------------------------------------
# Save integrated object
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, "integrated_analysis_1.h5Seurat")
cat(paste0("💾 Saving integrated object to: ", save_path, "\n"))
SaveH5Seurat(merged, filename = save_path, overwrite = TRUE)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write.csv(as.data.frame(Idents(merged)), file = file.path(output_dir, paste0("integrated_cluster_ids_1.csv")))
cat("✅ Cluster IDs saved\n")

png(file.path(output_dir, "integrated_umap_cell_cycle_reg_jeyoon_1.png"), width = 800, height = 600)
DimPlot(merged, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

script_end_time <- Sys.time()
elapsed <- difftime(script_end_time, script_start_time, units = "mins")
cat(paste0("⏱️ Completed in ", round(elapsed, 2), " minutes\n")) 