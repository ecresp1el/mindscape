#!/usr/bin/env Rscript

# mindscape_integrate_samples_cell_cycle_reg_jeyoon_day_30.R

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
# List all .h5Seurat files
h5_files <- list.files(input_base, pattern = "\\.h5Seurat$", full.names = TRUE, recursive = TRUE)

# Filter only 9853-MW-1 to 9853-MW-6 (typo, I know)
pattern <- paste0("9853-MW-(", paste(1:6, collapse = "|"), ")\\.h5Seurat$")
h5_files <- grep(pattern, h5_files, value = TRUE)

if (length(h5_files) == 0) stop("❌ No matching 9583-MW-1 through -6 .h5Seurat files found.")

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

# Parallelizing ScaleData
library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 50 * 1024^3)

cat("🔄 Scaling data...\n")
merged <- ScaleData(merged, vars.to.regress = "CC.Difference", features = rownames(merged))

# Deparallelizing following ScaleData
plan("sequential")

cat("✨ Finding variable features...\n")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 4000)

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
save_path <- file.path(output_dir, "integrated_analysis_jeyoon_day_30_try2.h5Seurat")
cat(paste0("💾 Preparing object for saving: ", save_path, "\n"))

# ✅ Ensure RNA assay has counts and data before saving
if ("RNA" %in% Assays(merged)) {
  DefaultAssay(merged) <- "RNA"

  rna_slots <- slotNames(merged[["RNA"]])

  if (!("counts" %in% rna_slots)) {
    cat("⚠️ RNA assay missing 'counts'; filling with 'data' slot\n")
    merged[["RNA"]]@counts <- GetAssayData(merged, slot = "data")
  }

  if (!("data" %in% rna_slots)) {
    cat("⚠️ RNA assay missing 'data'; filling with 'counts' slot\n")
    merged[["RNA"]]@data <- GetAssayData(merged, slot = "counts")
  }
}

cat(paste0("💾 Saving integrated object to: ", save_path, "\n"))
SaveH5Seurat(merged, filename = save_path, overwrite = TRUE)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write.csv(as.data.frame(Idents(merged)), file = file.path(output_dir, paste0("integrated_cluster_ids_jeyoon_day_30_only_6_try2.csv")))
cat("✅ Cluster IDs saved\n")

png(file.path(output_dir, "integrated_umap_cell_cycle_reg_jeyoon_day_30_only_6_try2.png"), width = 800, height = 600)
DimPlot(merged, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

script_end_time <- Sys.time()
elapsed <- difftime(script_end_time, script_start_time, units = "mins")
cat(paste0("⏱️ Completed in ", round(elapsed, 2), " minutes\n")) 