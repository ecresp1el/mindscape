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
  # Conversion to legacy h5Seurat not needed if files are already .h5Seurat
  obj <- LoadH5Seurat(f)
  # Normalization assumed already done, but if you want:
  # obj <- NormalizeData(obj)
  seurat_list[[basename(f)]] <- obj
}

# ------------------------------------------------------------------------------
# Integration using Seurat standard workflow
# ------------------------------------------------------------------------------

cat("🛠 Preparing objects (Find Variable Features)...\n")
seurat_list <- lapply(seurat_list, function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  return(x)
})

cat("🔗 Finding integration anchors...\n")
anchors <- FindIntegrationAnchors(object.list = seurat_list)

cat("🧬 Integrating data...\n")
integrated <- IntegrateData(anchorset = anchors)

# ------------------------------------------------------------------------------
# Downstream Seurat analysis
# ------------------------------------------------------------------------------
DefaultAssay(integrated) <- "integrated"

cat("✨ Finding variable features on integrated data...\n")
integrated <- FindVariableFeatures(integrated, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(integrated)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cat("🔄 Cell Cycle Scoring...\n")
integrated <- CellCycleScoring(integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cat ("🔄 Cell Cycle Scoring Difference...\n")
integrated$CC.Difference <- integrated$S.Score - integrated$G2M.Score

cat("🔄 Scaling data...\n")
integrated <- ScaleData(integrated, vars.to.regress = "CC.Difference", features = all.genes)

cat("🎯 Running PCA...\n")
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))

cat("🤝 Finding neighbors and clustering...\n")
integrated <- FindNeighbors(integrated, dims = 1:10)
integrated <- FindClusters(integrated, resolution = 0.5)

cat("🗺 Running UMAP...\n")
integrated <- RunUMAP(integrated, dims = 1:10)

cat(paste0("🔢 Number of cells in integrated object: ", ncol(integrated), "\n"))
cat("✅ Integration and analysis complete\n")

# ------------------------------------------------------------------------------
# Save integrated object
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, "integrated_analysis.h5Seurat")
cat(paste0("💾 Saving integrated object to: ", save_path, "\n"))
SaveH5Seurat(integrated, filename = save_path, overwrite = TRUE)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write.csv(as.data.frame(Idents(integrated)), file = file.path(output_dir, paste0("integrated_cluster_ids.csv")))
cat("✅ Cluster IDs saved\n")

png(file.path(output_dir, "integrated_umap_cell_cycle_ref_only_second.png"), width = 800, height = 600)
DimPlot(integrated, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

script_end_time <- Sys.time()
elapsed <- difftime(script_end_time, script_start_time, units = "mins")
cat(paste0("⏱️ Completed in ", round(elapsed, 2), " minutes\n"))
