# mindscape_test_real_data.R

#add print statements to the script
print("✅ Starting Seurat test script for real data")

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(ggplot2)
print("✅ Required libraries loaded successfully")

print("✅ Setting up environment variables")
# ------------------------------------------------------------------------------
# Dynamically identify the first 3 sample IDs from the Cell Ranger output folder
# ------------------------------------------------------------------------------
per_sample_base <- "/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs"
all_samples <- list.dirs(per_sample_base, recursive = FALSE, full.names = FALSE)
if (length(all_samples) < 3) {
  stop("❌ Fewer than 3 samples found in per_sample_outs.")
}
selected_samples <- all_samples[1:3]

for (sample_id in selected_samples) {
  print(paste0("🔄 Processing sample: ", sample_id))
  output_dir <- file.path("/nfs/turbo/umms-parent/Manny_test/mindscape_test_outputs", sample_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  print(paste0("✅ Output directory set to: ", output_dir))

  data_dir <- file.path(per_sample_base, sample_id, "count")
  cat(paste0("✅ Reading 10X matrix from ", data_dir, "\n"))

  feature_matrix_path <- file.path(data_dir, "sample_filtered_feature_bc_matrix")
  if (!dir.exists(feature_matrix_path)) {
    warning(paste0("⚠️ Skipping ", sample_id, ": missing expected feature matrix folder"))
    next
  }

  counts <- Read10X(data.dir = feature_matrix_path)

  # Create Seurat object and run the pipeline
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

  # Save outputs
  write.csv(as.data.frame(Idents(seurat_obj)), file = file.path(output_dir, paste0(sample_id, "_cluster_ids.csv")))
  cat("✅ Cluster IDs saved\n")

  png(file.path(output_dir, paste0(sample_id, "_umap.png")), width = 800, height = 600)
  DimPlot(seurat_obj, reduction = "umap", label = TRUE)
  dev.off()
  cat("✅ UMAP plot saved\n")

  # ------------------------------------------------------------------------------
  # Set the default assay and explicitly populate the RNA@data slot
  # This ensures that when the object is saved as .h5Seurat and later reloaded,
  # SeuratDisk will find the expected assay data structure. Without this,
  # loading or merging may fail due to missing 'counts' or 'data' slots.
  # ------------------------------------------------------------------------------
  DefaultAssay(seurat_obj) <- "RNA"
  # NOTE: Seurat v5 replaces slots with layers; use SetAssayData() for compatibility.
  seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = GetAssayData(seurat_obj, layer = "data"))

  # ------------------------------------------------------------------------------
  # Export Seurat object to .h5Seurat format
  # Note: We use overwrite = TRUE to allow repeated test runs to overwrite the same file.
  # This is important during development and testing but should be removed or set to FALSE
  # in production workflows to avoid accidental data loss.
  # ------------------------------------------------------------------------------
  SaveH5Seurat(
    seurat_obj,
    filename = file.path(output_dir, paste0(sample_id, ".h5Seurat")),
    overwrite = TRUE
  )
  cat("✅ Exported to .h5Seurat file (overwritten if already existed)\n")
}