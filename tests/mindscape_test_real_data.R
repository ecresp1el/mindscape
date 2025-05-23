# ------------------------------------------------------------------------------
# Track total script runtime
# ------------------------------------------------------------------------------
script_start_time <- Sys.time()
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
  # Per-sample timing
  sample_start_time <- Sys.time()
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
  # Write the 'data' layer explicitly back to the legacy @data slot
  # While Seurat v5 uses a layer-based model, SeuratDisk currently expects
  # the traditional 'data' slot to be present when reading an object from .h5Seurat.
  # Without this, LoadH5Seurat() will fail even if the data layer exists.
  # This step ensures compatibility with merging and reloading workflows.
  # ------------------------------------------------------------------------------
  DefaultAssay(seurat_obj) <- "RNA"
  assay_data <- GetAssayData(seurat_obj, layer = "data")
  # Use SetAssayData to assign data slot safely (now that RNA is a legacy Assay object)
  # First, convert to legacy Assay class, then assign the data slot via SetAssayData
  seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
  seurat_obj <- SetAssayData(seurat_obj, slot = "data", new.data = assay_data)
  cat("✅ Seurat data layer set — RNA@data confirmed for saving\n")
  print(Layers(seurat_obj[["RNA"]]))  # Should show "counts" and "data"

  # (Already converted to legacy Assay and set data slot above)

  # ------------------------------------------------------------------------------
  # Save the Seurat object to .h5Seurat with explicit assay and layer declarations
  # This ensures SeuratDisk saves the required 'counts' and 'data' layers from the RNA assay.
  # Without this, SeuratDisk may silently omit these layers, causing LoadH5Seurat() or merging to fail.
  # ------------------------------------------------------------------------------
  SaveH5Seurat(
    seurat_obj,
    filename = file.path(output_dir, paste0(sample_id, ".h5Seurat")),
    overwrite = TRUE,
    assays = "RNA",
    layers = c("counts", "data")
  )
  cat("✅ Exported to .h5Seurat file (overwritten if already existed)\n")

  # End per-sample timing
  sample_end_time <- Sys.time()
  elapsed <- as.numeric(difftime(sample_end_time, sample_start_time, units = "secs"))
  elapsed_formatted <- sprintf("%02d:%02d:%02d", elapsed %/% 3600, (elapsed %% 3600) %/% 60, round(elapsed %% 60))
  cat(paste0("⏱️ Sample ", sample_id, " completed in ", elapsed_formatted, " (hh:mm:ss)\n"))

  # ------------------------------------------------------------------------------
  # Reload the saved h5Seurat file to verify the layers were written correctly
  # This step confirms that downstream functions like LoadH5Seurat() will work as expected
  # and prevents surprises during later merging.
  # ------------------------------------------------------------------------------
  reloaded_obj <- LoadH5Seurat(file.path(output_dir, paste0(sample_id, ".h5Seurat")))
  cat("✅ Reloaded saved file. RNA assay layers:\n")
  print(Layers(reloaded_obj[["RNA"]]))
}

# ------------------------------------------------------------------------------
# Print total script runtime
# ------------------------------------------------------------------------------
script_end_time <- Sys.time()
total_time <- as.numeric(difftime(script_end_time, script_start_time, units = "secs"))
total_formatted <- sprintf("%02d:%02d:%02d", total_time %/% 3600, (total_time %% 3600) %/% 60, round(total_time %% 60))
cat(paste0("✅ Total script runtime: ", total_formatted, " (hh:mm:ss)\n"))