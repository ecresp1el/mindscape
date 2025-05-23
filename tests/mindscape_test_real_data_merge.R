# merge_seurat_objects.R

# ------------------------------------------------------------------------------
# This script loads multiple Seurat objects (saved as .h5Seurat files) generated
# by the MindScape test pipeline and merges them into a single combined object.
# ------------------------------------------------------------------------------

library(Seurat)
library(SeuratDisk)
library(ggplot2)

cat("✅ Starting merge script for MindScape Seurat objects\n")

# ------------------------------------------------------------------------------
# Define input and output directories
# ------------------------------------------------------------------------------
input_base <- "/nfs/turbo/umms-parent/Manny_test/mindscape_test_outputs"
output_path <- "/nfs/turbo/umms-parent/Manny_test/mindscape_test_outputs/merged"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Locate .h5Seurat files
# ------------------------------------------------------------------------------
h5_files <- list.files(input_base, pattern = "\\.h5Seurat$", recursive = TRUE, full.names = TRUE)
if (length(h5_files) < 2) {
  stop("❌ Need at least 2 .h5Seurat files to merge.")
}

cat(paste0("✅ Found ", length(h5_files), " .h5Seurat files to merge\n"))

# ------------------------------------------------------------------------------
# Load and merge
# ------------------------------------------------------------------------------
seurat_list <- list()
# Load each .h5Seurat file directly into memory
for (file in h5_files) {
  cat(paste0("→ Reading ", file, "\n"))

  # ------------------------------------------------------------------------------
  # Load each Seurat object safely with tryCatch
  # Some .h5Seurat files may be incomplete or corrupted (e.g., missing 'data' layer)
  # This block prevents the whole script from crashing and logs the failing files
  # so they can be debugged or regenerated. This is especially useful in testing
  # or mixed-quality batch datasets.
  #
  # IMPORTANT: This assumes that all input .h5Seurat files were saved using a legacy
  # Seurat v4-compatible Assay object (not Assay5), with both 'counts' and 'data' slots.
  # This is critical because SeuratDisk cannot load files that only use the Seurat v5 layered model.
  # ------------------------------------------------------------------------------
  tryCatch({
    seurat_obj <- LoadH5Seurat(file)
    seurat_list[[basename(file)]] <- seurat_obj
  }, error = function(e) {
    cat(paste0("❌ Failed to load ", file, " — ", e$message, "\n"))
  })
}

merged_seurat <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "merged_mindscape"
)

# ------------------------------------------------------------------------------
# Save merged object
# ------------------------------------------------------------------------------
# Save the merged Seurat object.
# The output file will use the Seurat v4 Assay structure (as long as input objects were converted properly),
# which ensures compatibility with downstream SeuratDisk tools and reloading workflows.
merged_path <- file.path(output_path, "merged_mindscape.h5Seurat")
SaveH5Seurat(merged_seurat, filename = merged_path, overwrite = TRUE)
cat("✅ Merged Seurat object saved to: ", merged_path, "\n")