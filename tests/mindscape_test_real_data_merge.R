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
for (file in h5_files) {
  cat(paste0("→ Reading ", file, "\n"))
  Convert(file, dest = "h5seurat", overwrite = TRUE)  # normalize format
  seurat_obj <- LoadH5Seurat(sub("\\.h5Seurat$", ".h5seurat", file))
  seurat_list[[basename(file)]] <- seurat_obj
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
merged_path <- file.path(output_path, "merged_mindscape.h5Seurat")
SaveH5Seurat(merged_seurat, filename = merged_path, overwrite = TRUE)
cat("✅ Merged Seurat object saved to: ", merged_path, "\n")