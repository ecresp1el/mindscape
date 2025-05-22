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
# Set a dedicated path for all MindScape test outputs
# This ensures all output files (plots, tables, h5Seurat objects) are written
# to a structured, shared location that is intentionally separated from any
# raw input data (e.g., Cell Ranger outputs). This promotes safety, clarity,
# and long-term reproducibility across different systems and users.
# ------------------------------------------------------------------------------
sample_id <- "10496-MW-1"
output_base <- "/nfs/turbo/umms-parent/Manny_test/mindscape_test_outputs"
output_dir <- file.path(output_base, sample_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
print(paste0("✅ Output directory set to: ", output_dir))

data_dir <- file.path("/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs", sample_id, "count")

cat(paste0("✅ Reading 10X matrix from ", data_dir, "\n"))

# ------------------------------------------------------------------------------
# Construct path to the actual expression matrix directory
# This should point to 'sample_filtered_feature_bc_matrix' inside the sample's count folder
# ------------------------------------------------------------------------------
feature_matrix_path <- file.path(data_dir, "sample_filtered_feature_bc_matrix")
if (!dir.exists(feature_matrix_path)) {
  stop("❌ Expected feature matrix directory not found: ", feature_matrix_path)
}

# Load the 10X gene expression matrix from the sample's filtered feature matrix directory
counts <- Read10X(data.dir = feature_matrix_path)

# ------------------------------------------------------------------------------
# Run a standard Seurat processing pipeline
# ------------------------------------------------------------------------------
# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)

# Normalize data
seurat_obj <- NormalizeData(seurat_obj)
# Identify highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj)
# Scale data
seurat_obj <- ScaleData(seurat_obj)
# Perform PCA dimensionality reduction
seurat_obj <- RunPCA(seurat_obj)
# Find nearest neighbors based on PCA
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
# Cluster cells
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# ------------------------------------------------------------------------------
# Save clustering output and visualization to output directory
# ------------------------------------------------------------------------------
# Save cluster IDs to CSV
write.csv(as.data.frame(Idents(seurat_obj)), file = file.path(output_dir, paste0(sample_id, "_cluster_ids.csv")))
cat("✅ Cluster IDs saved\n")

# Save UMAP plot as PNG
png(file.path(output_dir, paste0(sample_id, "_umap.png")), width = 800, height = 600)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

# ------------------------------------------------------------------------------
# Export Seurat object to h5Seurat format
# ------------------------------------------------------------------------------
SaveH5Seurat(seurat_obj, filename = file.path(output_dir, paste0(sample_id, ".h5Seurat")))
cat("✅ Exported to .h5Seurat file\n")