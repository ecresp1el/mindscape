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
# Define the output directory for storing plots, tables, and converted data
# ------------------------------------------------------------------------------
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR", unset = "/nfs/turbo/umms-parent/Manny_test/mindscape_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
print(paste0("✅ Output directory set to: ", output_dir))

# ------------------------------------------------------------------------------
# Define the path to the Cell Ranger 'count' output directory for the sample
# ------------------------------------------------------------------------------
sample_id <- "10496-MW-1"
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
# Export Seurat object to h5Seurat and h5ad formats for Python/Scanpy compatibility
# ------------------------------------------------------------------------------
SaveH5Seurat(seurat_obj, filename = file.path(output_dir, paste0(sample_id, ".h5Seurat")))
Convert(file.path(output_dir, paste0(sample_id, ".h5Seurat")), dest = "h5ad")
cat("✅ Exported to .h5ad for Scanpy\n")