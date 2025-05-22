# mindscape_test_real_data.R

#add print statements to the script
print("✅ Starting Seurat test script for real data")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
print("✅ Required libraries loaded successfully")

print("✅ Setting up environment variables")
# Define output directory
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR", unset = "/nfs/turbo/umms-parent/Manny_test/mindscape_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
print(paste0("✅ Output directory set to: ", output_dir))

# Path to one sample's matrix output
sample_id <- "10496-MW-1"
data_dir <- file.path("/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs", sample_id, "count")

cat(paste0("✅ Reading 10X matrix from ", data_dir, "\n"))

# Load matrix
counts <- Read10X(data.dir = data_dir)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)

# Standard pipeline
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Save cluster IDs
write.csv(as.data.frame(Idents(seurat_obj)), file = file.path(output_dir, paste0(sample_id, "_cluster_ids.csv")))
cat("✅ Cluster IDs saved\n")

# Save UMAP plot
png(file.path(output_dir, paste0(sample_id, "_umap.png")), width = 800, height = 600)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved\n")

# Export to h5ad
SaveH5Seurat(seurat_obj, filename = file.path(output_dir, paste0(sample_id, ".h5Seurat")))
Convert(file.path(output_dir, paste0(sample_id, ".h5Seurat")), dest = "h5ad")
cat("✅ Exported to .h5ad for Scanpy\n")