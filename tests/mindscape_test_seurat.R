# mindscape_test_seurat.R

# Load libraries
library(Seurat)
library(ggplot2)

cat("✅ Starting Seurat test script\n")

# Load built-in test data
data("pbmc_small")
cat("✅ pbmc_small loaded successfully\n")

# Basic Seurat pipeline
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
pbmc_small <- FindNeighbors(pbmc_small, dims = 1:10)
pbmc_small <- FindClusters(pbmc_small, resolution = 0.5)
pbmc_small <- RunUMAP(pbmc_small, dims = 1:10)

cat("✅ Analysis pipeline completed\n")

# Save clustering results to CSV
cluster_ids <- Idents(pbmc_small)
write.csv(as.data.frame(cluster_ids), file = "pbmc_small_cluster_ids.csv")
cat("✅ Cluster IDs saved to pbmc_small_cluster_ids.csv\n")

# Save UMAP plot
png("pbmc_small_umap.png", width = 800, height = 600)
DimPlot(pbmc_small, reduction = "umap", label = TRUE)
dev.off()
cat("✅ UMAP plot saved to pbmc_small_umap.png\n")