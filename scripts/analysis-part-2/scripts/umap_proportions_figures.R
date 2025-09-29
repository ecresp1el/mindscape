#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# User-defined cluster names (in the same order as seurat_clusters levels)
# ------------------------------------------------------------------------------
cluster_names <- c(
  "SST+/NPY+ Cortical Fated", 
  "Dividing cells", 
  "CRABP1+/PV Precursors", 
  "MGE Striatal/GP Fated", 
  "PV Precursors/Migrating cells/Cortical Fated", 
  "Stressed Cells", 
  "LHX8 + vMGE GABAergic Striatal/GP-fated 1", 
  "LHX8 + vMGE GABAergic Striatal/GP-fated 2", 
  "Pre-OPCs/OPCs", 
  "PV Precursors", 
  "Pre-Astrocytes/Astrocytes", 
  "vascular/ECM-producing"
)

# ------------------------------------------------------------------------------
# Retrieve environment variables
# ------------------------------------------------------------------------------
input_rds <- Sys.getenv("DEG_INPUT")
if (input_rds == "") stop("❌ DEG_INPUT is not set.")
if (!file.exists(input_rds)) stop(paste0("❌ Input file does not exist: ", input_rds))
cat(paste0("📥 Loading integrated Seurat object: ", input_rds, "\n"))

output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")
if (output_dir == "") stop("❌ MINDSCAPE_OUTPUT_DIR is not set.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("📁 Output directory: ", output_dir, "\n"))

# ------------------------------------------------------------------------------
# Load Seurat object
# ------------------------------------------------------------------------------
seu <- readRDS(input_rds)

# ------------------------------------------------------------------------------
# Detect which UMAP reduction exists
# ------------------------------------------------------------------------------
umap_reduction <- NULL
if ("umap" %in% names(seu@reductions)) {
  umap_reduction <- "umap"
} else if ("umap.integrated.cca" %in% names(seu@reductions)) {
  umap_reduction <- "umap.integrated.cca"
} else if (all(c("umap_1","umap_2") %in% colnames(seu[[]]))) {
  umap_reduction <- "umap"  # fallback if embeddings stored in metadata
} else {
  stop("❌ UMAP coordinates not found. Cannot match integration UMAP.")
}
cat(paste0("✅ Using UMAP reduction: ", umap_reduction, "\n"))

# ------------------------------------------------------------------------------
# Ensure clustering metadata exists
# ------------------------------------------------------------------------------
if (!"seurat_clusters" %in% colnames(seu[[]])) {
  stop("❌ seurat_clusters metadata not found. Cannot plot clusters.")
}
Idents(seu) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# Define color palette with cluster names
# ------------------------------------------------------------------------------
clusters <- levels(seu)
n_clusters <- length(clusters)
if (length(cluster_names) != n_clusters) {
  warning("Number of cluster names does not match number of clusters. Using default numbering.")
  cluster_names <- paste("Cluster", clusters)
}
names(cluster_names) <- clusters
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(n_clusters), clusters)

# ------------------------------------------------------------------------------
# UMAP plot
# ------------------------------------------------------------------------------
p1 <- DimPlot(seu, reduction = umap_reduction, group.by = "seurat_clusters", 
              cols = cluster_colors, label = TRUE) +
  scale_color_manual(values = cluster_colors, labels = cluster_names) +
  ggtitle("UMAP by cluster") +
  theme_classic()
ggsave(file.path(output_dir, "umap_clusters_rds_with_labels9.png"), p1, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# Bar plot of cluster proportions per sample
# ------------------------------------------------------------------------------
meta.data <- seu[[]]
sample_column <- "orig.ident"

counts <- meta.data %>%
  group_by(!!sym(sample_column), seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(!!sym(sample_column)) %>%
  mutate(proportion = count / sum(count))

p2 <- ggplot(counts, aes_string(x = sample_column, y = "proportion", fill = "seurat_clusters")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors, labels = cluster_names) +
  ylab("Proportion of cells") +
  xlab("Sample") +
  ggtitle("Cluster proportions per sample") +
  theme_classic()
ggsave(file.path(output_dir, "cluster_proportions_rds_with_labels9.png"), p2, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# Combined stacked plot
# ------------------------------------------------------------------------------
combined_plot <- p1 / p2

# Save in PNG, SVG, PDF
ggsave(file.path(output_dir, "umap_and_bar_combined_rds_with_labels9.png"),
       combined_plot, width = 7, height = 10, dpi = 300)
ggsave(file.path(output_dir, "umap_and_bar_combined_rds_with_labels9.svg"),
       combined_plot, width = 7, height = 10)
ggsave(file.path(output_dir, "umap_and_bar_combined_rds_with_labels9.pdf"),
       combined_plot, width = 7, height = 10)

cat("✅ UMAP, bar plot, and combined stacked plot (PNG, SVG, PDF) saved successfully.\n")
