#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Define cluster metadata
# ------------------------------------------------------------------------------

# Dictionary: cluster number → identity
cluster_identities <- c(
  "0"  = "Radial Glia",
  "1"  = "PV Neuron Precursors",
  "2"  = "SST + cIN",
  "3"  = "MGE Subpallial Neurons",
  "4"  = "Radial Glia",
  "5"  = "Inhibitory Progenitors (Not in Carmen/Miranda's)",
  "6"  = "Inhibitory Progenitors",
  "7"  = "OPCs (Not in Carmen/Miranda's)"
)

# Cluster order for stacked barplot (top → bottom)
cluster_order <- c(0, 4, 6, 2, 1, 3, 5, 7)        

# Generate color palette (linked by cluster number)
cluster_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(cluster_identities)),
  names(cluster_identities)
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
  umap_reduction <- "umap"
} else {
  stop("❌ UMAP coordinates not found.")
}
cat(paste0("✅ Using UMAP reduction: ", umap_reduction, "\n"))

# ------------------------------------------------------------------------------
# Ensure clustering metadata exists
# ------------------------------------------------------------------------------
if (!"seurat_clusters" %in% colnames(seu[[]])) {
  stop("❌ seurat_clusters metadata not found.")
}
Idents(seu) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# Map clusters → enforce ordering
# ------------------------------------------------------------------------------
seu$cluster_id <- as.character(seu$seurat_clusters)

# Control bar stacking order
seu$cluster_id <- factor(seu$cluster_id, levels = as.character(cluster_order))

# ------------------------------------------------------------------------------
# UMAP plot (labels = cluster numbers, legend = names)
# ------------------------------------------------------------------------------
p1 <- DimPlot(
  seu, reduction = umap_reduction, group.by = "cluster_id",
  cols = cluster_colors, label = TRUE, label.size = 4
) +
  scale_color_manual(
    values = cluster_colors,
    labels = cluster_identities[levels(seu$cluster_id)],   # legend shows names
    breaks = levels(seu$cluster_id)                       # keep number order
  ) +
  ggtitle("UMAP by cluster") +
  theme_classic() +
  guides(color = guide_legend(title = "Cell Type"))

ggsave(file.path(output_dir, "umap_clusters_day30_res_0.2.png"), p1, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# Bar plot of cluster proportions per sample
# ------------------------------------------------------------------------------
meta.data <- seu[[]]
sample_column <- "orig.ident"

counts <- meta.data %>%
  group_by(!!sym(sample_column), cluster_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(!!sym(sample_column)) %>%
  mutate(proportion = count / sum(count))

p2 <- ggplot(counts, aes_string(x = sample_column, y = "proportion", fill = "cluster_id")) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = cluster_colors,
    labels = cluster_identities[levels(seu$cluster_id)],   # legend = names
    breaks = levels(seu$cluster_id)
  ) +
  ylab("Proportion of cells") +
  xlab("Sample") +
  ggtitle("Cluster proportions per sample") +
  theme_classic() +
  guides(fill = guide_legend(title = "Cell Type"))

ggsave(file.path(output_dir, "cluster_proportions_day30_res_0.2.png"), p2, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# Combined stacked plot
# ------------------------------------------------------------------------------
combined_plot <- p1 / p2
ggsave(file.path(output_dir, "umap_and_bar_combined_day30_res_0.2.png"), combined_plot, width = 7, height = 10, dpi = 300)
ggsave(file.path(output_dir, "umap_and_bar_combined_day30_res_0.2.svg"), combined_plot, width = 7, height = 10)
ggsave(file.path(output_dir, "umap_and_bar_combined_day30_res_0.2.pdf"), combined_plot, width = 7, height = 10)

cat("✅ UMAP, bar plot, and combined stacked plot saved successfully (Day 30).\n")
