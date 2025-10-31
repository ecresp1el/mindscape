#!/usr/bin/env Rscript

# ==============================================================================
# MindScape Step 5 of 5 – UMAP and Cell-Type Proportion Figures (Day 30)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Environment setup
# ------------------------------------------------------------------------------

input_rds  <- Sys.getenv("DEG_INPUT")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")

if (input_rds == "" || output_dir == "") {
  stop("❌ Environment variables DEG_INPUT and/or MINDSCAPE_OUTPUT_DIR are not set.")
}

cat("📂 Input file: ", input_rds, "\n")
cat("💾 Output directory: ", output_dir, "\n")

# ------------------------------------------------------------------------------
# Load Seurat object
# ------------------------------------------------------------------------------

cat("📖 Loading Seurat object...\n")
seu <- readRDS(input_rds)
cat("✅ Loaded successfully.\n")

# ------------------------------------------------------------------------------
# Cluster identities and palette (Day 30)
# ------------------------------------------------------------------------------

cluster_identities <- c(
  "0" = "Radial Glia",
  "1" = "PV Neuron Precursors",
  "2" = "SST + cIN",
  "3" = "MGE Subpallial Neurons",
  "4" = "Radial Glia",
  "5" = "Inhibitory Progenitors (Not in Carmen/Miranda's)",
  "6" = "Inhibitory Progenitors",
  "7" = "OPCs (Not in Carmen/Miranda's)"
)

cluster_order <- as.character(c(0, 4, 6, 2, 1, 3, 5, 7))

palette_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(cluster_identities)),
  names(cluster_identities)
)

# ------------------------------------------------------------------------------
# Ensure Seurat metadata is consistent
# ------------------------------------------------------------------------------

if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("❌ 'seurat_clusters' not found in metadata. Clustering must be performed first.")
}
Idents(seu) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# Append new metadata column linking cluster number and identity
# ------------------------------------------------------------------------------

cat("\n=== Adding new metadata column 'cluster_number_name' ===\n")

seu$cluster_id <- as.character(seu$seurat_clusters)
seu$cluster_number_name <- paste0(
  seu$cluster_id, " - ", cluster_identities[seu$cluster_id]
)

# ------------------------------------------------------------------------------
# Verification of metadata append
# ------------------------------------------------------------------------------

cat("\n=== Verifying metadata append ===\n")

meta_cols <- colnames(seu@meta.data)
if (!"cluster_number_name" %in% meta_cols) {
  stop("❌ 'cluster_number_name' column not found in metadata.")
} else {
  cat("✅ Column 'cluster_number_name' successfully added.\n")
}

if (any(is.na(seu$cluster_number_name)) || any(seu$cluster_number_name == "")) {
  stop("❌ Some entries in 'cluster_number_name' are missing or blank.")
} else {
  cat("✅ All entries in 'cluster_number_name' are non-empty.\n")
}

# ------------------------------------------------------------------------------
# Entry-by-entry validation
# ------------------------------------------------------------------------------

cat("\n=== Verifying each entry format in 'cluster_number_name' ===\n")

expected_values <- paste0(names(cluster_identities), " - ", cluster_identities)
invalid_entries <- which(!seu$cluster_number_name %in% expected_values)

if (length(invalid_entries) > 0) {
  cat("❌ Detected invalid entries in 'cluster_number_name':\n")
  print(head(seu$cluster_number_name[invalid_entries], 10))
  stop(paste0("❌ Verification failed: ", length(invalid_entries), " invalid entries found."))
} else {
  cat("✅ Every entry in 'cluster_number_name' correctly matches expected mapping.\n")
}

cat("\n--- Mapping Summary ---\n")
print(table(seu$cluster_number_name, useNA = "ifany"))
cat("---------------------------------\n")

# ------------------------------------------------------------------------------
# Save verified Seurat object
# ------------------------------------------------------------------------------

save_path <- file.path(output_dir, "clustered_day30_with_cluster_names_1.rds")
saveRDS(seu, save_path)
cat("💾 Verified Seurat object saved to:\n  ", save_path, "\n")

# ------------------------------------------------------------------------------
# Detect UMAP reduction
# ------------------------------------------------------------------------------

cat("\n=== Detecting UMAP reduction ===\n")

if ("umap" %in% names(seu@reductions)) {
  umap_reduction <- "umap"
} else if ("umap.integrated.cca" %in% names(seu@reductions)) {
  umap_reduction <- "umap.integrated.cca"
} else if (all(c("umap_1", "umap_2") %in% colnames(seu[[]]))) {
  umap_reduction <- "umap"
} else {
  stop("❌ No UMAP reduction found in Seurat object.")
}
cat("✅ Using UMAP reduction: ", umap_reduction, "\n")

# ------------------------------------------------------------------------------
# Plot 1 – UMAP with clusters labeled by number, legend by cell-type name
# ------------------------------------------------------------------------------

umap_plot <- DimPlot(
  seu, reduction = umap_reduction, group.by = "seurat_clusters",
  label = TRUE, label.size = 4, repel = TRUE, order = cluster_order
) +
  scale_color_manual(values = palette_colors, labels = cluster_identities) +
  ggtitle("Day 30 – UMAP by Cluster (res = 0.2)") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ------------------------------------------------------------------------------
# Plot 2 – Cluster proportion bar plot per sample
# ------------------------------------------------------------------------------

if (!"orig.ident" %in% colnames(seu@meta.data)) {
  stop("❌ 'orig.ident' not found in metadata. Sample origin must exist for proportion plots.")
}

prop_df <- seu@meta.data %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  group_by(orig.ident) %>%
  mutate(prop = count / sum(count))

bar_plot <- ggplot(prop_df, aes(x = orig.ident, y = prop, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = palette_colors, labels = cluster_identities) +
  theme_bw(base_size = 12) +
  labs(
    title = "Cluster Proportions by Sample (Day 30 res 0.2)",
    x = "Sample", y = "Proportion"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ------------------------------------------------------------------------------
# Combine and save figures
# ------------------------------------------------------------------------------

combined <- umap_plot + bar_plot + plot_layout(ncol = 1, heights = c(2, 1))

out_prefix <- file.path(output_dir, "day30_res_0.2_2")

ggsave(paste0(out_prefix, "_combined.png"), combined, width = 9, height = 11, dpi = 300)
ggsave(paste0(out_prefix, "_combined.svg"), combined, width = 9, height = 11)
ggsave(paste0(out_prefix, "_combined.pdf"), combined, width = 9, height = 11)

cat("✅ All figures saved successfully.\n")
cat("🎉 UMAP + Proportion visualization for Day 30 complete.\n")
