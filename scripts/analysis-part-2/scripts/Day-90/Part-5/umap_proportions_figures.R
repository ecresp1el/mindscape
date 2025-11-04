#!/usr/bin/env Rscript

# ==============================================================================
# MindScape Step 5 of 5 – UMAP and Cell-Type Proportion Figures (Day 90)
# ==============================================================================

# ------------------------------------------------------------------------------
# INPUTS / OUTPUTS
# ------------------------------------------------------------------------------
# INPUTS:
#   - Environment variable DEG_INPUT:
#       Path to clustered Seurat object (.rds) containing UMAP reduction and cluster metadata.
#   - Environment variable MINDSCAPE_OUTPUT_DIR:
#       Output directory for saving the verified Seurat object and plots.
#
# OUTPUTS:
#   - clustered_day90_with_cluster_names_2.rds:
#       Verified Seurat object with appended cluster number–name mappings.
#   - day90_res_0.5_3_combined.{png,svg,pdf}:
#       Combined UMAP + cluster proportion figure in multiple formats.
#   - Console log messages:
#       Detailed verification and completion output confirming all steps.
# ------------------------------------------------------------------------------

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
# Cluster identities and palette (Day 90)
# ------------------------------------------------------------------------------

cluster_identities <- c(
  "0" = "MGE Striatal/GP Fated",
  "1" = "SST+, NPY +, Cortical Fated",
  "2" = "CRABP1+/PV Precursors",
  "3" = "PV precursors/Migrating cells/Cortical-fated",
  "4" = "Pre-Astrocytes/Astrocytes 1",
  "5" = "LHX8+ vMGE GABergic Striatal/GP fated 1",
  "6" = "Stressed Cells",
  "7" = "Stressed Cells",
  "8" = "LHX8+ vMGE GABergic Striatal/GP fated 2",
  "9" = "Pre-OPCs/OPCs",
  "10" = "Pre-Astrocytes/Astrocytes 2",
  "11" = "PV Precursors",
  "12" = "Dividing cells"
)

cluster_order <- as.character(0:12)

palette_colors <- brewer.pal(12, "Paired")
palette_colors <- c(palette_colors, "#8B0000")  # 13th distinct color

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

# Convert cluster IDs to characters
seu$cluster_id <- as.character(seu$seurat_clusters)

# Safely append non-intrusive new column
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

save_path <- file.path(output_dir, "clustered_day90_with_cluster_names_2.rds")
saveRDS(seu, save_path)
cat("💾 Verified Seurat object saved to:\n  ", save_path, "\n")

# ------------------------------------------------------------------------------
# UMAP verification and plotting
# ------------------------------------------------------------------------------

if (!"umap" %in% names(seu@reductions)) {
  stop("❌ No UMAP reduction found in the Seurat object.")
}
cat("✅ UMAP reduction verified.\n")

# ------------------------------------------------------------------------------
# Plot 1 – UMAP with clusters labeled by number, legend by cell-type name
# ------------------------------------------------------------------------------

umap_plot <- DimPlot(
  seu, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, label.size = 4, repel = TRUE, order = cluster_order
) +
  scale_color_manual(values = palette_colors, labels = cluster_identities) +
  ggtitle("Day 90 – UMAP by Cluster (res = 0.5)") +
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
    title = "Cluster Proportions by Sample (Day 90 res 0.5)",
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

out_prefix <- file.path(output_dir, "day90_res_0.5_3")

ggsave(paste0(out_prefix, "_combined.png"), combined, width = 9, height = 11, dpi = 300)
ggsave(paste0(out_prefix, "_combined.svg"), combined, width = 9, height = 11)
ggsave(paste0(out_prefix, "_combined.pdf"), combined, width = 9, height = 11)

cat("✅ All figures saved successfully.\n")
cat("🎉 UMAP + Proportion visualization for Day 90 complete.\n")
