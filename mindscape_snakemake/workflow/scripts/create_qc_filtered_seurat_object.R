#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# Debug mode if not run by Snakemake
if (!exists("snakemake")) {
  message("Running in debug mode (no snakemake object detected)")
  sample_id <- "10496-MW-1"
  input_rds <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/10496-MW-1.rds"
  output_rds <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/qc_filtered/10496-MW-1_qc_filtered.rds"
  plot_file <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/RidgePlots_QC_nFeature_nCounts_percent_Mito_sample/10496-MW-1_QC_ridgeplot.png"
} else {
  sample_id <- snakemake@params[["sample"]]
  input_rds <- snakemake@input[["raw_rds"]]
  output_rds <- snakemake@output[["qc_rds"]]
  plot_file <- snakemake@output[["ridgeplot"]]
}

message("📥 Reading input Seurat RDS: ", input_rds)
seurat_obj <- readRDS(input_rds)

# RidgePlot before filtering
plot_dir <- dirname(plot_file)
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

png(filename = plot_file, width = 1600, height = 600)
print(
  RidgePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    ggtitle(paste("QC Metrics for Sample:", sample_id), subtitle = "nFeature_RNA / nCount_RNA / percent.mt")
)
dev.off()

# Apply filtering thresholds
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 500 &
           nFeature_RNA < 5000 &
           nCount_RNA < 15000 &
           percent.mt < 5
)

# Save filtered Seurat object
qc_dir <- dirname(output_rds)
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE)
}

saveRDS(seurat_obj, file = output_rds)
message("✅ Saved QC-filtered Seurat object for sample: ", sample_id)
