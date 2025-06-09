#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# Fallback mode for testing outside Snakemake
if (!exists("snakemake")) {
  message("Running in debug mode (no snakemake object detected)")

  matrix_dir <- "/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs/10496-MW-1/count/sample_filtered_feature_bc_matrix"
  output_rds <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/10496-MW-1_debug.rds"
  sample_id <- "debug"
} else {
  matrix_dir <- snakemake@input[["matrix_dir"]]
  output_rds <- snakemake@output[["rds"]]
  sample_id <- snakemake@wildcards[["sample_id"]]
}

message("Reading 10X matrix from: ", matrix_dir)
data <- Read10X(data.dir = matrix_dir)

seurat_obj <- CreateSeuratObject(counts = data, project = sample_id, assay = "RNA")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Save QC RidgePlot before filtering
plot_dir <- file.path(dirname(output_rds), "RidgePlots_QC_nFeature_nCounts_percent_Mito_sample")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
plot_file <- file.path(plot_dir, paste0(sample_id, "_RidgePlot_nFeature_nCounts_percentMito.png"))
png(filename = plot_file, width = 1600, height = 600)
print(
  RidgePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    ggtitle(paste0("QC Metrics: ", sample_id), subtitle = "RidgePlot of nFeature_RNA, nCount_RNA, percent.mt")
)
dev.off()

# Apply QC subsetting
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 500 &
           nFeature_RNA < 5000 &
           nCount_RNA < 15000 &
           percent.mt < 5
)

# Save to output
out_dir <- dirname(output_rds)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(seurat_obj, file = output_rds)
message("✅ Done: RDS file created for ", sample_id)
