#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Fallback mode for testing outside Snakemake
if (!exists("snakemake")) {
  message("Running in debug mode (no snakemake object detected)")

  matrix_dir <- "/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs/10496-MW-1/count/sample_filtered_feature_bc_matrix"
  output_rds <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/10496-MW-1_debug.rds"
} else {
  matrix_dir <- snakemake@input[["matrix_dir"]]
  output_rds <- snakemake@output[["rds"]]
}

message("Reading 10X matrix from: ", matrix_dir)
message("Saving Seurat object to: ", output_rds)

# Create output directory if needed
out_dir <- dirname(output_rds)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load, create, and save Seurat object
data <- Read10X(data.dir = matrix_dir)
seurat_obj <- CreateSeuratObject(counts = data)
saveRDS(seurat_obj, file = output_rds)

message("✅ Done: RDS file created.")