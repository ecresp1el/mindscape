#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Fallback mode for testing outside Snakemake
if (!exists("snakemake")) {
  message("Running in debug mode (no snakemake object detected)")

  #matrix_dir <- "/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs/10496-MW-1/count/sample_filtered_feature_bc_matrix"
  #output_rds <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/10496-MW-1_debug.rds"
  #sample_id <- "debug"
  matrix_dir <- "/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs/10496-MW-2/count/sample_filtered_feature_bc_matrix"
  output_rds <- "/nfs/turbo/umms-parent/MindscapeProjects/10496-MW-per-sample-rds/seurat_rds/10496-MW-2_debug.rds"
  sample_id <- "10496-MW-2"

} else {
  matrix_dir <- snakemake@input[["matrix_dir"]]
  output_rds <- snakemake@output[["rds"]]

  sample_id <- NULL
  if (!is.null(snakemake@wildcards[["sample"]])) {
    sample_id <- snakemake@wildcards[["sample"]]
    message("✅ sample_id obtained from wildcards: ", sample_id)
  } else if (!is.null(snakemake@params[["sample"]])) {
    sample_id <- snakemake@params[["sample"]]
    message("✅ sample_id obtained from params: ", sample_id)
  } else {
    stop("❌ Neither snakemake@wildcards[['sample']] nor snakemake@params[['sample']] is defined.")
  }

  message("🔎 sample_id = ", sample_id)
  message("📁 matrix_dir = ", matrix_dir)
  message("💾 output_rds = ", output_rds)
}

message("Reading 10X matrix from: ", matrix_dir)

# ✅ Sanity check: matrix directory exists
if (!dir.exists(matrix_dir)) {
  stop("❌ matrix_dir does not exist: ", matrix_dir)
}
# load matrix data
data <- Read10X(data.dir = matrix_dir)

# ✅ Sanity check: matrix is non-empty
if (ncol(data) == 0) {
  stop("❌ Read10X returned an empty matrix for: ", matrix_dir)
}

message("Matrix dimensions: ", paste(dim(data), collapse = " x "))
print(head(colnames(data)))

#proceed if matrix is valid
seurat_obj <- tryCatch({
  CreateSeuratObject(counts = data, project = sample_id, assay = "RNA")
}, error = function(e) {
  stop("❌ Failed to create Seurat object: ", e$message)
})

# Save Seurat object (no QC performed here)
out_dir <- dirname(output_rds)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(seurat_obj, file = output_rds)
message("✅ Done: raw Seurat object saved for ", sample_id)