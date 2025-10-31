#!/usr/bin/env Rscript

# mindscape_process_sample_normalized.R
# -------------------------------------------------------------------------------
# PURPOSE:
# This script is part of the MindScape pipeline and is specifically designed for
# analyzing ventralized S-phase organoid samples processed using 10x Genomics
# FLEX scRNA-seq chemistry. It assumes the data was demultiplexed via
# `cellranger multi` and organized under:
#
#   /nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs/<sample_id>/count/
#
# This script is intended to be executed as part of a SLURM array job. It reads
# one sample at a time (via the MINDSCAPE_SINGLE_SAMPLE environment variable),
# processes it using a standard Seurat pipeline, and writes outputs to the path
# defined in MINDSCAPE_OUTPUT_DIR.
#
# OUTPUTS (per sample):
# - Normalized Seurat object (.rds, native Seurat v5 format)
# -------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting Seurat normalization for a single sample\n")

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
})
cat("✅ Required libraries loaded successfully\n")

# ------------------------------------------------------------------------------
# Retrieve environment variables set by SLURM wrapper script
# ------------------------------------------------------------------------------
sample_id <- Sys.getenv("MINDSCAPE_SINGLE_SAMPLE")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")

if (sample_id == "" || output_dir == "") {
  stop("❌ MINDSCAPE_SINGLE_SAMPLE or MINDSCAPE_OUTPUT_DIR is not set.")
}

cat(paste0("🔄 Processing sample: ", sample_id, "\n"))
cat(paste0("📁 Output directory: ", output_dir, "\n"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load sample matrix
# ------------------------------------------------------------------------------
per_sample_base <- "/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs"
data_dir <- file.path(per_sample_base, sample_id, "count")
feature_matrix_path <- file.path(data_dir, "sample_filtered_feature_bc_matrix")

if (!dir.exists(feature_matrix_path)) {
  stop(paste0("❌ Matrix directory not found: ", feature_matrix_path))
}

cat(paste0("📥 Reading 10X matrix from: ", feature_matrix_path, "\n"))
counts <- Read10X(data.dir = feature_matrix_path)

# ------------------------------------------------------------------------------
# Minimal Seurat pipeline for normalization and variable features
# ------------------------------------------------------------------------------
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# ------------------------------------------------------------------------------
# QC filtering with strict validation
# ------------------------------------------------------------------------------
qc_min_features <- 500
qc_max_features <- 5000
qc_max_counts   <- 15000
qc_max_mt       <- 5

n_before <- ncol(seurat_obj)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > qc_min_features &
           nFeature_RNA < qc_max_features &
           nCount_RNA   < qc_max_counts &
           percent.mt   < qc_max_mt
)

n_after <- ncol(seurat_obj)

if (n_after == 0) {
  stop("❌ QC filtering removed all cells — check thresholds!")
}
if (n_after == n_before) {
  warning("⚠️ QC filtering did not remove any cells — thresholds may be too lenient.")
}

violations <- which(
  seurat_obj$nFeature_RNA <= qc_min_features |
  seurat_obj$nFeature_RNA >= qc_max_features |
  seurat_obj$nCount_RNA   >= qc_max_counts |
  seurat_obj$percent.mt   >= qc_max_mt
)

if (length(violations) > 0) {
  stop("❌ QC validation failed — some cells outside thresholds remain.")
}

cat(paste0("✅ QC filtering complete: ", n_after, " of ", n_before, " cells retained\n"))

# ------------------------------------------------------------------------------
# Normalization
# ------------------------------------------------------------------------------
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
cat("✅ Normalization complete\n")

# ------------------------------------------------------------------------------
# Cell Cycle Scoring
# ------------------------------------------------------------------------------
cat("🧬 Cell Cycle Scoring...\n")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score

# ------------------------------------------------------------------------------
# Save to native .rds (Seurat v5 format)
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, paste0(sample_id, ".rds"))
cat("💾 Saving Seurat object to:", save_path, "\n")
saveRDS(seurat_obj, file = save_path)

# ------------------------------------------------------------------------------
# Strict verification helper: reload and compare
# ------------------------------------------------------------------------------
verify_rds_integrity <- function(original_obj, saved_path) {
  reloaded <- readRDS(saved_path)
  identical_check <- identical(original_obj, reloaded)
  if (!identical_check) {
    stop("❌ Verification failed: saved and reloaded objects are NOT identical.")
  }
  cat("✅ Verification passed: saved and reloaded objects are strictly identical\n")
  invisible(TRUE)
}

verify_rds_integrity(seurat_obj, save_path)

# ------------------------------------------------------------------------------
# Runtime reporting
# ------------------------------------------------------------------------------
script_end_time <- Sys.time()
elapsed <- as.numeric(difftime(script_end_time, script_start_time, units = "secs"))
elapsed_formatted <- sprintf("%02d:%02d:%02d", elapsed %/% 3600, (elapsed %% 3600) %/% 60, round(elapsed %% 60))
cat(paste0("⏱️ Completed in ", elapsed_formatted, " (hh:mm:ss)\n"))