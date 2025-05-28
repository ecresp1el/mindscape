# mindscape_process_sample_normalized.R
#
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
# - Normalized Seurat object (.h5Seurat) for downstream integration
#
# NOTE:
# This script converts Seurat v5 objects to legacy Assay format for compatibility
# with SeuratDisk::LoadH5Seurat and downstream merging.
# -------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting Seurat normalization for a single sample\n")

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(Seurat)
library(SeuratDisk)
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
# Minimal Seurat pipeline for normalization and inidvidual variable features
# ------------------------------------------------------------------------------
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Optional QC filtering (adjust thresholds if needed)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat(paste0("✅ QC filtering complete: ", ncol(seurat_obj), " cells retained\n"))

# Normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
cat("✅ Normalization and variable feature selection complete\n")

# ------------------------------------------------------------------------------
# Convert to legacy Seurat v4 Assay and save for SeuratDisk compatibility
# ------------------------------------------------------------------------------
DefaultAssay(seurat_obj) <- "RNA"
assay_data <- GetAssayData(seurat_obj, layer = "data")
seurat_obj[["RNA"]] <- as(seurat_obj[["RNA"]], "Assay")
seurat_obj <- SetAssayData(seurat_obj, slot = "data", new.data = assay_data)
cat("✅ Seurat data layer set — RNA@data confirmed for saving\n")
print(Layers(seurat_obj[["RNA"]]))

# ------------------------------------------------------------------------------
# Save to .h5Seurat with explicit layers
# ------------------------------------------------------------------------------
SaveH5Seurat(
  seurat_obj,
  filename = file.path(output_dir, paste0(sample_id, ".h5Seurat")),
  overwrite = TRUE,
  assays = "RNA",
  layers = c("counts", "data")
)

# ------------------------------------------------------------------------------
# Runtime reporting
# ------------------------------------------------------------------------------
script_end_time <- Sys.time()
elapsed <- as.numeric(difftime(script_end_time, script_start_time, units = "secs"))
elapsed_formatted <- sprintf("%02d:%02d:%02d", elapsed %/% 3600, (elapsed %% 3600) %/% 60, round(elapsed %% 60))
cat(paste0("⏱️ Completed in ", elapsed_formatted, " (hh:mm:ss)\n"))