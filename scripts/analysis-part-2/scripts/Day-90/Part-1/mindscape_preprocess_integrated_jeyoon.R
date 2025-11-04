#!/usr/bin/env Rscript
# mindscape_preprocess_integrated_jeyoon_day_30.R
# -------------------------------------------------------------------------------
# PURPOSE:
# Load normalized Seurat objects (.rds), merge, regress out cell cycle,
# scale, find variable features, run PCA, integrate (CCA),
# then save checkpoint RDS with verification (jhelper functions included).
# This is step 1 of 5 steps in a pipeline to generate a UMAP and cell-type proportions figure for a day 90 timepoint. This step focuses on preprocessing and integrating data. 
# -------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# INPUTS / OUTPUTS
# ------------------------------------------------------------------------------
# INPUTS:
#   - Environment variable MINDSCAPE_INPUT_DIR:
#       Directory containing normalized Seurat .rds files (per sample)
#   - Each file must be a valid Seurat v5 object
#
# OUTPUTS:
#   - Environment variable MINDSCAPE_OUTPUT_DIR:
#       Target directory for integrated Seurat object
#   - Output file:
#       preprocessed_integrated_jeyoon_test_3.rds
#   - Verification log messages confirming round-trip integrity
# ------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting preprocessing + (CCA) integration\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

# ---------------------------
# jhelper functions
# ---------------------------
jhelper_check_seurat <- function(min_version = "5.0.0") {
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat not installed in this R environment.")
  }
  if (packageVersion("Seurat") < package_version(min_version)) {
    stop(paste0("Seurat >= ", min_version, " required. Found: ", packageVersion("Seurat")))
  }
  message("Seurat version: ", packageVersion("Seurat"))
  invisible(TRUE)
}

jhelper_assert_dir_exists <- function(path) {
  if (!dir.exists(path)) stop(paste0("Directory not found: ", path))
  invisible(TRUE)
}

jhelper_verify_rds_roundtrip <- function(obj, save_path) {
  if (!file.exists(save_path)) stop("Saved file not found for verification.")
  reloaded <- readRDS(save_path)
  if (identical(obj, reloaded)) {
    cat("✅ Verification passed: reloaded identical to saved\n")
  } else {
    stop("❌ Verification failed: objects differ after reload")
  }
  invisible(TRUE)
}

# ---------------------------
# Env and inputs
# ---------------------------
input_base <- Sys.getenv("MINDSCAPE_INPUT_DIR")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")
if (input_base == "" || output_dir == "") {
  stop("❌ MINDSCAPE_INPUT_DIR or MINDSCAPE_OUTPUT_DIR not set.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# check packages / versions
jhelper_check_seurat()

cat(paste0("📂 Input directory: ", input_base, "\n"))
cat(paste0("📁 Output directory: ", output_dir, "\n"))

# ---------------------------
# Find sample .rds files (specific samples)
# ---------------------------
rds_files <- list.files(input_base, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(rds_files) == 0) stop("❌ No .rds files found in input directory.")
cat(paste0("🔎 Found ", length(rds_files), " .rds files\n"))
# If you want to restrict to particular sample names, modify below (example commented)
# sample_pattern <- "9853-MW-(1|2|3)\\.rds$"
# rds_files <- grep(sample_pattern, rds_files, value = TRUE)

# ---------------------------
# Load Seurat objects
# ---------------------------
seurat_list <- lapply(rds_files, function(f) {
  cat("📥 Loading:", f, "\n")
  obj <- readRDS(f)
  if (!inherits(obj, "Seurat")) stop(paste0("Not a Seurat object: ", f))
  obj
})
names(seurat_list) <- basename(rds_files)

# ---------------------------
# Merge + preprocess
# ---------------------------
cat("🛠 Merging Seurat Objects...\n")
merged <- Reduce(function(x, y) merge(x, y), seurat_list)

cat("🧬 Re-normalizing (safe) ...\n")
merged <- NormalizeData(merged)

# Use parallel for scaling, but be conservative
plan("multisession", workers = 8)
options(future.globals.maxSize = 50 * 1024^3)

cat("🔄 Scaling data (regressing CC.Difference)...\n")
if (!"CC.Difference" %in% colnames(merged@meta.data)) {
  warning("CC.Difference not found in metadata — ScaleData will run without regression.")
  merged <- ScaleData(merged, features = rownames(merged))
} else {
  merged <- ScaleData(merged, vars.to.regress = "CC.Difference", features = rownames(merged))
}
plan("sequential")

cat("✨ Finding variable features...\n")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 4000)

cat("🎯 Running PCA...\n")
merged <- RunPCA(merged, features = VariableFeatures(merged))

cat("🤝 Integrating layers with CCA...\n")
# IntegrateLayers + CCAIntegration are Seurat v5 constructs
merged <- IntegrateLayers(
  object = merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# ---------------------------
# Save checkpoint + verify
# ---------------------------
save_path <- file.path(output_dir, "preprocessed_integrated_jeyoon_test_3.rds")
cat("💾 Saving checkpoint:", save_path, "\n")
saveRDS(merged, save_path)

# run verification
jhelper_verify_rds_roundtrip(merged, save_path)

cat("✅ Preprocessing + integration complete\n")
cat("⏱ Total time: ", round(difftime(Sys.time(), script_start_time, units = "mins"), 2), " mins\n")
