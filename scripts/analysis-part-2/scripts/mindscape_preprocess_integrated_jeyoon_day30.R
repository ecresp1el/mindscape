#!/usr/bin/env Rscript

# mindscape_preprocess_integrated_jeyoon_day_30.R
# -------------------------------------------------------------------------------
# PURPOSE:
# Load normalized Seurat objects (.rds), merge, regress out cell cycle,
# scale, find variable features, run PCA, integrate (CCA),
# then save checkpoint RDS with verification.
# -------------------------------------------------------------------------------

script_start_time <- Sys.time()
cat("✅ Starting preprocessing + integration (before clustering)\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

# ------------------------------------------------------------------------------
# Retrieve environment variables
# ------------------------------------------------------------------------------
input_base <- Sys.getenv("MINDSCAPE_INPUT_DIR")
output_dir <- Sys.getenv("MINDSCAPE_OUTPUT_DIR")

if (input_base == "" || output_dir == "") {
  stop("❌ MINDSCAPE_INPUT_DIR or MINDSCAPE_OUTPUT_DIR not set.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load .rds normalized Seurat objects
# ------------------------------------------------------------------------------
rds_files <- list.files(input_base, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
pattern <- paste0("9853-MW-(", paste(1:6, collapse = "|"), ")\\.rds$")
rds_files <- grep(pattern, rds_files, value = TRUE)

if (length(rds_files) == 0) stop("❌ No matching RDS files found.")

seurat_list <- lapply(rds_files, function(f) {
  obj <- readRDS(f)
  if (!inherits(obj, "Seurat")) stop(paste("❌ Not a Seurat object:", f))
  obj
})

# ------------------------------------------------------------------------------
# Preprocessing
# ------------------------------------------------------------------------------
cat("🛠 Merging Seurat Objects...\n")
merged <- merge(x = seurat_list[[1]], y = seurat_list[-1])

cat("🧬 Normalize again...\n")
merged <- NormalizeData(merged)

# Parallel ScaleData
plan("multisession", workers = 8)
options(future.globals.maxSize = 50 * 1024^3)
cat("🔄 Scaling data (regressing CC.Difference)...\n")
merged <- ScaleData(merged, vars.to.regress = "CC.Difference", features = rownames(merged))
plan("sequential")

cat("✨ Finding variable features...\n")
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 4000)

cat("🎯 Running PCA...\n")
merged <- RunPCA(merged, features = VariableFeatures(merged))

cat("🤝 Integrating data (CCA)...\n")
merged <- IntegrateLayers(
  object = merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# Save checkpoint with verification
# ------------------------------------------------------------------------------
verify_rds_roundtrip <- function(obj, save_path) {
  reloaded <- readRDS(save_path)
  if (identical(obj, reloaded)) {
    cat("✅ Verification passed: reloaded identical to saved\n")
  } else {
    stop("❌ Verification failed: objects differ")
  }
}

save_path <- file.path(output_dir, "preprocessed_integrated_jeyoon_day_30.rds")
cat("💾 Saving checkpoint:", save_path, "\n")
saveRDS(merged, save_path)

verify_rds_roundtrip(merged, save_path)

cat("✅ Preprocessing + integration complete\n")
cat("⏱ Total time: ", round(difftime(Sys.time(), script_start_time, units="mins"), 2), " mins\n")