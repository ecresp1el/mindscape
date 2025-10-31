#!/usr/bin/env Rscript
# ==============================================================================
# CompareClusterPvals_DirectThreshold_Summary_WithSharedProp.R
#
# PURPOSE:
#   To compare p-values and utilize this comparison to compare scRNA-seq datasets. 
#   The *proportion of shared genes* between the clusters must
#   exceed a minimum cutoff for the mapping to be considered valid.
#
#   This ensures that even if all shared genes meet p-value thresholds, the
#   match will only be accepted if enough genes overlap between clusters.
#
# INPUT (env vars):
#   REFERENCE_DEGS  - CSV with (gene, cluster, p_val_adj OR p_val)
#   NEW_DEGS        - CSV with (gene, cluster, p_val_adj OR p_val)
#   OUT_DIR         - output folder
#
# PARAMETERS (env vars; defaults):
#   PVALUE_DIFF_THRESH   (default 0.05)  - max |p_ref - p_new| allowed per gene
#   MIN_SHARED           (default 10)    - minimum shared genes for comparison
#   MIN_SHARED_PROP      (default 1.0)   - minimum proportion of genes shared
#                                          between clusters (1.0 = 100%)
#
# OUTPUT FILES:
#   - pairwise_pval_similarity.csv   : all pairwise stats
#   - matched_clusters.csv           : pairs passing thresholds + shared prop
#   - cluster_match_summary.csv      : final summary with unmatched clusters
# ==============================================================================
# This script should be run whenever it is necessary to determine if runs of an
# analysis or other sets of clusters are similar.  


# ------------------------------- packages ------------------------------------
suppressPackageStartupMessages({
  library(dplyr)   # For data manipulation: filter, mutate, select, arrange, inner_join
  library(readr)   # For reading/writing CSVs efficiently
  library(tidyr)   # For data tidying (not heavily used here but helpful for future expansions)
})

# ------------------------------- helpers -------------------------------------
# Function to identify which p-value column to use in the dataset
get_pcol <- function(df) {
  if ("p_val_adj" %in% colnames(df)) return("p_val_adj")  # Prefer adjusted p-values
  if ("p_val" %in% colnames(df)) return("p_val")          # Fallback to raw p-values
  stop("CSV must contain either 'p_val' or 'p_val_adj'.")  # Stop if neither exists
}

# Safe CSV reading function with existence check
safe_read <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)  # Prevent silent errors
  read_csv(path, show_col_types = FALSE)                  # Read CSV without printing column types
}

# ------------------------------- env vars ------------------------------------
# Read environment variables (input/output paths)
REF_PATH <- Sys.getenv("REFERENCE_DEGS")   # Reference DEGs CSV path
NEW_PATH <- Sys.getenv("NEW_DEGS")         # New DEGs CSV path
OUT_DIR  <- Sys.getenv("OUT_DIR", "./cluster_pval_match_output")  # Default output folder

# Read numeric thresholds from env vars or use defaults
PVALUE_DIFF_THRESH <- as.numeric(Sys.getenv("PVALUE_DIFF_THRESH", "0.05"))
MIN_SHARED <- as.integer(Sys.getenv("MIN_SHARED", "10"))
MIN_SHARED_PROP <- as.numeric(Sys.getenv("MIN_SHARED_PROP", "1.0"))

# Ensure required input paths are set
if (REF_PATH == "" || NEW_PATH == "") stop("Set REFERENCE_DEGS and NEW_DEGS.")

# Create output directory if it doesn't exist
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Print input/output info and thresholds for logging
cat("📥 Reference:", REF_PATH, "\n")
cat("📥 New:", NEW_PATH, "\n")
cat("📤 Output:", OUT_DIR, "\n")
cat("🔧 PVALUE_DIFF_THRESH =", PVALUE_DIFF_THRESH,
    "| MIN_SHARED =", MIN_SHARED,
    "| MIN_SHARED_PROP =", MIN_SHARED_PROP, "\n\n")

# ------------------------------- load data -----------------------------------
# Load CSVs safely
ref_df <- safe_read(REF_PATH)
new_df <- safe_read(NEW_PATH)

# Determine which p-value column to use
ref_pcol <- get_pcol(ref_df)
new_pcol <- get_pcol(new_df)

# Standardize p-value column to "p_val" and ensure cluster is character
ref_df <- ref_df %>% mutate(p_val = !!sym(ref_pcol), cluster = as.character(cluster))
new_df <- new_df %>% mutate(p_val = !!sym(new_pcol), cluster = as.character(cluster))

# Get sorted unique clusters
ref_clusters <- sort(unique(ref_df$cluster))
new_clusters <- sort(unique(new_df$cluster))

# Log clusters present in each dataset
cat("🔹 Reference clusters:", paste(ref_clusters, collapse = ", "), "\n")
cat("🔹 New clusters:", paste(new_clusters, collapse = ", "), "\n\n")

# ------------------------------- comparison ----------------------------------
# Initialize list to store pairwise comparison results
pair_list <- list()

# Loop over all reference clusters
for (r in ref_clusters) {
  # Subset reference DEGs for this cluster
  rsub <- ref_df %>% filter(cluster == r) %>% select(gene, p_val)
  
  # Loop over all new clusters
  for (n in new_clusters) {
    # Subset new DEGs for this cluster
    nsub <- new_df %>% filter(cluster == n) %>% select(gene, p_val)

    # Identify shared genes between the two clusters
    shared <- intersect(rsub$gene, nsub$gene)
    n_shared <- length(shared)
    
    # Skip if number of shared genes < MIN_SHARED
    if (n_shared < MIN_SHARED) next

    # --- NEW: shared proportion condition ------------------------------------
    # Calculate proportion of shared genes relative to total unique genes in both clusters
    total_genes <- length(unique(c(rsub$gene, nsub$gene)))
    shared_prop <- n_shared / total_genes
    
    # Skip if shared proportion < MIN_SHARED_PROP
    if (shared_prop < MIN_SHARED_PROP) next
    # -------------------------------------------------------------------------

    # Merge reference and new cluster p-values for only shared genes
    merged <- inner_join(
      rsub %>% filter(gene %in% shared),
      nsub %>% filter(gene %in% shared),
      by = "gene", suffix = c("_ref", "_new")
    )

    # Calculate absolute p-value difference and whether it is within threshold
    merged <- merged %>%
      mutate(abs_diff = abs(p_val_ref - p_val_new),
             within_thresh = abs_diff <= PVALUE_DIFF_THRESH)

    # Compute fraction of shared genes within threshold
    prop_within <- mean(merged$within_thresh, na.rm = TRUE)
    # Determine if all shared genes meet threshold
    all_within <- all(merged$within_thresh, na.rm = TRUE)

    # Store result in the pair_list
    pair_list[[length(pair_list) + 1]] <- tibble(
      ref_cluster = r,
      new_cluster = n,
      n_shared = n_shared,
      shared_prop = shared_prop,
      prop_within_thresh = prop_within,
      all_within_thresh = all_within,
      mean_abs_diff = mean(merged$abs_diff, na.rm = TRUE),
      median_abs_diff = median(merged$abs_diff, na.rm = TRUE)
    )
  }
}

# Combine all pairwise comparisons into a single dataframe
pairwise_df <- bind_rows(pair_list)

# Save pairwise results
pairwise_csv <- file.path(OUT_DIR, "pairwise_pval_similarity.csv")
write.csv(pairwise_df, pairwise_csv, row.names = FALSE)
cat("💾 Wrote pairwise results to:", pairwise_csv, "\n")

# ------------------------------- matches -------------------------------------
# Filter only pairs where all shared genes are within threshold and proportion of shared genes meets minimum
matched_df <- pairwise_df %>%
  filter(prop_within_thresh >= 1.0, shared_prop >= MIN_SHARED_PROP) %>%
  arrange(as.numeric(ref_cluster), as.numeric(new_cluster))

# Save matched clusters
matched_csv <- file.path(OUT_DIR, "matched_clusters.csv")
write.csv(matched_df, matched_csv, row.names = FALSE)
cat("💾 Wrote >=100%-similar + prop-shared matches to:", matched_csv, "\n")

# ------------------------------- summary -------------------------------------
# Determine which clusters were matched
ref_with_match <- unique(matched_df$ref_cluster)
new_with_match <- unique(matched_df$new_cluster)

# Identify unmatched clusters
ref_unmatched <- setdiff(ref_clusters, ref_with_match)
new_unmatched <- setdiff(new_clusters, new_with_match)

# Create summary tibbles for matched and unmatched clusters
summary_list <- list(
  tibble(
    type = "Matched Pair",
    ref_cluster = as.character(matched_df$ref_cluster),
    new_cluster = as.character(matched_df$new_cluster),
    n_shared = as.integer(matched_df$n_shared),
    shared_prop = as.numeric(matched_df$shared_prop)
  ),
  tibble(
    type = "Unmatched Reference",
    ref_cluster = as.character(ref_unmatched),
    new_cluster = NA_character_,
    n_shared = NA_integer_,
    shared_prop = NA_real_
  ),
  tibble(
    type = "Unmatched New",
    ref_cluster = NA_character_,
    new_cluster = as.character(new_unmatched),
    n_shared = NA_integer_,
    shared_prop = NA_real_
  )
)

# Combine into a single summary dataframe
summary_df <- bind_rows(summary_list)

# Save summary
summary_csv <- file.path(OUT_DIR, "cluster_match_summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)

# Print summary info
cat("\n📊 Summary:\n")
cat("  ✔ Matched pairs:", nrow(matched_df), "\n")
cat("  ❌ Unmatched reference clusters:", paste(ref_unmatched, collapse = ", "), "\n")
cat("  ❌ Unmatched new clusters:", paste(new_unmatched, collapse = ", "), "\n")
cat("💾 Wrote summary to:", summary_csv, "\n")
cat("✅ Completed p-value similarity comparison (with shared prop constraint).\n")
