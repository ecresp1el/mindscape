#!/usr/bin/env Rscript

# ==============================================================================
# CompareClusterPvals_DirectThreshold_UniqueMatching.R
#
# PURPOSE:
#   Pairwise compare p-values for every (ref_cluster, new_cluster).
#   Require:
#     - at least MIN_SHARED shared genes,
#     - at least MIN_SHARED_PROP proportion of genes shared (shared / union),
#     - at least MIN_PROP_WITHIN proportion of shared genes having |p_ref-p_new| <= PVALUE_DIFF_THRESH.
#   Produce a UNIQUE 1:1 mapping by maximizing prop_within (Hungarian solve_LSAP).
#
# INPUT (env vars):
#   REFERENCE_DEGS, NEW_DEGS, OUT_DIR
#
# PARAMETERS (env vars; defaults):
#   PVALUE_DIFF_THRESH  - numeric (default 1e-12) absolute p-value difference tolerance
#   MIN_SHARED          - integer (default 10)
#   MIN_SHARED_PROP     - numeric fraction (default 1.0)
#   MIN_PROP_WITHIN     - numeric fraction of shared genes required within PVALUE_DIFF_THRESH (default 1.0)
#
# OUTPUT:
#   pairwise_pval_similarity.csv
#   allowed_pairs.csv
#   hungarian_assignment_raw.csv
#   final_matches.csv
#   cluster_match_summary.csv
#   per-match pval difference CSVs for final matches
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(clue)   # solve_LSAP
  library(tibble)
})

# ------------------------------- env vars ------------------------------------

REF_PATH <- Sys.getenv("REFERENCE_DEGS")
NEW_PATH <- Sys.getenv("NEW_DEGS")
OUT_DIR  <- Sys.getenv("OUT_DIR", "./cluster_pval_match_output")

PVALUE_DIFF_THRESH <- as.numeric(Sys.getenv("PVALUE_DIFF_THRESH", "1e-12"))
MIN_SHARED <- as.integer(Sys.getenv("MIN_SHARED", "10"))
MIN_SHARED_PROP <- as.numeric(Sys.getenv("MIN_SHARED_PROP", "1.0"))
MIN_PROP_WITHIN <- as.numeric(Sys.getenv("MIN_PROP_WITHIN", "1.0"))

if (REF_PATH == "" || NEW_PATH == "") stop("Set REFERENCE_DEGS and NEW_DEGS environment variables.")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("📥 Reference:", REF_PATH, "\n")
cat("📥 New:", NEW_PATH, "\n")
cat("📤 Output:", OUT_DIR, "\n")
cat("🔧 PVALUE_DIFF_THRESH =", PVALUE_DIFF_THRESH,
    "| MIN_SHARED =", MIN_SHARED,
    "| MIN_SHARED_PROP =", MIN_SHARED_PROP,
    "| MIN_PROP_WITHIN =", MIN_PROP_WITHIN, "\n\n")

# ------------------------------- helpers -------------------------------------

get_pcol <- function(df) {
  if ("p_val_adj" %in% colnames(df)) return("p_val_adj")
  if ("p_val" %in% colnames(df)) return("p_val")
  stop("CSV must contain either 'p_val' or 'p_val_adj'.")
}
safe_read <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read_csv(path, show_col_types = FALSE)
}

# ------------------------------- load data -----------------------------------

ref_df <- safe_read(REF_PATH)
new_df <- safe_read(NEW_PATH)

ref_pcol <- get_pcol(ref_df)
new_pcol <- get_pcol(new_df)

ref_df <- ref_df %>% mutate(p_val = !!sym(ref_pcol), cluster = as.character(cluster))
new_df <- new_df %>% mutate(p_val = !!sym(new_pcol), cluster = as.character(cluster))

ref_clusters <- sort(unique(ref_df$cluster))
new_clusters <- sort(unique(new_df$cluster))

cat("🔹 Reference clusters:", paste(ref_clusters, collapse = ", "), "\n")
cat("🔹 New clusters:", paste(new_clusters, collapse = ", "), "\n\n")

# ------------------------------- pairwise ------------------------------------

pair_list <- list()

for (r in ref_clusters) {
  rsub <- ref_df %>% filter(cluster == r) %>% select(gene, p_val)
  for (n in new_clusters) {
    nsub <- new_df %>% filter(cluster == n) %>% select(gene, p_val)

    shared <- intersect(rsub$gene, nsub$gene)
    n_shared <- length(shared)
    n_ref_total <- nrow(rsub)
    n_new_total <- nrow(nsub)

    # compute union size and prop_shared
    total_genes <- length(unique(c(rsub$gene, nsub$gene)))
    prop_shared <- if (total_genes > 0) n_shared / total_genes else 0

    # Skip pairs that don't meet minimum shared count
    if (n_shared < MIN_SHARED) {
      pair_list[[length(pair_list) + 1]] <- tibble(
        ref_cluster = r, new_cluster = n, n_ref_genes = n_ref_total, n_new_genes = n_new_total,
        n_shared = n_shared, prop_shared = prop_shared,
        prop_within = NA_real_, all_within = NA, mean_abs_diff = NA_real_, median_abs_diff = NA_real_
      )
      next
    }

    # Skip pairs that don't meet shared proportion
    if (prop_shared < MIN_SHARED_PROP) {
      pair_list[[length(pair_list) + 1]] <- tibble(
        ref_cluster = r, new_cluster = n, n_ref_genes = n_ref_total, n_new_genes = n_new_total,
        n_shared = n_shared, prop_shared = prop_shared,
        prop_within = NA_real_, all_within = NA, mean_abs_diff = NA_real_, median_abs_diff = NA_real_
      )
      next
    }

    # Compare p-values for shared genes
    merged <- inner_join(
      rsub %>% filter(gene %in% shared) %>% rename(p_ref = p_val),
      nsub %>% filter(gene %in% shared) %>% rename(p_new = p_val),
      by = "gene"
    )

    merged <- merged %>%
      mutate(abs_diff = abs(p_ref - p_new),
             within = abs_diff <= PVALUE_DIFF_THRESH)

    prop_within <- mean(merged$within, na.rm = TRUE)
    all_within <- all(merged$within, na.rm = TRUE)

    pair_list[[length(pair_list) + 1]] <- tibble(
      ref_cluster = r, new_cluster = n, n_ref_genes = n_ref_total, n_new_genes = n_new_total,
      n_shared = n_shared, prop_shared = prop_shared,
      prop_within = prop_within, all_within = all_within,
      mean_abs_diff = mean(merged$abs_diff, na.rm = TRUE),
      median_abs_diff = median(merged$abs_diff, na.rm = TRUE)
    )
  }
}

pairwise_df <- bind_rows(pair_list)

pairwise_csv <- file.path(OUT_DIR, "pairwise_pval_similarity.csv")
write.csv(pairwise_df, pairwise_csv, row.names = FALSE)
cat("💾 Wrote pairwise results to:", pairwise_csv, "\n")

# -------------------------- allowed pairs ------------------------------------

allowed_df <- pairwise_df %>%
  filter(!is.na(prop_within) & n_shared >= MIN_SHARED & prop_shared >= MIN_SHARED_PROP & prop_within >= MIN_PROP_WITHIN)

allowed_csv <- file.path(OUT_DIR, "allowed_pairs.csv")
write.csv(allowed_df, allowed_csv, row.names = FALSE)
cat("💾 Wrote allowed pairs (passed thresholds) to:", allowed_csv, "\n")
cat("🔹 Allowed pair count:", nrow(allowed_df), "\n\n")

if (nrow(allowed_df) == 0) {
  cat("❌ No pairs passed thresholds — exiting.\n")
  # still create empty outputs for downstream
  write.csv(tibble(), file.path(OUT_DIR, "hungarian_assignment_raw.csv"), row.names = FALSE)
  write.csv(tibble(), file.path(OUT_DIR, "final_matches.csv"), row.names = FALSE)
  write.csv(tibble(), file.path(OUT_DIR, "cluster_match_summary.csv"), row.names = FALSE)
  quit(status = 0)
}

# -------------------------- unique 1:1 assignment -----------------------------

# Build cost matrix minimizing (1 - prop_within) for allowed pairs; disallow others with a big cost
ref_levels <- sort(unique(pairwise_df$ref_cluster))
new_levels <- sort(unique(pairwise_df$new_cluster))
cost_mat <- matrix(Inf, nrow = length(ref_levels), ncol = length(new_levels),
                   dimnames = list(ref_levels, new_levels))

# initialize with large cost
big_cost <- 1e6
cost_mat[,] <- big_cost

for (i in seq_len(nrow(allowed_df))) {
  r <- as.character(allowed_df$ref_cluster[i])
  n <- as.character(allowed_df$new_cluster[i])
  # cost: lower if prop_within higher; prop_within in [0,1], so cost in [0,1]
  cost_mat[r, n] <- 1 - allowed_df$prop_within[i]
}

# Solve Hungarian for square matrix: if different sizes, solve_LSAP expects square; expand if needed.
nrowm <- nrow(cost_mat); ncolm <- ncol(cost_mat)
if (nrowm != ncolm) {
  # make square by padding with dummy rows/cols filled with big_cost
  msize <- max(nrowm, ncolm)
  mat2 <- matrix(big_cost, nrow = msize, ncol = msize)
  rownames(mat2) <- c(rownames(cost_mat), paste0("pad_ref_", seq_len(msize - nrowm)))
  colnames(mat2) <- c(colnames(cost_mat), paste0("pad_new_", seq_len(msize - ncolm)))
  mat2[rownames(cost_mat), colnames(cost_mat)] <- cost_mat
  cost_mat_sq <- mat2
} else {
  cost_mat_sq <- cost_mat
}

assignment <- solve_LSAP(cost_mat_sq)
assigned_new_all <- colnames(cost_mat_sq)[assignment]

# Build assignment dataframe and then filter to real rows/cols (no pads)
assigned_df_full <- data.frame(
  ref_cluster = rownames(cost_mat_sq),
  new_cluster = assigned_new_all,
  cost = mapply(function(r, n) cost_mat_sq[r, n], rownames(cost_mat_sq), assigned_new_all, SIMPLIFY = TRUE),
  stringsAsFactors = FALSE
)

# Keep only assignments between real clusters (not padded)
assigned_df <- assigned_df_full %>%
  filter(grepl("^pad_ref_", ref_cluster) == FALSE & grepl("^pad_new_", new_cluster) == FALSE)

hungarian_csv <- file.path(OUT_DIR, "hungarian_assignment_raw.csv")
write.csv(assigned_df, hungarian_csv, row.names = FALSE)
cat("💾 Wrote raw Hungarian assignment to:", hungarian_csv, "\n")

# -------------------------- finalize matches ---------------------------------

# Keep only assignments that correspond to allowed pairs (i.e., cost < big_cost)
final_assigned <- assigned_df %>%
  left_join(allowed_df %>% select(ref_cluster, new_cluster, prop_within, prop_shared, n_shared), by = c("ref_cluster", "new_cluster")) %>%
  filter(!is.na(prop_within) & cost < big_cost)

final_csv <- file.path(OUT_DIR, "final_matches.csv")
write.csv(final_assigned, final_csv, row.names = FALSE)
cat("💾 Wrote final unique matches to:", final_csv, "\n")
cat("🔹 Final unique matches count:", nrow(final_assigned), "\n")

# Save per-match gene-level diffs for final matches
for (i in seq_len(nrow(final_assigned))) {
  r <- final_assigned$ref_cluster[i]
  n <- final_assigned$new_cluster[i]
  merged <- inner_join(
    ref_df %>% filter(cluster == r) %>% select(gene, p_val) %>% rename(p_ref = p_val),
    new_df %>% filter(cluster == n) %>% select(gene, p_val) %>% rename(p_new = p_val),
    by = "gene"
  ) %>% mutate(abs_diff = abs(p_ref - p_new), within = abs_diff <= PVALUE_DIFF_THRESH)
  write.csv(merged, file.path(OUT_DIR, paste0("per_gene_pval_ref", r, "_new", n, ".csv")), row.names = FALSE)
}

# -------------------------- summary of unmatched -----------------------------

ref_matched <- final_assigned$ref_cluster
new_matched <- final_assigned$new_cluster

ref_unmatched <- setdiff(ref_clusters, ref_matched)
new_unmatched <- setdiff(new_clusters, new_matched)

summary_list <- list(
  tibble(type = "Matched Pair",
         ref_cluster = final_assigned$ref_cluster,
         new_cluster = final_assigned$new_cluster,
         n_shared = final_assigned$n_shared,
         prop_shared = final_assigned$prop_shared,
         prop_within = final_assigned$prop_within),
  tibble(type = "Unmatched Reference",
         ref_cluster = as.character(ref_unmatched),
         new_cluster = NA_character_,
         n_shared = NA_integer_,
         prop_shared = NA_real_,
         prop_within = NA_real_),
  tibble(type = "Unmatched New",
         ref_cluster = NA_character_,
         new_cluster = as.character(new_unmatched),
         n_shared = NA_integer_,
         prop_shared = NA_real_,
         prop_within = NA_real_)
)

summary_df <- bind_rows(summary_list)

summary_csv <- file.path(OUT_DIR, "cluster_match_summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
cat("💾 Wrote cluster match summary to:", summary_csv, "\n")

cat("\n📊 Summary:\n")
cat("  ✔ Final unique matches:", nrow(final_assigned), "\n")
cat("  ❌ Unmatched reference clusters:", paste(ref_unmatched, collapse = ", "), "\n")
cat("  ❌ Unmatched new clusters:", paste(new_unmatched, collapse = ", "), "\n")
cat("✅ Completed unique p-value-based cluster matching.\n")
