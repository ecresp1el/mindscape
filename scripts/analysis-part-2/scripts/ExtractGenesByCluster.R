#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Environment variables
# ------------------------------------------------------------------------------
deg_csv   <- Sys.getenv("DEG_INPUT_CSV")
out_dir   <- Sys.getenv("DEG_OUTPUT")
genes     <- unlist(strsplit(Sys.getenv("TARGET_GENES"), ","))
identities <- unlist(strsplit(Sys.getenv("TARGET_IDENTITIES"), ","))

if (deg_csv == "" || out_dir == "" || length(genes) == 0) {
  stop("❌ DEG_INPUT_CSV, DEG_OUTPUT, or TARGET_GENES not set.")
}
if (length(genes) != length(identities)) {
  stop("❌ TARGET_GENES and TARGET_IDENTITIES must have the same length.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load DEG CSV
# ------------------------------------------------------------------------------
cat("📥 Reading DEG table from:", deg_csv, "\n")
deg_data <- read.csv(deg_csv, stringsAsFactors = FALSE)

# Ensure cluster is always character
deg_data$cluster <- as.character(deg_data$cluster)

required_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster")
missing_cols <- setdiff(required_cols, colnames(deg_data))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns in DEG file: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Helper: Strict check that filtered rows exist in original CSV
# ------------------------------------------------------------------------------
check_rows_exist_in_original <- function(filtered_df, original_df) {
  non_na_rows <- filtered_df[!apply(is.na(filtered_df), 1, all), ]
  if (nrow(non_na_rows) == 0) return(TRUE)

  # Only check overlap on required DEG columns
  common_cols <- intersect(names(original_df), names(non_na_rows))
  merged <- merge(
    non_na_rows[, common_cols, drop = FALSE],
    original_df[, common_cols, drop = FALSE],
    by = common_cols,
    all.x = TRUE
  )

  if (nrow(merged) != nrow(non_na_rows)) {
    stop("❌ Some filtered non-NA rows do NOT exist in the original DEG CSV before appending identities!")
  } else {
    cat("✅ All filtered non-NA rows exist in original CSV.\n")
  }
}

# ------------------------------------------------------------------------------
# Helper: Validate saved CSV
# ------------------------------------------------------------------------------
validate_saved_csv <- function(saved_file, df_original) {
  cat("📥 Reloading saved CSV for validation...\n")
  loaded_df <- read.csv(saved_file, stringsAsFactors = FALSE)

  # Ensure same column order
  loaded_df <- loaded_df[, names(df_original), drop = FALSE]

  # Convert everything to character for comparison
  df_compare1 <- data.frame(lapply(df_original, as.character), stringsAsFactors = FALSE)
  df_compare2 <- data.frame(lapply(loaded_df, as.character), stringsAsFactors = FALSE)

  # Replace "" with NA for fairness
  df_compare1[df_compare1 == ""] <- NA
  df_compare2[df_compare2 == ""] <- NA

  if (!identical(df_compare1, df_compare2)) {
    cat("❌ Reloaded CSV does not match the dataframe written to disk!\n")
    mismatch_file1 <- sub("\\.csv$", "_original_dump.csv", saved_file)
    mismatch_file2 <- sub("\\.csv$", "_reloaded_dump.csv", saved_file)
    write.csv(df_compare1, mismatch_file1, row.names = FALSE)
    write.csv(df_compare2, mismatch_file2, row.names = FALSE)
    stop("Saved mismatch debug dumps to: ", mismatch_file1, " and ", mismatch_file2)
  } else {
    cat("✅ Saved CSV matches the dataframe.\n")
  }
}

# ------------------------------------------------------------------------------
# Filter genes and append identities
# ------------------------------------------------------------------------------
results <- list()
missing_genes <- c()

for (i in seq_along(genes)) {
  g <- genes[i]
  identity <- identities[i]

  # Grab all rows for this gene
  df_all <- deg_data %>% filter(gene == g)

  if (nrow(df_all) == 0) {
    # Gene not present in DEG CSV at all
    cat("⚠️ Gene not found in DEG CSV:", g, "\n")
    df_all <- data.frame(
      p_val = NA_real_,
      avg_log2FC = NA_real_,
      pct.1 = NA_real_,
      pct.2 = NA_real_,
      p_val_adj = NA_real_,
      gene = g,
      cluster = NA_character_,
      status = "missing",
      stringsAsFactors = FALSE
    )
    missing_genes <- c(missing_genes, g)
  } else {
    # Split into passing and failing clusters
    df_pass <- df_all %>%
      filter(!is.na(p_val) & p_val <= 0.05) %>%
      mutate(status = "pass")
    df_fail <- df_all %>%
      filter(is.na(p_val) | p_val > 0.05) %>%
      mutate(
        p_val = NA_real_,
        avg_log2FC = NA_real_,
        pct.1 = NA_real_,
        pct.2 = NA_real_,
        p_val_adj = NA_real_,
        status = "NaP"
      )

    df_all <- bind_rows(df_pass, df_fail)

    check_rows_exist_in_original(df_pass, deg_data)
  }

  # Add identity
  df_all$cell_identity <- rep(identity, nrow(df_all))

  # Add color-coded identity column (don’t color rows with NaP/missing)
  df_all$cluster_identity_colored <- ifelse(
    df_all$status != "pass", identity,
    ifelse(df_all$avg_log2FC < 1,
           paste0("<span style='background-color:red'>", identity, "</span>"),
           ifelse(df_all$avg_log2FC > 1,
                  paste0("<span style='background-color:green'>", identity, "</span>"),
                  paste0("<span style='background-color:yellow'>", identity, "</span>")
           )
    )
  )

  results[[length(results) + 1]] <- df_all
}

final_df <- bind_rows(results) %>% arrange(cluster, gene)

# ------------------------------------------------------------------------------
# Save final CSV
# ------------------------------------------------------------------------------
out_file <- file.path(out_dir, "DEG_entries_selected_genes_res_0.2.csv")
write.csv(final_df, out_file, row.names = FALSE)

# Validate saved CSV
validate_saved_csv(out_file, final_df)

cat("✅ Saved", nrow(final_df), "rows for", length(results), "genes to", out_file, "\n")

if (length(missing_genes) > 0) {
  cat("⚠️ Genes missing entirely (NA rows inserted, status='missing'):", paste(missing_genes, collapse = ", "), "\n")
}
