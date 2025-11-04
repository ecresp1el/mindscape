#!/usr/bin/env Rscript

# ==============================================================================
# MindScape_DEG_FindExtractAttach.R
#
# PURPOSE:
# Uses a reference Seurat object to annotate a new Seurat object
# ==============================================================================

# ==============================================================================
# INPUTS / OUTPUTS
# ==============================================================================
# INPUTS:
#   - ENV VAR: DEG_INPUT ........ Path to reference Seurat RDS file for DEG analysis.
#   - ENV VAR: DEG_OUTPUT ....... Output directory for all results.
#   - ENV VAR: TARGET_GENES ..... Comma-separated list of genes of interest.
#   - ENV VAR: TARGET_IDENTITIES  Comma-separated list of corresponding identities.
#   - ENV VAR: CLUSTER_IDENTITIES (optional) Custom cluster name mapping string.
#
# OUTPUTS:
#   - DEGs_cluster_<id>.csv ............ DEG table per cluster.
#   - DEGs_all_clusters.csv ............ Combined DEG results for all clusters.
#   - DEG_entries_selected_genes_extracted.csv ... Filtered DEGs for target genes, annotated.
#   - *_dump.csv ........................ Debug CSVs written only if verification fails.
#
# NOTES:
#   • All CSVs are verified with strict integrity checks (NFS-safe).
#   • The R script can run standalone if environment variables are pre-set.
# ==============================================================================

# Load packages
suppressPackageStartupMessages({  # 👉 Suppresses startup messages for clean logs
  library(Seurat)                 # 👉 Main single-cell analysis package
  library(dplyr)                  # 👉 For data wrangling and manipulation
  library(ggplot2)                # 👉 For visualization
  library(RColorBrewer)           # 👉 Provides color palettes
  library(patchwork)              # 👉 Used for combining ggplot2 figures
  library(readr)                  # 👉 Reliable CSV reading/writing functions
})

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

# Load in environment variables
input_rds  <- Sys.getenv("DEG_INPUT")             # 👉 Input Seurat RDS file path
output_dir <- Sys.getenv("DEG_OUTPUT")            # 👉 Output directory for results
target_genes <- unlist(strsplit(Sys.getenv("TARGET_GENES"), ","))           # 👉 Genes of interest (comma-separated)
target_identities <- unlist(strsplit(Sys.getenv("TARGET_IDENTITIES"), ",")) # 👉 Identities of interest (comma-separated)

# Check for the necessary i/o being present
if (input_rds == "" || output_dir == "")
  stop("❌ DEG_INPUT or DEG_OUTPUT not set.")     # 👉 Stops execution if inputs missing
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)  # 👉 Creates output dir if it doesn’t exist

# Print information for output and targets
cat("📂 Input:", input_rds, "\n💾 Output directory:", output_dir, "\n")      # 👉 Log current file paths
if (length(target_genes) > 0) cat("🎯 Target genes:", paste(target_genes, collapse = ", "), "\n") # 👉 Log targets

# ==============================================================================
# LOAD SEURAT OBJECT
# ==============================================================================

cat("\n📖 Loading Seurat object...\n")             # 👉 Log start of object load
seu <- readRDS(input_rds)                         # 👉 Load RDS into memory
if (!inherits(seu, "Seurat")) stop("❌ Input is not a Seurat object.")  # 👉 Validate object type
DefaultAssay(seu) <- "RNA"                        # 👉 Set default assay to RNA for downstream DEGs
seu <- JoinLayers(seu)                            # 👉 Merge any split RNA layers (e.g., counts, data)
cat("✅ Seurat object loaded and RNA layers joined.\n") # 👉 Confirm success

if (!"seurat_clusters" %in% colnames(seu@meta.data))  # 👉 Check clustering column presence
  stop("❌ 'seurat_clusters' not found in metadata.")
Idents(seu) <- "seurat_clusters"                  # 👉 Use cluster assignments as active identity
clusters <- sort(unique(Idents(seu)))             # 👉 Extract sorted cluster IDs
cat("🔹 Found clusters:", paste(clusters, collapse = ", "), "\n")  # 👉 Print summary of clusters

# ==============================================================================
# HELPER: VERIFY CSV INTEGRITY (Strict / Fatal / NFS-safe)
# ==============================================================================

verify_csv_integrity <- function(df, csv_path) {
  if (!file.exists(csv_path))
    stop(paste0("❌ Verification failed — file missing: ", csv_path)) # 👉 Abort if file missing
  Sys.sleep(0.5)  # 👉 NFS-safe delay to ensure write buffer flush before reading file

  reloaded <- tryCatch(  # 👉 Safely reload the CSV file to check it matches in-memory version
    read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("❌ Could not reload CSV for verification: ", e$message)
  )

  # Define columns shared between unverified and reloaded
  common_cols <- intersect(colnames(df), colnames(reloaded)) # 👉 Find overlapping columns between two dataframes

  if (length(common_cols) == 0) {
    stop("❌ Verification failed — no matching columns in ", csv_path) # 👉 Abort if mismatch
  } else {
    message(sprintf("✅ %d shared columns found between in-memory and reloaded CSV.", length(common_cols))) # 👉 Inform success
  }

  df2 <- df[, common_cols, drop = FALSE]   # 👉 Keep only shared columns from in-memory df
  rel  <- reloaded[, common_cols, drop = FALSE] # 👉 Keep shared columns from reloaded df

  # Normalize factors/characters and safely coerce numeric-like strings
  normalize_df <- function(x) {
    x[] <- lapply(x, function(col) {                 # 👉 Iterate columns, normalize formats
      if (is.factor(col)) col <- as.character(col)   # 👉 Convert factors → characters
      if (is.character(col)) {
        col[col == ""] <- NA                         # 👉 Treat empty strings as NA
        if (all(grepl("^-?[0-9.]+$", col[!is.na(col)]))) {  # 👉 Detect if column is numeric-like
          suppressWarnings(col <- as.numeric(col))   # 👉 Convert to numeric when safe
        }
      }
      col                                            # 👉 Return cleaned column
    })
    x                                                # 👉 Return normalized dataframe
  }

  # Normalize both dataframes
  df2 <- normalize_df(df2)
  rel <- normalize_df(rel)

  # Numeric columns: allow tiny float or integer/string mismatches
  num_cols <- intersect(names(df2)[sapply(df2, is.numeric)],
                        names(rel)[sapply(rel, is.numeric)]) # 👉 Identify numeric columns shared by both
  if (length(num_cols) > 0) {
    for (col in num_cols) {
      if (!isTRUE(all.equal(df2[[col]], rel[[col]], tolerance = 1e-6))) { # 👉 Check numeric equality within tolerance
        stop(sprintf("❌ Verification failed: numeric mismatch in %s (%s)",
                     col, basename(csv_path)))                             # 👉 Stop on numeric mismatch
      }
    }
  }

  # Final equality check ignoring attributes but allowing type casting
  if (!isTRUE(all.equal(df2, rel, check.attributes = FALSE))) {  # 👉 Perform full data comparison ignoring minor attrs
    dbg1 <- sub("\\.csv$", "_original_dump.csv", csv_path)       # 👉 Dump debug files if mismatch found
    dbg2 <- sub("\\.csv$", "_reloaded_dump.csv", csv_path)
    write.csv(df2, dbg1, row.names = FALSE, quote = FALSE)
    write.csv(rel, dbg2, row.names = FALSE, quote = FALSE)
    stop(sprintf("❌ Verification failed for %s. Debug dumps written to:\n  %s\n  %s",
                 basename(csv_path), dbg1, dbg2))
  }

  message("✅ Verification passed: ", basename(csv_path))  # 👉 Confirm CSV integrity success
  TRUE                                                    # 👉 Return success flag
}

# ==============================================================================
# DEG ANALYSIS (STEP 3)
# ==============================================================================

# DEGs for reference object
cat("\n📊 Running DEG analysis for all clusters...\n")
all_markers <- list()  # 👉 Initialize empty list to store DEGs per cluster

for (cl in clusters) {  # 👉 Iterate over each cluster ID
  cat("🧪 Finding markers for cluster", cl, "...\n")
  # Find markers in an error wrapper with Wilcoxon test
  markers <- tryCatch(
    FindMarkers(seu, ident.1 = cl, only.pos = TRUE,        # 👉 Find DEGs vs all other clusters
                min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox"),
    error = function(e) { warning("⚠️ Cluster ", cl, " failed: ", e$message); return(NULL) } # 👉 Graceful fail
  )
  if (!is.null(markers) && nrow(markers) > 0) {  # 👉 Only continue if DEGs were found
    markers$gene <- rownames(markers)            # 👉 Add gene names as column
    markers$cluster <- as.character(cl)          # 👉 Record cluster ID
    out_csv <- file.path(output_dir, paste0("DEGs_cluster_", cl, ".csv"))  # 👉 Define per-cluster CSV path
    write.csv(markers, out_csv, row.names = FALSE)                          # 👉 Save DEGs to CSV
    verify_csv_integrity(markers, out_csv)                                  # 👉 Verify file correctness
    all_markers[[cl]] <- markers                                            # 👉 Store in combined list
  }
}

if (length(all_markers) == 0) stop("❌ No DEGs detected in any cluster.")   # 👉 Abort if no DEGs at all
deg_all <- bind_rows(all_markers)                                          # 👉 Combine all clusters into one table
deg_csv <- file.path(output_dir, "DEGs_all_clusters.csv")                  # 👉 Define combined DEG path
write.csv(deg_all, deg_csv, row.names = FALSE)                             # 👉 Write combined CSV
verify_csv_integrity(deg_all, deg_csv)                                     # 👉 Verify combined CSV correctness
cat("💾 Combined DEGs saved:", deg_csv, "\n")                              # 👉 Log success

# ==============================================================================
# RETRIEVE EXISTING CLUSTER IDENTITIES OR ENV VARIABLE
# ==============================================================================

cat("\n🧬 Using cluster identities from Seurat object or environment...\n")

if ("cluster_number_name" %in% colnames(seu@meta.data)) {
  cat("✅ Found existing 'cluster_number_name' in Seurat metadata.\n")  # 👉 Prefer metadata-based labels
  Idents(seu) <- "cluster_number_name"
} else {
  cluster_identity_string <- Sys.getenv("CLUSTER_IDENTITIES")           # 👉 Fallback to env variable if no metadata
  if (cluster_identity_string == "") {
    stop("❌ No 'cluster_number_name' in Seurat object and CLUSTER_IDENTITIES not provided.")  # 👉 Hard fail if missing both
  }
  cluster_pairs <- unlist(strsplit(cluster_identity_string, ","))       # 👉 Split env string into cluster=name pairs
  cluster_names <- sub("^[0-9]+=", "", cluster_pairs)                   # 👉 Extract cluster names
  cluster_ids <- sub("=.*", "", cluster_pairs)                          # 👉 Extract numeric cluster IDs
  cluster_identities <- setNames(cluster_names, cluster_ids)            # 👉 Map IDs → names
  cat("✅ Parsed CLUSTER_IDENTITIES from environment:\n")
  print(cluster_identities)

  Idents(seu) <- "cluster_number_name"                                  # 👉 (Placeholder) maintain consistent metadata ref
}

identities <- levels(seu)                                                # 👉 List all current active identities
cat("🔹 Active identities:\n")                                           # 👉 Log identities to console
print(identities)

# ==============================================================================
# AUTO-EXTRACT DEG ENTRIES FROM INTERNAL CSV (STEP 5b)
# ==============================================================================

cat("\n📄 Extracting desired marker entries from internally generated DEG CSV...\n")
if (length(target_genes) == 0) {  # 👉 Begin filtering block if target genes are specified
  cat("⚠️ No TARGET_GENES specified. Skipping extraction step.\n")
} else {

  # Helper functions
  check_rows_exist_in_original <- function(filtered_df, original_df) {
    # 💬 Remove rows that are entirely NA to ensure we only check meaningful entries
    non_na_rows <- filtered_df[!apply(is.na(filtered_df), 1, all), ] # Store value-containing rows
    if (nrow(non_na_rows) == 0) return(TRUE) # 💬 Nothing to check — safe to exit early
    common_cols <- intersect(names(original_df), names(non_na_rows))  # 💬 Columns shared between filtered and original
    
    # 💬 Merge ensures each filtered row matches an existing original one
    merged <- merge( 
      non_na_rows[, common_cols, drop = FALSE],
      original_df[, common_cols, drop = FALSE],
      by = common_cols,
      all.x = TRUE
    )
    
    # 💬 If merged rows < filtered rows, then some filtered entries were not found
    if (nrow(merged) != nrow(non_na_rows)) {
      stop("❌ Some filtered non-NA rows do NOT exist in the original DEG CSV before appending identities!")
    } else {
      cat("✅ All filtered non-NA rows exist in original CSV.\n")
    }
  }

  validate_saved_csv <- function(saved_file, df_original) {
    cat("📥 Reloading saved CSV for validation...\n")
    loaded_df <- read.csv(saved_file, stringsAsFactors = FALSE) # 💬 Read file freshly from disk to validate
    loaded_df <- loaded_df[, names(df_original), drop = FALSE]  # 💬 Ensure column order matches original
    
    # 💬 Convert both dataframes to all-character form for direct structural comparison
    df_compare1 <- data.frame(lapply(df_original, as.character), stringsAsFactors = FALSE)
    df_compare2 <- data.frame(lapply(loaded_df, as.character), stringsAsFactors = FALSE)
    df_compare1[df_compare1 == ""] <- NA # 💬 Empty string normalization
    df_compare2[df_compare2 == ""] <- NA
    
    # 💬 Full equality test ignoring type differences (since coerced to character)
    if (!identical(df_compare1, df_compare2)) {
      cat("❌ Reloaded CSV does not match the dataframe written to disk!\n")
      # 💬 Write mismatch debug dumps for inspection
      mismatch_file1 <- sub("\\.csv$", "_original_dump.csv", saved_file)
      mismatch_file2 <- sub("\\.csv$", "_reloaded_dump.csv", saved_file)
      write.csv(df_compare1, mismatch_file1, row.names = FALSE)
      write.csv(df_compare2, mismatch_file2, row.names = FALSE)
      stop("Saved mismatch debug dumps to: ", mismatch_file1, " and ", mismatch_file2)
    } else {
      cat("✅ Saved CSV matches the dataframe.\n")
    }
  }

  # Load internal DEGs CSV
  cat("📥 Reading internally generated DEGs from:", deg_csv, "\n")
  deg_data <- read.csv(deg_csv, stringsAsFactors = FALSE)
  deg_data$cluster <- as.character(deg_data$cluster) # 💬 Ensure clusters are treated as string IDs

  # 💬 Verify required columns exist in the DEG dataset
  required_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster")
  missing_cols <- setdiff(required_cols, colnames(deg_data))
  if (length(missing_cols) > 0) { 
    stop("❌ Missing required columns in DEG file: ", paste(missing_cols, collapse = ", "))
  }

  # 💬 Ensure gene and identity vectors align
  if (length(target_genes) != length(target_identities)) {
    stop("❌ TARGET_GENES and TARGET_IDENTITIES must have the same length.")
  }

  results <- list()          # 💬 Will hold filtered gene-specific DEG tables
  missing_genes <- c()       # 💬 Track genes not found in the CSV

  # Filter for all target genes
  for (i in seq_along(target_genes)) {
    g <- target_genes[i]          # 💬 Current gene symbol
    identity <- target_identities[i]  # 💬 Matched identity annotation
    cat("\n🔍 Processing gene:", g, "for identity:", identity, "\n")

    # 💬 Filter DEG table for the current gene
    df_all <- deg_data %>% filter(gene == g)

    if (nrow(df_all) == 0) { # 💬 If not found, record as missing placeholder
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
      # 💬 Split gene entries into significant vs. non-significant (p ≤ 0.05)
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
          status = "NaP" # 💬 NaP = Not passing threshold
        )
      
      # 💬 Combine both for full traceability
      df_all <- bind_rows(df_pass, df_fail)
      # 💬 Confirm all passing rows exist in the original dataset
      check_rows_exist_in_original(df_pass, deg_data)
    }

    # 💬 Add identity metadata and HTML-colored tags for up/downregulation
    df_all$cell_identity <- rep(identity, nrow(df_all))
    df_all$cluster_identity_colored <- ifelse(
      df_all$status != "pass", identity,  # 💬 Non-significant rows just keep identity
      ifelse(df_all$avg_log2FC < 1,
             paste0("<span style='background-color:red'>", identity, "</span>"),  # 💬 Downregulated
             ifelse(df_all$avg_log2FC > 1,
                    paste0("<span style='background-color:green'>", identity, "</span>"), # 💬 Upregulated
                    paste0("<span style='background-color:yellow'>", identity, "</span>") # 💬 Neutral
             )
      )
    )
    results[[length(results) + 1]] <- df_all  # 💬 Append to results list
  }

  # 💬 Merge all per-gene tables together
  final_df <- bind_rows(results) %>% arrange(cluster, gene)

  # 💬 Write combined extracted DEGs for target genes
  out_file <- file.path(output_dir, "DEG_entries_selected_genes_extracted.csv")
  write.csv(final_df, out_file, row.names = FALSE)
  validate_saved_csv(out_file, final_df)  # 💬 Integrity check for saved CSV

  cat("✅ Saved", nrow(final_df), "rows for", length(results), "genes to", out_file, "\n")
  if (length(missing_genes) > 0) {
    cat("⚠️ Genes missing entirely (NA rows inserted, status='missing'):", 
        paste(missing_genes, collapse = ", "), "\n")
  }
}

cat("\n🎉 Complete: DEG + marker extraction using internally generated DEGs and existing identities.\n")

# ==============================================================================
# STEP 6 — IDENTIFY TOP UNIQUE UPREGULATED GENES PER CLUSTER (FLEXIBLE, REFERENCE)
# ==============================================================================

cat("\n🌟 Identifying top unique upregulated genes per cluster (flexible threshold) — REFERENCE...\n")

final_extracted_path <- file.path(output_dir, "DEG_entries_selected_genes_extracted.csv")
if (!file.exists(final_extracted_path)) {
  stop("❌ File not found: ", final_extracted_path, 
       "\nCannot compute unique genes without extracted DEG entries.")
}

# Read the file into deg_extracted
deg_extracted <- read.csv(final_extracted_path, stringsAsFactors = FALSE)

# Confirmation message
if (nrow(deg_extracted) == 0) {
  warning("⚠️ File exists but contains no rows!")
} else {
  cat("🌟 SUCCESS: DEG entries loaded with", nrow(deg_extracted), "rows and", ncol(deg_extracted), "columns.\n")
}

# 💬 Verify necessary columns exist before proceeding
required_cols <- c("gene", "cluster", "avg_log2FC", "p_val", "status")
missing_cols <- setdiff(required_cols, colnames(deg_extracted))
if (length(missing_cols) > 0) {
  stop("❌ Missing required columns in DEG_entries_selected_genes_extracted.csv: ", 
       paste(missing_cols, collapse = ", "))
}

# 💬 Keep only valid, upregulated, significant genes
deg_filtered <- deg_extracted %>%
  filter(status == "pass" & !is.na(avg_log2FC) & avg_log2FC > 0 &
         !is.na(gene) & !is.na(cluster))

if (nrow(deg_filtered) == 0) {
  stop("❌ No upregulated significant genes found for unique-gene analysis.")
}

# -------------------- FLEXIBLE UNIQUENESS PARAMETER (shared) -------------------
# 💬 Parameters are read from environment for reproducibility across datasets
max_clusters_allowed <- as.numeric(Sys.getenv("MAX_CLUSTER_UNIQUENESS", "1"))
n_top_per_cluster  <- as.integer(Sys.getenv("TOP_N_PER_CLUSTER", "5"))

# 💬 Sanity checks for environment parameters
if (length(max_clusters_allowed) != 1 || is.na(max_clusters_allowed)) {
  stop("❌ Environment variable MAX_CLUSTER_UNIQUENESS must be a single numeric value (e.g. 1, 3, 6, 12).")
}
if (length(n_top_per_cluster) != 1 || is.na(n_top_per_cluster) || n_top_per_cluster < 1) {
  stop("❌ Environment variable TOP_N_PER_CLUSTER must be a single positive integer (e.g. 5).")
}

cat("ℹ️ Selecting genes appearing in at most", max_clusters_allowed, "cluster(s)...\n")
cat("ℹ️ Selecting top", n_top_per_cluster, "genes per cluster for summary.\n")
cat("ℹ️ Using MAX_CLUSTER_UNIQUENESS =", max_clusters_allowed, "for REFERENCE dataset.\n")

# 💬 Count in how many clusters each gene is significant (to measure specificity)
gene_cluster_counts <- deg_filtered %>%
  distinct(gene, cluster) %>%
  count(gene, name = "cluster_count") %>%
  arrange(gene)

cat("✅ gene_cluster_counts (unique clusters per gene):\n")
print(gene_cluster_counts)
cat("\n")

# 💬 Keep genes expressed significantly in ≤ N clusters
selected_genes <- gene_cluster_counts %>%
  filter(cluster_count <= max_clusters_allowed) %>%
  pull(gene) %>%
  unique()

cat("✅ selected_genes (genes appearing in <= ", max_clusters_allowed, " clusters):\n", sep = "")
print(selected_genes)
cat("\n")

deg_selected <- deg_filtered %>% filter(gene %in% selected_genes)

if (nrow(deg_selected) == 0) {
  cat("⚠️ No genes found with the specified uniqueness threshold for reference.\n")
} else {
  cat("✅ Found", length(unique(deg_selected$gene)), 
      "genes passing the uniqueness threshold across clusters (reference).\n")

  # 🧠 For each cluster, take top N genes (highest avg_log2FC) among those passing uniqueness.
  # Deterministic ordering: ensures consistent ordering across runs.
  top_per_cluster <- deg_selected %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), gene, .by_group = TRUE) %>% # sort by effect size (descending)
    slice_head(n = n_top_per_cluster) %>%  # keep top N per cluster
    ungroup()

  head(top_per_cluster, 20) # 🧠 preview first 20 entries for debugging

  # 🧠 Create individual CSV files for each cluster’s top genes
  clusters_unique <- sort(unique(top_per_cluster$cluster))
  for (cl in clusters_unique) {
    df_cl <- top_per_cluster %>% filter(cluster == cl)
    out_path <- file.path(output_dir, paste0("TopGenes_cluster_", cl, ".csv"))
    write.csv(df_cl, out_path, row.names = FALSE)
    cat("💾 Saved top genes for cluster", cl, "->", out_path, "\n")
  }

  # 🧠 Save combined deterministic summary (used later by matching step)
  combined_path <- file.path(output_dir, "TopUniqueGenes_all_clusters.csv")
  # ensure deterministic column order
  write.csv(top_per_cluster[order(as.numeric(top_per_cluster$cluster), top_per_cluster$gene),
                            , drop = FALSE], combined_path, row.names = FALSE)
  cat("✅ Combined reference unique summary saved:", combined_path, "\n")
}

# ==============================================================================
# STEP 7 — MATCH NEW DATASET CLUSTERS TO REFERENCE IDENTITIES (WITH CSV OUTPUT)
# ==============================================================================

cat("\n🧬 Matching new dataset clusters to reference identities using unique genes...\n")

# 🧠 Load and prepare the "new" Seurat dataset to be annotated
new_input_rds <- Sys.getenv("NEW_INPUT_RDS")
if (new_input_rds == "") {
  cat("⚠️ NEW_INPUT_RDS not provided. Skipping cluster matching step.\n")
} else if (!file.exists(new_input_rds)) {
  stop("❌ NEW_INPUT_RDS not found:", new_input_rds)
} else {
  cat("📥 Loading new Seurat object:", new_input_rds, "\n")
  new_seu <- readRDS(new_input_rds)
  if (!inherits(new_seu, "Seurat")) stop("❌ NEW_INPUT_RDS is not a Seurat object.")
  DefaultAssay(new_seu) <- "RNA"
  new_seu <- JoinLayers(new_seu)  # 🧠 join RNA layers for analysis
  Idents(new_seu) <- "seurat_clusters"  # 🧠 set cluster identities
  new_clusters <- sort(unique(Idents(new_seu)))
  cat("🔹 Found clusters in new dataset:", paste(new_clusters, collapse = ", "), "\n")
}

# -------------------- Run DEG analysis for new dataset ---------------------
  cat("\n📊 Running DEG analysis for new dataset...\n")
  new_markers <- list()
  for (cl in new_clusters) {
    cat("🧪 Finding markers for new cluster", cl, "...\n")
    mk <- tryCatch(
      FindMarkers(new_seu, ident.1 = cl, only.pos = TRUE,
                  min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox"),
      error = function(e) { warning("⚠️ New cluster ", cl, " failed: ", e$message); return(NULL) }
    )
    if (!is.null(mk) && nrow(mk) > 0) {
      mk$gene <- rownames(mk)  # 🧠 record gene names
      mk$cluster <- as.character(cl)
      new_markers[[cl]] <- mk
    }
  }
  new_deg_all <- bind_rows(new_markers)
  cat("✅ DEG analysis complete for new dataset (", nrow(new_deg_all), "rows)\n")

  # -------------------- Filter for TARGET_GENES if provided -----------------
if (length(target_genes) == 0) {
  cat("⚠️ No TARGET_GENES specified. Using all genes for further analysis.\n")
  deg_new_filtered <- new_deg_all %>%
    filter(!is.na(p_val) & p_val <= 0.05 & avg_log2FC > 0) %>% # 🧠 keep sig. upregulated
    mutate(status = "pass")
} else {
  deg_new_filtered <- list()
  missing_genes_new <- c()

  # 🧠 Mirror the filtering logic used for the reference dataset
  for (i in seq_along(target_genes)) {
    g <- target_genes[i]
    identity <- target_identities[i]
    cat("\n🔍 Processing target gene:", g, "for identity:", identity, "\n")
    
    df_gene <- new_deg_all %>% filter(gene == g)
    
    if (nrow(df_gene) == 0) {
      # 🧠 Gene missing entirely from DEGs → mark with NA row
      cat("⚠️ Gene not found in new dataset:", g, "\n")
      df_gene <- data.frame(
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
      missing_genes_new <- c(missing_genes_new, g)
    } else {
      # 🧠 Split into passing (p ≤ 0.05) and failing (non-significant) rows
      df_pass <- df_gene %>% 
        filter(!is.na(p_val) & p_val <= 0.05) %>%
        mutate(status = "pass")
      
      df_fail <- df_gene %>%
        filter(is.na(p_val) | p_val > 0.05) %>%
        mutate(
          p_val = NA_real_,
          avg_log2FC = NA_real_,
          pct.1 = NA_real_,
          pct.2 = NA_real_,
          p_val_adj = NA_real_,
          status = "NaP"
        )
      
      df_gene <- bind_rows(df_pass, df_fail)
      # check_rows_exist_in_original(df_pass, new_deg_all) # 🧠 optional consistency check
    }

    deg_new_filtered[[length(deg_new_filtered) + 1]] <- df_gene
  }

  deg_new_filtered <- bind_rows(deg_new_filtered)
  
  if (length(missing_genes_new) > 0) {
    cat("⚠️ Target genes missing entirely in new dataset:", paste(missing_genes_new, collapse = ", "), "\n")
  }
}

  # -------------------- Identify top unique genes per cluster (NEW dataset; uses identical logic) -----------------
  cat("\n🌟 Identifying top unique upregulated genes per cluster (flexible threshold) — NEW dataset...\n")

  # 🧠 Safety check: skip uniqueness analysis if no data available
  if (!exists("deg_new_filtered") || nrow(deg_new_filtered) == 0) {
    cat("⚠️ No DEG rows for new dataset to analyze for uniqueness.\n")
    top_new_unique <- tibble::tibble()
  } else {
    # 🧠 Normalize and filter using same significance and direction thresholds
    deg_new_filtered <- deg_new_filtered %>%
      mutate(cluster = as.character(cluster))

    deg_new_pass_filtered <- deg_new_filtered %>%
      filter(status == "pass" & !is.na(avg_log2FC) & avg_log2FC > 0 &
             !is.na(gene) & !is.na(cluster))

    # 🧠 Count number of clusters each gene appears in (for uniqueness)
    new_gene_cluster_counts <- deg_new_pass_filtered %>%
      distinct(gene, cluster) %>%
      count(gene, name = "cluster_count") %>%
      arrange(gene)
  }
  
# ---------------------------------------------------------------------------
# 🧠 Print summary of how many clusters each gene appears in
# ---------------------------------------------------------------------------
cat("✅ new_gene_cluster_counts (unique clusters per gene):\n")
print(new_gene_cluster_counts)
cat("\n")

# 🧠 Select genes that appear in ≤ threshold number of clusters
# Here, genes are considered "unique" if they appear in no more than
# `max_clusters_allowed` clusters (helps identify cluster-specific markers)
new_unique_genes <- new_gene_cluster_counts %>%
  filter(cluster_count <= max_clusters_allowed) %>%  # keep genes meeting the uniqueness constraint
  pull(gene) %>%                                    # extract just the gene column
  unique()                                          # ensure distinct gene names

# Print a list of genes passing the uniqueness threshold
cat("✅ new_unique_genes (genes appearing in <= ", max_clusters_allowed, " clusters):\n", sep = "")
print(new_unique_genes)
cat("\n")

# 🧠 Filter the DEG table to include only these cluster-specific genes
# This restricts the differential expression results to cluster-unique genes
deg_new_unique <- deg_new_pass_filtered %>%
  filter(gene %in% new_unique_genes)

# Handle the case where no genes pass the uniqueness filter
if (nrow(deg_new_unique) == 0) {
  cat("⚠️ No genes found in new dataset with the specified uniqueness threshold.\n\n")
  top_new_unique <- tibble::tibble()  # create empty tibble for compatibility downstream
} else {
  # Report how many unique genes remain after filtering
  cat("✅ Found", length(unique(deg_new_unique$gene)), 
      "genes passing the uniqueness threshold across new dataset clusters.\n\n")

  # 🧠 Recreate deterministic per-cluster ranking (same as reference)
  # Rank genes within each cluster by descending log fold-change, breaking ties by gene name
  top_new_unique <- deg_new_unique %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), gene, .by_group = TRUE) %>%  # ensure consistent order
    slice_head(n = n_top_per_cluster) %>%                  # keep only top N per cluster
    ungroup()

  # Display the top unique genes table
  cat("✅ top_new_unique (top", n_top_per_cluster, "genes per cluster):\n")
  print(top_new_unique)
  cat("\n")

  # 🧠 Compact textual summary for quick inspection
  # Generate an easy-to-read string of top genes per cluster
  cat("📌 Top genes per cluster (compact view, new dataset):\n")
  top_new_unique %>%
    group_by(cluster) %>%
    summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop") %>%
    print()

  # Save new-dataset combined summary (distinct filename)
  # This file contains all top unique genes across clusters in the new dataset
  combined_new_path <- file.path(output_dir, "TopUniqueGenes_new_all_clusters.csv")
  write.csv(top_new_unique[order(as.numeric(top_new_unique$cluster), top_new_unique$gene),
                           , drop = FALSE],
            combined_new_path, row.names = FALSE)
  cat("✅ Combined new-dataset unique summary saved:", combined_new_path, "\n")
}
  # Close of conditional block

# ---------------------------------------------------------------------------
# -------------------- Compare with reference -------------------------------
# ---------------------------------------------------------------------------

# Construct expected reference file path and confirm its existence
ref_path <- file.path(output_dir, "TopUniqueGenes_all_clusters.csv")
file.exists(ref_path)

# Stop if reference file cannot be found — it's required for comparison
if (!file.exists(ref_path)) stop("❌ Reference unique-gene file not found:", ref_path)
ref_unique <- read.csv(ref_path, stringsAsFactors = FALSE)

# Validate that required columns exist in the reference dataset
if (!"cluster" %in% names(ref_unique))
  stop("❌ 'cluster' column missing in reference unique-gene file.")
if (!"gene" %in% names(ref_unique))
  stop("❌ 'gene' column missing in reference unique-gene file.")

# Split reference and new gene lists by cluster for overlap comparison
ref_groups <- split(ref_unique$gene, ref_unique$cluster)
new_groups <- split(top_new_unique$gene, top_new_unique$cluster)

# Prepare an empty list to store overlap statistics between cluster pairs
match_summary <- list()

# Iterate over each new cluster and compare to every reference cluster
for (new_cl in names(new_groups)) {
  new_genes <- new_groups[[new_cl]]
  row_entries <- list()
  for (ref_cl in names(ref_groups)) {
    ref_genes <- ref_groups[[ref_cl]]
    overlap_genes <- intersect(ref_genes, new_genes)      # shared genes between clusters
    overlap_n <- length(overlap_genes)                    # count overlap
    row_entries[[ref_cl]] <- list(
      new_cluster = new_cl,
      ref_cluster = ref_cl,
      overlap_count = overlap_n,
      overlap_genes = paste(overlap_genes, collapse = ";") # store overlap genes as string
    )
  }
  match_summary[[new_cl]] <- bind_rows(row_entries)        # combine rows per new cluster
}

# ---------------------------------------------------------------------------
# -------------------- Build and Save Enhanced Cluster Match Summary --------
# ---------------------------------------------------------------------------

# Combine all per-cluster comparison tables into one full data frame
match_df_full <- bind_rows(match_summary)

# Retrieve reference cluster identity names if available in metadata
if ("cluster_number_name" %in% colnames(seu@meta.data)) {
  # Use the mapping from the original Seurat object to attach descriptive names
  ref_map <- seu@meta.data %>%
    distinct(seurat_clusters, cluster_number_name) %>%
    rename(ref_cluster = seurat_clusters, ref_cluster_name = cluster_number_name)
} else {
  # Fallback: build a minimal map with numeric cluster IDs only
  ref_map <- data.frame(
    ref_cluster = unique(ref_unique$cluster),
    ref_cluster_name = as.character(unique(ref_unique$cluster)),
    stringsAsFactors = FALSE
  )
}

# Sanity check that the reference map was created successfully
if (nrow(ref_map) == 0) {
  stop("❌ ref_map is empty — could not extract cluster names from reference Seurat object.")
}

# Add readable cluster names to the match summary via join
match_df_full_named <- match_df_full %>%
  left_join(ref_map, by = "ref_cluster")

# Fallback if the name column is missing after the join
if (!"ref_cluster_name" %in% colnames(match_df_full_named)) {
  warning("⚠️ 'ref_cluster_name' column missing after join — creating placeholder.")
  match_df_full_named$ref_cluster_name <- as.character(match_df_full_named$ref_cluster)
}

# Read overlap threshold from environment variable (default = 3)
overlap_gene_threshold <- as.numeric(Sys.getenv("OVERLAP_GENE_THRESHOLD", "3"))
if (is.na(overlap_gene_threshold) || overlap_gene_threshold < 1) {
  overlap_gene_threshold <- 3
  cat("⚠️ Invalid or missing OVERLAP_GENE_THRESHOLD; defaulting to 3.\n")
} else {
  cat("🔧 Using overlap gene threshold:", overlap_gene_threshold, "\n")
}

# Label formatting and column cleanup
# This ensures clusters are labeled as "id - name" when needed
match_df_full_named <- match_df_full_named %>%
  rowwise() %>%
  mutate(
    ref_cluster_label = if (grepl(paste0("^", ref_cluster), ref_cluster_name)) {
      ref_cluster_name
    } else {
      paste0(ref_cluster, " - ", ref_cluster_name)
    },
    overlap_gene_count = overlap_count
  ) %>%
  ungroup() %>%
  select(
    new_cluster,
    ref_cluster,
    ref_cluster_label,
    overlap_gene_count,
    overlap_genes
  ) %>%
  arrange(as.numeric(new_cluster), desc(overlap_gene_count))

# Apply threshold filtering to keep only cluster pairs with enough shared genes
match_df_filtered <- match_df_full_named %>%
  dplyr::filter(overlap_gene_count >= overlap_gene_threshold)

# Warn if no matched pairs pass the threshold
if (nrow(match_df_filtered) == 0) {
  warning(paste0("⚠️ No matched cluster pairs passed the ≥ ", 
                 overlap_gene_threshold, " shared genes threshold."))
  cat("ℹ️ Proceeding with empty mapping — no clusters will be assigned verified_cluster_identity.\n")
} else {
  cat("✅ Matched cluster pairs passing threshold:\n")
  print(match_df_filtered)
}

# Save the comprehensive cluster match summary for reference
csv_out <- file.path(output_dir, "Cluster_Match_Summary.csv")
write.csv(match_df_full_named, csv_out, row.names = FALSE)
cat("💾 Enhanced cluster match summary saved:", csv_out, "\n")

# ------------------------------------------------------------------------------
# ✅ Append verified cluster identity metadata (safe, deterministic mapping)
# ------------------------------------------------------------------------------

cat("🧬 Attaching verified cluster identities to Seurat object (safe mapping)...\n")

# Check that the new Seurat object has the expected cluster metadata
if (!"seurat_clusters" %in% colnames(new_seu@meta.data)) {
  stop("❌ seurat_clusters metadata missing — cannot map cluster identities.")
}

# Ensure match_df_filtered is non-empty before proceeding
if (!exists("match_df_filtered") || nrow(match_df_filtered) == 0) {
  warning("⚠️ match_df_filtered is empty or missing — no mappings to apply. Skipping mapping.")
  # create NA column so downstream scripts don't break
  new_seu$verified_cluster_identity <- NA_character_
} else {
  # Convert columns to safe types for joining
  match_df_filtered <- match_df_filtered %>%
    mutate(
      new_cluster = as.character(new_cluster),
      ref_cluster = as.character(ref_cluster),
      ref_cluster_label = as.character(ref_cluster_label),
      overlap_gene_count = as.integer(overlap_gene_count)
    )

  # Choose best-matching reference cluster per new cluster (highest overlap)
  best_matches <- match_df_filtered %>%
    group_by(new_cluster) %>%
    arrange(new_cluster, desc(overlap_gene_count), ref_cluster_label) %>%
    slice_head(n = 1) %>%
    ungroup()

  # Build a named mapping: names(new_cluster) → ref_cluster_label
  mapping <- setNames(best_matches$ref_cluster_label, best_matches$new_cluster)

  # Apply mapping across all cells by their cluster assignment
  current_clusters <- as.character(new_seu$seurat_clusters)
  unique_current_clusters <- sort(unique(current_clusters))

  # Print mapping overview for transparency
  clusters_with_mapping <- intersect(unique_current_clusters, names(mapping))
  clusters_without_mapping <- setdiff(unique_current_clusters, names(mapping))
  cat("ℹ️ Clusters in new_seu:", paste(unique_current_clusters, collapse = ", "), "\n")
  cat("ℹ️ Will map the following clusters:", paste(clusters_with_mapping, collapse = ", "), "\n")
  if (length(clusters_without_mapping) > 0) {
    warning("⚠️ The following clusters in new_seu have no mapping in match_df_filtered and will be NA: ",
            paste(clusters_without_mapping, collapse = ", "))
  }

  # Perform actual mapping to create verified_cluster_identity
  new_verified <- mapping[current_clusters]   # NA assigned automatically for missing keys

  # Validate mapping length consistency
  if (length(new_verified) != length(current_clusters)) {
    stop("❌ Unexpected error: length mismatch when creating verified_cluster_identity.")
  }

  # Attach the new verified identity column to metadata
  new_seu$verified_cluster_identity <- as.character(new_verified)

  # Report number of mapped vs total cells
  n_mapped_cells <- sum(!is.na(new_seu$verified_cluster_identity))
  n_total_cells  <- length(new_seu$verified_cluster_identity)
  cat(sprintf("✅ Applied mapping: %d / %d cells assigned verified_cluster_identity.\n",
              n_mapped_cells, n_total_cells))

  # Warn if mapping exists but produced no assigned cells
  if (n_mapped_cells == 0) {
    warning("⚠️ Mapping applied but no cells received a verified identity (check match_df_filtered/new cluster labels).")
  }
}

# ------------------------------------------------------------------------------
# ✅ Save updated Seurat object and print metadata for verification
# ------------------------------------------------------------------------------

# Save the Seurat object that now includes the verified cluster identity metadata
out_rds_verified <- file.path(output_dir, "clustered_with_verified_identities.rds")
saveRDS(new_seu, out_rds_verified)
cat(paste0("💾 Verified Seurat object saved to: ", out_rds_verified, "\n\n"))

# Display a compact preview of metadata to verify mapping results
cat("📊 Metadata verification (showing new columns):\n")
meta_preview <- new_seu@meta.data %>%
  dplyr::select(seurat_clusters, verified_cluster_identity) %>%
  dplyr::distinct() %>%
  dplyr::arrange(as.numeric(as.character(seurat_clusters)))

print(meta_preview, row.names = FALSE)
cat("\n✅ Metadata successfully appended and verified.\n")
