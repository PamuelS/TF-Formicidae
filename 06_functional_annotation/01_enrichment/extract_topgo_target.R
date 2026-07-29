#!/usr/bin/env Rscript

library(tidyverse)

# ──────────────────────────────────────────────
# 1. SETTINGS & PATHS
# ──────────────────────────────────────────────
base_input_dir <- "/DATASMALL/samuel.pederzini/TF-Formicidae/05_pgls_analysis/01_pgls_results"

# Map input result folders to output analysis folders
folder_map <- list(
  list(input_dir = file.path(base_input_dir, "OG_signif_0.25"), output_dir = "01_enriched_0.25"),
  list(input_dir = file.path(base_input_dir, "OG_signif_0.40"), output_dir = "02_enriched_0.40"),
  list(input_dir = file.path(base_input_dir, "OG_signif_0.50"), output_dir = "03_enriched_0.50"),
  list(input_dir = file.path(base_input_dir, "OG_signif_0.70"), output_dir = "04_enriched_0.70")
)

# ──────────────────────────────────────────────
# 2. HELPER FUNCTION: DIRECT SPLIT BY SLOPE SIGN
# ──────────────────────────────────────────────
extract_and_save_targets <- function(sig_data, output_subfolder, analysis_label) {
  
  if (!dir.exists(output_subfolder)) {
    dir.create(output_subfolder, recursive = TRUE)
  }

  # Identify OG/Gene column automatically
  og_col <- case_when(
    "OG" %in% colnames(sig_data) ~ "OG",
    "Orthogroup" %in% colnames(sig_data) ~ "Orthogroup",
    "Gene" %in% colnames(sig_data) ~ "Gene",
    TRUE ~ colnames(sig_data)[1]
  )

  # 1. Extract ALL unique OGs with a POSITIVE slope (> 0)
  polymorphic_OGs <- sig_data %>%
    filter(Regression_Coefficient > 0) %>%
    pull(all_of(og_col)) %>% 
    unique()

  # 2. Extract ALL unique OGs with a NEGATIVE slope (< 0)
  monomorphic_OGs <- sig_data %>%
    filter(Regression_Coefficient < 0) %>%
    pull(all_of(og_col)) %>% 
    unique()

  poly_path <- file.path(output_subfolder, "topGO_target_Polymorphic_OGs.txt")
  mono_path <- file.path(output_subfolder, "topGO_target_Monomorphic_OGs.txt")

  # Write plain text target files
  write.table(polymorphic_OGs, poly_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(monomorphic_OGs, mono_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

  message("  [", analysis_label, "] Saved -> ", output_subfolder, "/")
  message("      - Total Rows Processed : ", nrow(sig_data))
  message("      - Positive OGs (>0)    : ", length(polymorphic_OGs))
  message("      - Negative OGs (<0)    : ", length(monomorphic_OGs))
}

# ──────────────────────────────────────────────
# 3. MAIN LOOP ACROSS ALL FOLDERS
# ──────────────────────────────────────────────
for (item in folder_map) {
  in_dir  <- item$input_dir
  out_dir <- item$output_dir

  message("\n==================================================")
  message(" Reading from: ", in_dir)
  message(" Outputting to: ", out_dir)
  message("==================================================")

  if (!dir.exists(in_dir)) {
    message(" WARNING: Input directory does not exist: ", in_dir)
    next
  }

  # 1. Process P-value Analysis -> 00_pvalue_analysis
  pval_file <- file.path(in_dir, "pgls5_significant_OGs_pval.tsv")
  if (file.exists(pval_file)) {
    df_pval <- read.delim(pval_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    out_pval_dir <- file.path(out_dir, "00_pvalue_analysis")
    extract_and_save_targets(df_pval, out_pval_dir, "P-VALUE")
  } else {
    message("  [P-VALUE] File not found: ", pval_file)
  }

  # 2. Process FDR Analysis -> 01_fdr_analysis
  fdr_file <- file.path(in_dir, "pgls5_significant_OGs_fdr.tsv")
  if (file.exists(fdr_file)) {
    df_fdr <- read.delim(fdr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    out_fdr_dir <- file.path(out_dir, "01_fdr_analysis")
    extract_and_save_targets(df_fdr, out_fdr_dir, "FDR")
  } else {
    message("  [FDR] File not found: ", fdr_file)
  }
}

message("\nDirect OG extraction complete!")
