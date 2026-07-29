#!/usr/bin/env Rscript

library(tidyverse)

# ──────────────────────────────────────────────
# 0. PARSE COMMAND-LINE ARGUMENTS
# ──────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: Please specify the R^2 threshold.\nUsage: Rscript filter_pgls.R <r2_thresh>\nExample: Rscript filter_pgls.R 0.70")
}

# Capture R^2 threshold from command line
r2_input  <- args[1]
r2_thresh <- as.numeric(r2_input)

if (is.na(r2_thresh)) {
  stop("Error: The provided R^2 threshold '", r2_input, "' is not a valid number.")
}

# Set P-value & FDR defaults (or optionally override via 2nd and 3rd arguments)
p_thresh   <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.05)
fdr_thresh <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.05)

# Dynamically set output directory name based on R^2 input value
output_dir <- paste0("OG_signif_", r2_input)
input_file <- "merged_pgls5_castes.tsv"

message("==================================================")
message(" Running PGLS Filter")
message(" Target R^2 Threshold : > ", r2_thresh)
message(" Target P-value       : < ", p_thresh)
message(" Target FDR           : < ", fdr_thresh)
message(" Output Directory     : ", output_dir)
message("==================================================")

# ──────────────────────────────────────────────
# 1. READ MERGED FILE
# ──────────────────────────────────────────────
if (!file.exists(input_file)) {
  stop("Input file '", input_file, "' not found! Please check file path.")
}

df_all <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Handle column naming if last column is 'motif' instead of 'motivo'
if ("motif" %in% colnames(df_all) && !"motivo" %in% colnames(df_all)) {
  df_all <- df_all %>% rename(motivo = motif)
}

# ──────────────────────────────────────────────
# 2. FILTER OUT FAIL AND NA VALUES
# ──────────────────────────────────────────────
df_clean <- df_all %>%
  filter(Status == "success") %>%
  filter(!is.na(P_value), !is.na(FDR), !is.na(R_squared_adj))

# ──────────────────────────────────────────────
# 3. FILTER BY P-VALUE & FDR (COMBINED WITH R^2)
# ──────────────────────────────────────────────
df_sig_pval <- df_clean %>%
  filter(P_value < p_thresh, R_squared_adj > r2_thresh) %>%
  arrange(motivo, P_value)

df_sig_fdr <- df_clean %>%
  filter(FDR < fdr_thresh, R_squared_adj > r2_thresh) %>%
  arrange(motivo, FDR)

message("Significant OGs by P-value (< ", p_thresh, "): ", nrow(df_sig_pval))
message("Significant OGs by FDR     (< ", fdr_thresh, "): ", nrow(df_sig_fdr))

# ──────────────────────────────────────────────
# 4. OUTPUT GENERATION
# ──────────────────────────────────────────────
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

file_pval_out <- file.path(output_dir, "pgls5_significant_OGs_pval.tsv")
file_fdr_out  <- file.path(output_dir, "pgls5_significant_OGs_fdr.tsv")

write.table(df_sig_pval, file = file_pval_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_sig_fdr,  file = file_fdr_out,  sep = "\t", quote = FALSE, row.names = FALSE)

generate_summary <- function(df, filename) {
  if (nrow(df) == 0) return()
  
  summary_df <- df %>%
    group_by(motivo) %>%
    summarise(
      n_OG_sig           = n(),
      median_R2_adj      = round(median(R_squared_adj), 4),
      median_pval        = round(median(P_value), 4),
      median_fdr         = round(median(FDR), 4),
      n_positive_coef    = sum(Regression_Coefficient > 0),
      n_negative_coef    = sum(Regression_Coefficient < 0),
      direction_dominant = ifelse(n_positive_coef >= n_negative_coef, "polymorphic", "monomorphic"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_OG_sig))

  write.table(summary_df, file = file.path(output_dir, filename), sep = "\t", quote = FALSE, row.names = FALSE)
}

generate_summary(df_sig_pval, "pgls5_summary_by_motif_pval.tsv")
generate_summary(df_sig_fdr,  "pgls5_summary_by_motif_fdr.tsv")

message("\nProcess completed successfully!")
message("→ Results saved in: ", output_dir, "/")
