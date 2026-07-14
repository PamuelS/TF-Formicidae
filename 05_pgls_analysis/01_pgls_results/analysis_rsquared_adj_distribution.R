library(ggplot2)

# 1. Set your path
main_dir <- "."

# 2. Get the folders
motif_dirs <- list.dirs(path = main_dir, full.names = TRUE, recursive = FALSE)
exclude_folders <- c("00_OG_significant", "01_OG_FDR_significant", "tsv_version")

# Remove excluded folders from the list
motif_dirs <- motif_dirs[!basename(motif_dirs) %in% exclude_folders]

cat("========================================\n")
cat("Starting Analysis on", length(motif_dirs), "folders...\n")
cat("========================================\n\n")

# 3. Loop through each folder
for (motif_path in motif_dirs) {

  motif_name <- basename(motif_path)

  # PRINT TO TERMINAL: Let the user know R entered the folder
  cat(" Entering folder [", motif_name, "]\n", sep="")

  # Explicitly target only 'pgls5_castes.csv'
  file_path <- file.path(motif_path, "pgls5_castes.csv")

  # Loudly warn if the specific file is missing
  if (!file.exists(file_path)) {
    cat("    SKIPPED: 'pgls5_castes.csv' not found inside this folder.\n\n")
    next
  }

  cat("    Found file:", basename(file_path), "\n")

  # Read the file with its header
  current_data <- read.csv(file_path, header = TRUE)

  # Force the 7th column to numeric, handling any non-numeric values
  current_data[, 7] <- as.numeric(as.character(current_data[, 7]))
  plot_data <- data.frame(r_squared_adj = current_data[, 7])
  plot_data <- na.omit(plot_data)

  # =========================================================================
  # PART 1: CREATE & SAVE THE GGPLOT PNG
  # =========================================================================
  p <- ggplot(plot_data, aes(x = r_squared_adj)) +
    geom_histogram(binwidth = 0.25, fill = "steelblue", color = "white", boundary = 0) +
    geom_vline(xintercept = 0.25, linetype = "dashed", color = "red", size = 0.8) +
    labs(title = paste("R² Adj Distribution:", motif_name), x = "R² Adj", y = "Frequency") +
    theme_minimal() + 
    xlim(-1, 1)

  plot_save_path <- file.path(motif_path, paste0(motif_name, "_r_squared_adj_distribution.png"))
  ggsave(filename = plot_save_path, plot = p, width = 6, height = 4, dpi = 300)

  # =========================================================================
  # PART 2: FILTER & SAVE SIGNIFICANT RESULTS (> 0.25)
  # =========================================================================
  significant_data <- current_data[which(current_data[, 7] > 0.25), ]
  results_save_path <- file.path(motif_path, paste0(motif_name, "_r_squared_adj_significant_results.csv"))
  write.csv(significant_data, file = results_save_path, row.names = FALSE)

  # =========================================================================
  # PART 3: PRINT ASCII HISTOGRAM TO MONITOR & SAVE TEXT FILE
  # =========================================================================
  text_report_path <- file.path(motif_path, paste0(motif_name, "_r_squared_adj_ascii_distribution.txt"))
  report_lines <- c(
    "----------------------------------------",
    paste("  R² ADJ DISTRIBUTION FOR:", motif_name),
    "----------------------------------------"
  )

  if (nrow(plot_data) > 0) {
    bins <- cut(plot_data$r_squared_adj, breaks = seq(-1, 1, by = 0.2), include.lowest = TRUE)
    counts <- table(bins)
    max_count <- max(counts)
    scale_factor <- if (max_count > 40) 40 / max_count else 1

    for (bin_name in names(counts)) {
      count_val <- counts[bin_name]
      bars <- paste0(rep("█", round(count_val * scale_factor)), collapse = "")
      report_lines <- c(report_lines, sprintf("%-14s | %s (%d)", bin_name, bars, count_val))
    }
  } else {
    report_lines <- c(report_lines, "No data available to plot.")
  }

  # Print the distribution to the terminal
  cat(paste(report_lines, collapse = "\n"), "\n")

  # Save text copy
  writeLines(report_lines, con = text_report_path)

  cat("    SUCCESS: 3 new files generated inside this folder.\n\n")
}

cat("========================================\n")
cat(" ALL FOLDERS PROCESSED SUCCESSFULLY!\n")
cat("========================================\n")
