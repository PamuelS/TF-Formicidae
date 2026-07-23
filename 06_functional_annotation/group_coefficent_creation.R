# Load required library
library(dplyr)

# 1. Load your significant PGLS results
# (Adjust the column names if yours are slightly different)
sig_data <- read.delim("pgls5_significant_OGs.tsv")

# Set your strictness threshold (e.g., 0.80 means 80% of a motif's OGs must point in the same direction)
bias_threshold <- 0.80

# 2. Calculate the bias for each motif
motif_summary <- sig_data %>%
  group_by(motivo) %>%
  summarise(
    Total_OGs = n(),
    Pos_OGs = sum(Regression_Coefficient > 0),
    Neg_OGs = sum(Regression_Coefficient < 0),
    Polymorphic_Ratio = Pos_OGs / Total_OGs,
    .groups = 'drop'
  )

# 3. Identify the Extreme Motifs
polymorphic_motifs <- motif_summary %>% 
  filter(Polymorphic_Ratio >= bias_threshold) %>% 
  pull(motivo)

monomorphic_motifs <- motif_summary %>% 
  filter(Polymorphic_Ratio <= (1 - bias_threshold)) %>% 
  pull(motivo)

# 4. Extract the unique OGs for each group
# We also ensure we only grab the OGs going in the correct direction for those motifs
polymorphic_OGs <- sig_data %>%
  filter(motivo %in% polymorphic_motifs & Regression_Coefficient > 0) %>%
  pull(OG) %>%
  unique()

monomorphic_OGs <- sig_data %>%
  filter(motivo %in% monomorphic_motifs & Regression_Coefficient < 0) %>%
  pull(OG) %>%
  unique()

# 5. Print a summary to the console
cat("=== MOTIF BIAS SUMMARY ===\n")
cat("Total Motifs Analyzed:", nrow(motif_summary), "\n")
cat("Threshold Used:", bias_threshold * 100, "%\n\n")

cat("POLYMORPHIC GROUP (Far Left):\n")
cat("- Number of Highly Biased Motifs:", length(polymorphic_motifs), "\n")
cat("- Total Unique OGs to test in topGO:", length(polymorphic_OGs), "\n\n")

cat("MONOMORPHIC GROUP (Far Right):\n")
cat("- Number of Highly Biased Motifs:", length(monomorphic_motifs), "\n")
cat("- Total Unique OGs to test in topGO:", length(monomorphic_OGs), "\n\n")

# 6. Save the OG lists as plain text files for your topGO analysis
write.table(polymorphic_OGs, "topGO_target_Polymorphic_OGs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(monomorphic_OGs, "topGO_target_Monomorphic_OGs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
