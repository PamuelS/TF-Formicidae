#!/usr/bin/env Rscript

library(topGO)
library(tidyverse)

# ──────────────────────────────────────────────
# 1. PATHS & SETTINGS
# ──────────────────────────────────────────────
# UPDATE THIS PATH TO YOUR ACTUAL MAPPING FILE!
mapping_file_path <- "../00_gene_universe/go_back.tsv" 

# Read PGLS input data from 05_pgls_analysis
base_pgls_dir <- "/DATASMALL/samuel.pederzini/TF-Formicidae/05_pgls_analysis/01_pgls_results"

# Output directory for enrichment results
base_enrichment_dir <- "/DATASMALL/samuel.pederzini/TF-Formicidae/06_functional_annotation/01_enrichment"

master_enrichment_file <- "master_topGO_enrichment_results.tsv"

# R2 to output folder prefix mapping
r2_to_folder <- c(
  "0.25" = "01_enriched_0.25",
  "0.40" = "02_enriched_0.40",
  "0.50" = "03_enriched_0.50",
  "0.70" = "04_enriched_0.70"
)

# ──────────────────────────────────────────────
# 2. LOAD DATA
# ──────────────────────────────────────────────
# 2.1 Load Gene-to-GO Mapping (Gene Universe)
gene2GO_mapping <- readMappings(file = mapping_file_path)
geneUniverse    <- names(gene2GO_mapping)

# 2.2 Load Master Enrichment Results
if (!file.exists(master_enrichment_file)) {
  stop("Master enrichment file not found: ", master_enrichment_file)
}
enrich_df <- read.delim(master_enrichment_file, sep = "\t", stringsAsFactors = FALSE)

# Deduce Ontology (BP, CC, MF) from the Source_File column
enrich_df <- enrich_df %>%
  mutate(Ontology = case_when(
    grepl("_BP\\.txt", Source_File) ~ "BP",
    grepl("_CC\\.txt", Source_File) ~ "CC",
    grepl("_MF\\.txt", Source_File) ~ "MF",
    TRUE ~ "BP"
  ))

# Identify all unique combinations to process
combinations <- enrich_df %>%
  select(R_squared, Analysis_Type, Target_Group) %>%
  distinct()

# ──────────────────────────────────────────────
# 3. MAIN EXTRACTION LOOP
# ──────────────────────────────────────────────
for (i in 1:nrow(combinations)) {
  
  raw_r2       <- combinations$R_squared[i]
  r2           <- sprintf("%.2f", as.numeric(raw_r2)) # e.g. "0.40"
  ana_type     <- as.character(combinations$Analysis_Type[i])
  target_group <- as.character(combinations$Target_Group[i])
  
  message("\n==================================================")
  message(sprintf(" Processing: R2=%s | %s | %s Bias", r2, ana_type, target_group))
  message("==================================================")
  
  # Determine input PGLS file
  pgls_folder <- file.path(base_pgls_dir, paste0("OG_signif_", r2))
  
  is_pval <- grepl("pval|p-value", tolower(ana_type))
  if (is_pval) {
    pgls_file <- file.path(pgls_folder, "pgls5_significant_OGs_pval.tsv")
    sub_analysis_folder <- "00_pvalue_analysis"
  } else {
    pgls_file <- file.path(pgls_folder, "pgls5_significant_OGs_fdr.tsv")
    sub_analysis_folder <- "01_fdr_analysis"
  }
  
  if (!file.exists(pgls_file)) {
    message("  [WARNING] PGLS file missing, skipping: ", pgls_file)
    next
  }
  
  # Determine exact output path inside 06_functional_annotation/01_enrichment
  parent_enrich_folder <- r2_to_folder[[r2]]
  if (is.null(parent_enrich_folder)) {
    message("  [WARNING] Unknown R2 value, skipping: ", r2)
    next
  }
  
  out_dir <- file.path(base_enrichment_dir, parent_enrich_folder, sub_analysis_folder)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Load target OGs for this specific test
  pgls_data <- read.delim(pgls_file, sep = "\t", stringsAsFactors = FALSE)
  
  if (tolower(target_group) == "polymorphic") {
    genesOfInterest <- pgls_data$OG[pgls_data$Regression_Coefficient > 0]
  } else {
    genesOfInterest <- pgls_data$OG[pgls_data$Regression_Coefficient < 0]
  }
  genesOfInterest <- as.character(genesOfInterest)
  
  # Create logical factor list for topGO
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # Filter enrichment table for current combination
  sub_enrich <- enrich_df %>% 
    filter(R_squared == raw_r2, Analysis_Type == ana_type, Target_Group == target_group)
  
  out_rows <- list()
  
  # Loop over ontologies present in this subset
  for (ont in unique(sub_enrich$Ontology)) {
    ont_enrich <- sub_enrich %>% filter(Ontology == ont)
    
    suppressMessages({
      topGO_obj <- tryCatch(
        new("topGOdata",
            ontology = ont,
            allGenes = geneList,
            annot    = annFUN.gene2GO,
            gene2GO  = gene2GO_mapping),
        error = function(e) { message("  [Error] Failed to build topGO for ", ont); return(NULL) }
      )
    })
    
    if (is.null(topGO_obj)) next
    
    ont_enrich$Mapped_OGs <- NA_character_
    ont_enrich$N_Mapped_OGs <- NA_integer_
    
    for (j in 1:nrow(ont_enrich)) {
      go_id <- ont_enrich$GO.ID[j]
      
      mapped_genes_list <- tryCatch(genesInTerm(topGO_obj, go_id), error = function(e) list())
      
      if (go_id %in% names(mapped_genes_list)) {
        mapped_genes    <- mapped_genes_list[[go_id]]
        intersected_ogs <- intersect(mapped_genes, genesOfInterest)
        
        ont_enrich$Mapped_OGs[j]   <- paste(intersected_ogs, collapse = ",")
        ont_enrich$N_Mapped_OGs[j] <- length(intersected_ogs)
      } else {
        ont_enrich$Mapped_OGs[j]   <- "None"
        ont_enrich$N_Mapped_OGs[j] <- 0
      }
    }
    out_rows[[ont]] <- ont_enrich
  }
  
  if (length(out_rows) > 0) {
    final_combo_df <- bind_rows(out_rows)
    
    clean_dir    <- tolower(target_group)
    out_filename <- sprintf("GO_Mapped_OGs_%s.tsv", clean_dir)
    out_filepath <- file.path(out_dir, out_filename)
    
    write.table(final_combo_df, out_filepath, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  -> Saved mapping to: ", out_filepath)
  }
}

message("\nAll extractions complete!")
