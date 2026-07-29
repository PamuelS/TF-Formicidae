# to launch it properly:

# for folder in 01_enriched_0.25 02_enriched_0.40 03_enriched_0.55 04_enriched_0.70; do
#     echo "Processing $folder..."
#     Rscript GO_enrichment_R2_bounds.R "$folder"
# done



library(tidyverse)
library(topGO)

# 1. Accept the target directory as a command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a target folder path. Usage: Rscript GO_enrichment_R2_bounds.R <folder_path>")
}
target_dir <- args[1]

# Set the working directory to the target folder so outputs save there
setwd(target_dir)

# 2. Load your GO mapping dictionary 
# (Since we are inside the 01_enriched_X folder, "../00_gene_universe/go_back.tsv" will correctly point to the parent's sibling folder)
gene_universe_path <- "../../../00_gene_universe/go_back.tsv"
if(!file.exists(gene_universe_path)) {
  stop(paste("Cannot find gene universe at", gene_universe_path, "- Please check your folder structure."))
}
gene_universe <- readMappings(file = gene_universe_path)
geneUniverse <- names(gene_universe)

# 3. Dynamically find target files in the current folder
# This finds all .txt files but strictly excludes previously generated topGO output files
all_txt_files <- list.files(pattern = "\\.txt$")
target_files <- all_txt_files[!grepl("^topGOe_", all_txt_files)]

if (length(target_files) == 0) {
  stop(paste("No valid target .txt files found in", target_dir))
}

# 4. Automatically load all found target files into a named list
target_list <- list()
for (file in target_files) {
  # Strip the ".txt" extension to use the filename as the trait_name
  trait_name <- tools::file_path_sans_ext(file)
  target_list[[trait_name]] <- read.table(file, header = FALSE)
}

# 5. Your TopGO Functions (Unchanged, just slightly safer character conversion for p-values)
GOenrichment <- function(trait, trait_name) {
  genesOfInterest <- as.character(trait$V1)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse

  print(paste("Running analysis for:", trait_name))

  ontology_values = c("BP", "MF", "CC")

  GOdata_list <- lapply(ontology_values, function(ontology_value) {
    GOdata_name <- paste("GOdata_", ontology_value, sep = "")
    assign(GOdata_name, new("topGOdata", ontology=ontology_value, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene_universe))
  })

  elim_list <- lapply(seq_along(ontology_values), function(i) {
    elim_name <- paste("elim_", ontology_values[i], sep="")
    assign(elim_name, runTest(GOdata_list[[i]], algorithm="elim", statistic="fisher"))
  })

  results_elim <- function(GO_data, elim_data) {
    resulte <- GenTable(GO_data, Classic_Fisher = elim_data, orderBy = "Classic_Fisher", topNodes=1000, numChar=1000)
    # Small fix: topGO sometimes outputs "< 1e-30" as text, which breaks as.numeric. gsub fixes this.
    resulte$Classic_Fisher <- as.numeric(resulte$Classic_Fisher)
    resulte <- subset(resulte, Classic_Fisher < 0.05)
    return(resulte)
  }

  results_elim_list <- lapply(seq_along(ontology_values), function(i) {
    resulte_name <- paste("resulte_", ontology_values[i], sep="")
    assign(resulte_name, envir = .GlobalEnv, results_elim(GOdata_list[[i]], elim_list[[i]]))
  })

  write_elim_results <- function(result, ontology_value, raw_trait_name) {
    trait <- gsub("_", "", gsub("^s", "", raw_trait_name))
    table_name <- paste("./topGOe_", trait, "_", ontology_value, ".txt", sep="")
    write.table(result, file=table_name, quote=F, sep = "\t", row.names = F)
  }

  lapply(seq_along(ontology_values), function(i) {
    write_elim_results(results_elim_list[[i]], ontology_values[i], trait_name)
  })
}

GO_enrichment <- function(list) {
  lapply(seq_along(list), function(i) {
    GOenrichment(list[[i]], names(list)[i])
  })
}

# 6. Execute!
GO_enrichment(target_list)
