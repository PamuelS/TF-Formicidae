#!/usr/bin/env Rscript

library(caper)
library(ape)
library(dplyr)
library(tidyr)
library(purrr)
library(phytools)
library(pbmcapply)

cat("Starting script ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# ──────────────────────────────────────────────
#0. SETTINGS ARGUMENTS
# ──────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

tree_path <- args[1]
pheno_path <- args[2]
motifs_path <- args[3]
cores_path <- as.numeric(args[4])
output_path <- args[5]

# ──────────────────────────────────────────────
# 1. TREE ----
# ──────────────────────────────────────────────
tree <- read.tree(tree_path)
#tree <- read.tree("abbreviative_name_rooted_Leptsp_tree.txt")

# ──────────────────────────────────────────────
# 2. PHENOTYPE TABLE ----
#    Columns used:
#      Abbreviation         → species ID (must match tree tip labels)
#      Number of queens     → pgls6a / pgls6b
#      Castes               → pgls5
#      Pasassitism          → pgls1
#      Queen wings          → pgls4a / pgls4b
#      Larval hemolymph feeding → pgls2
#      Cocoon               → pgls3
# ──────────────────────────────────────────────
pheno_raw <- read.table(
  pheno_path,
  #"updated_dataset.tsv",
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE,
  encoding    = "UTF-8"
)

# Keep only columns we need and rename them for clarity
pheno <- pheno_raw %>%
  select(
    species       = "Abbreviation",
    queens_raw    = "Number of queens",
    castes_raw    = "Castes",
    parasit_raw   = "Pasassitism",
    wings_raw     = "Queen wings",
    larval_raw    = "Larval hemolymph feeding",
    cocoon_raw    = "Cocoon"
  ) %>%
  mutate(across(everything(), ~ trimws(as.character(.))))

# ── Italian → English translation + factor coding ──────────────────────────
translate_binary <- function(x) {
  dplyr::case_when(
    grepl("^[Aa]ss",   x) ~ "absent",
    grepl("^[Pp]res",  x) ~ "present",
    TRUE ~ NA_character_
  )
}

translate_wings <- function(x) {
  dplyr::case_when(
    grepl("[Aa]lata",      x) ~ "winged",
    grepl("[Ee]rgatoide",  x) ~ "ergatoid",
    grepl("[Bb]rachi",     x) ~ "brachypterous",
    TRUE ~ NA_character_
  )
}

translate_castes <- function(x) {
  dplyr::case_when(
    grepl("[Mm]ono",  x) ~ "monomorphic",   
    grepl("[Dd]imor", x) ~ "polymorphic",   
    grepl("[Pp]oli",  x) ~ "polymorphic",   
    grepl("[Ii]nes",  x) ~ NA_character_,   
    TRUE ~ NA_character_
  )
}

translate_queens <- function(x) {
  dplyr::case_when(
    grepl("[Mm]onogi",    x) ~ "monogyne",
    grepl("[Pp]oligi",    x) ~ "polygyne",
    grepl("[Ff]acol",     x) ~ "facultative",
    grepl("[Gg]amergate", x) ~ "monogyne",   
    TRUE ~ NA_character_
  )
}

pheno <- pheno %>%
  mutate(
    parasit  = translate_binary(parasit_raw),
    larval   = translate_binary(larval_raw),
    cocoon   = translate_binary(cocoon_raw),
    wings    = translate_wings(wings_raw),
    castes   = translate_castes(castes_raw),
    queens   = translate_queens(queens_raw)
  ) |> 
  select(-c(parasit_raw, larval_raw, cocoon_raw,
           wings_raw, castes_raw, queens_raw))

# ──────────────────────────────────────────────
# 3. OUTLIER HANDLING & NORMALISATION ----
# ──────────────────────────────────────────────
handle_outliers <- function(y) {
  x        <- y[!is.na(y)]
  med      <- median(x)
  #' MAD -> outlier-resistant measure of variability
  mad_val  <- median(abs(x - med))
  mean_ad_val <- mean(abs(x - med))
  threshold <- if (mad_val != 0) 2.5 * mad_val else 2.5 * mean_ad_val
  #' Winsorizing. Statistical data-cleaning technique that limits
  #' the impact of extreme outliers by capping them to the most extreme value
  #' within the nearest non-extreme values
  #' Observation are valid even if outliers
  x[x > med + threshold] <- med + threshold
  x[x < med - threshold] <- med - threshold
  y[!is.na(y)] <- x
  return(y)
}

normalize <- function(y) {
  x <- y[!is.na(y)]
  x <- (x - min(x)) / (max(x) - min(x))
  y[!is.na(y)] <- x
  return(y)
}

# ──────────────────────────────────────────────
# 4. CORE PGLS RUNNER ----
# ──────────────────────────────────────────────
bounds_scenarios <- list(c(1e-05, 1), c(1e-03, 1), c(1e-01, 1), c(1, 1))

run_pgls_with_bound <- function(data, bounds, tree,
                                species_pheno, phenotype_col) {
  OG     <- rownames(data)
  #' Make scores as a vertical df not horizontal
  scores_df <- data.frame(species = colnames(data),
                          score = as.numeric(data))
  #' Add the phenotype
  merged_data <- merge(scores_df, species_pheno, by = "species", all.x = TRUE)
  #' Remove levels that are no longer present among ours
  merged_data[[phenotype_col]] <- droplevels(merged_data[[phenotype_col]])
  
  #' Create the comparative data object that merge phy with the phenotype and
  #' score df
  comparative_data <- comparative.data(
    phy       = tree,
    data      = merged_data,
    names.col = species,
    vcv       = TRUE
  )

  #' Create the formula based on the df used
  formula_pgls <- as.formula(paste("score ~", phenotype_col))
  #' Run pgls for each bound given
  for (bound in bounds) {
    pgls_result <- tryCatch(
      pgls(formula_pgls, data = comparative_data,
           lambda = "ML", bounds = list(lambda = bound)),
      error = function(e) NULL
    )
    #' If success write extract values and write the output df
    if (!is.null(pgls_result)) {
      s <- summary(pgls_result)
      r_sq     <- s$r.squared
      n        <- nrow(merged_data)
      #' p is the number of contrast obtained as nlevels - 1
      p        <- length(levels(merged_data[[phenotype_col]])) - 1
      r_sq_adj <- 1 - ((1 - r_sq) * (n - 1) / (n - p - 1))
      coef_tbl <- s$coefficients

      coef_rows <- coef_tbl[-1, , drop = FALSE]

      return(data.frame(
        OG                    = OG,
        Regression_Coefficient = coef_rows[, 1],
        P_value               = coef_rows[, 4],
        Lambda                = s$param[2],
        Bound_used            = paste(bound[1], "to", bound[2]),
        R_squared             = r_sq,
        R_squared_adj         = r_sq_adj,
        n                     = n,
        Status                = "success",
        Contrast              = rownames(coef_rows),
        stringsAsFactors      = FALSE
      ))
    }
  }

  return(data.frame(
    OG = OG, Regression_Coefficient = NA, P_value = NA,
    Lambda = NA, Bound_used = "None",
    R_squared = NA, R_squared_adj = NA, n = NA,
    Status = "fail", Contrast = NA,
    stringsAsFactors = FALSE
  ))
}

# ──────────────────────────────────────────────
# 5. PROCESS ONE MOTIF FILE FOR ONE CONTRAST ----
# ──────────────────────────────────────────────
process_motif_file <- function(file, tree, species_pheno, phenotype_col,
                               min_per_group, cores) {
  #' Read the motif table and make Og the rowname
  data_raw <- read.table(file, header = TRUE, sep = "\t")
  rownames(data_raw) <- data_raw$index
  data_raw$index     <- NULL

  #' Probably not needed at least when elaborating DISCO OGs
  #data_raw_elaborated[] <- lapply(data_raw, function(col)
  #  if (is.numeric(col)) handle_outliers(col) else col)
  #' min-maz normalisation
  data_raw[] <- lapply(data_raw, function(col)
    if (is.numeric(col)) normalize(col) else col)
  data_raw <- as.data.frame(data_raw)

  #' Check if species in the dataframe match species in the
  #' phenotype vector (neeed for castes)
  valid_species <- intersect(
    colnames(data_raw),
    species_pheno$species[!is.na(species_pheno[[phenotype_col]])]
  )
  data_raw <- data_raw[, valid_species, drop = FALSE]

  #' Count nspecies per OG grouped by phenotype level
  levels_vec <- levels(species_pheno[[phenotype_col]])
  counts_per_group <- lapply(levels_vec, function(lv) {
    sp_lv <- species_pheno$species[
      !is.na(species_pheno[[phenotype_col]]) &
        species_pheno[[phenotype_col]] == lv
    ]
    sp_lv <- intersect(sp_lv, colnames(data_raw))
    rowSums(!is.na(data_raw[, sp_lv, drop = FALSE]))
  })
  names(counts_per_group) <- levels_vec

  #' Check if OG matches min threshold of phenotype presence
  keep <- Reduce("&", mapply(function(cnt, minv) cnt >= minv,
                             counts_per_group, min_per_group,
                             SIMPLIFY = FALSE))

  data_filt <- data_raw[keep, , drop = FALSE]

  data_filt <- data_filt[
    apply(data_filt, 1, function(x) var(na.omit(x)) > 0), ,
    drop = FALSE
  ]

  if (nrow(data_filt) == 0) {
    message("No rows passed filters for phenotype: ", phenotype_col)
    return(data.frame())
  }

  #' Apply the pgls function using parallelisation based on row indices not col
  results <- pbmclapply(seq_len(nrow(data_filt)), function(i) {
    #' Extract single OG
    single_og <- data_filt[i, , drop = FALSE]
    run_pgls_with_bound(
      single_og, bounds_scenarios, tree, species_pheno, phenotype_col)
  }, ignore.interactive = getOption("ignore.interactive", T), mc.cores = cores)
  
  results <- bind_rows(results)
}

# ──────────────────────────────────────────────
# 6. DEFINE THE 8 CONTRASTS ----
# ──────────────────────────────────────────────
build_sp_pheno <- function(species_vec, values_vec, col_name,
                           levels_order, ordered = FALSE) {
  df <- data.frame(
    species   = species_vec,
    phenotype = factor(values_vec, levels = levels_order, ordered = ordered),
    stringsAsFactors = FALSE
  )
  colnames(df)[2] <- col_name
  df <- df[!is.na(df[[col_name]]), ]
  df
}

#sp_pgls1  <- build_sp_pheno(pheno$species, pheno$parasit, "pgls1_parasitism",
#                            levels_order = c("absent", "present"))
#sp_pgls2  <- build_sp_pheno(pheno$species, pheno$larval, "pgls2_larval",
#                            levels_order = c("absent", "present"))
#sp_pgls3  <- build_sp_pheno(pheno$species, pheno$cocoon, "pgls3_cocoon", 
#                            levels_order = c("absent", "present"))
#sp_pgls4a <- build_sp_pheno(pheno$species, pheno$wings, "pgls4a_wings",
#                            levels_order = c("winged", "ergatoid"))
#sp_pgls4b <- build_sp_pheno(pheno$species, pheno$wings, "pgls4b_wings",
#                            levels_order = c("winged", "brachypterous"))
sp_pgls5  <- build_sp_pheno(pheno$species, pheno$castes, "pgls5_castes",
                            levels_order = c("monomorphic", "polymorphic"))

#sp_pgls6a <- build_sp_pheno(pheno$species, pheno$queens, "pgls6a_queens",
#                            levels_order = c("monogyne", "polygyne"))
#sp_pgls6a <- sp_pgls6a[sp_pgls6a$pgls6a_queens %in% c("monogyne", "polygyne"), ]
#sp_pgls6a$pgls6a_queens <- droplevels(sp_pgls6a$pgls6a_queens)

#sp_pgls6b <- build_sp_pheno(pheno$species, pheno$queens, "pgls6b_queens",
#                            levels_order = c(
#                              "monogyne", "facultative", "polygyne"),
#                            ordered = TRUE)

# ──────────────────────────────────────────────
# 7. MINIMUM SPECIES PER GROUP ----
# ──────────────────────────────────────────────
min_binary    <- c(absent = 3, present = 2)
min_wings4a   <- c(winged = 3, ergatoid = 2)
min_wings4b   <- c(winged = 3, brachypterous = 2)
min_castes    <- c(monomorphic = 3, polymorphic = 3)
min_queens6a  <- c(monogyne = 3, polygyne = 2)
min_queens6b  <- c(monogyne = 2, facultative = 2, polygyne = 2)

# ──────────────────────────────────────────────
# 8. INPUT MOTIF FILE (from Snakemake) ----
# ──────────────────────────────────────────────
motif_file <- motifs_path
#motif_file <- "totalscore_MA0010.2.tsv"

# ──────────────────────────────────────────────
# 9. RUN ALL 8 ANALYSES ----
# ──────────────────────────────────────────────
cores <- cores_path

#message("Running pgls1 : Parasitism")
#res_pgls1 <- process_motif_file(motif_file, tree, sp_pgls1,
#                                "pgls1_parasitism", min_binary)

#message("Running pgls2 : Larval haemolymph feeding")
#res_pgls2 <- process_motif_file(motif_file, tree, sp_pgls2,
#                                "pgls2_larval", min_binary)

#message("Running pgls3 : Cocoon")
#res_pgls3 <- process_motif_file(motif_file, tree, sp_pgls3,
#                                "pgls3_cocoon", min_binary)

#message("Running pgls4a : Queen wings (winged vs ergatoid)")
#res_pgls4a <- process_motif_file(motif_file, tree, sp_pgls4a,
#                                 "pgls4a_wings", min_wings4a)

#message("Running pgls4b : Queen wings (winged vs brachypterous)")
#res_pgls4b <- process_motif_file(motif_file, tree, sp_pgls4b,
#                                 "pgls4b_wings", min_wings4b)

message("Running pgls5 : Castes")
res_pgls5 <- process_motif_file(motif_file, tree, sp_pgls5,
                                "pgls5_castes", min_castes, cores)

#message("Running pgls6a : Number of queens (monogyne vs polygyne)")
#res_pgls6a <- process_motif_file(motif_file, tree, sp_pgls6a,
#                                 "pgls6a_queens", min_queens6a)

#message(paste0("Running pgls6b : Number of queens ",
#        "(ordered: monogyne < facultative < polygyne)"))
#res_pgls6b <- process_motif_file(motif_file, tree, sp_pgls6b,
#                                 "pgls6b_queens", min_queens6b)

# ──────────────────────────────────────────────
# 10. WRITE OUTPUT FILES ----
# ──────────────────────────────────────────────
#write.csv(res_pgls1,  snakemake@output[["pgls1"]],  row.names = FALSE)
#write.csv(res_pgls2,  snakemake@output[["pgls2"]],  row.names = FALSE)
#write.csv(res_pgls3,  snakemake@output[["pgls3"]],  row.names = FALSE)
#write.csv(res_pgls4a, snakemake@output[["pgls4a"]], row.names = FALSE)
#write.csv(res_pgls4b, snakemake@output[["pgls4b"]], row.names = FALSE)
write.csv(res_pgls5,  output_path,  row.names = FALSE)
#write.csv(res_pgls6a, snakemake@output[["pgls6a"]], row.names = FALSE)
#write.csv(res_pgls6b, snakemake@output[["pgls6b"]], row.names = FALSE)

message("All 8 PGLS analyses complete.")
cat("Ending script ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
