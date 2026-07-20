#!/usr/bin/env Rscript


#this script allows to eliminate and filter every OGs presenti in every motif pgls file that hasn't reached the significant limit of the p-value and the R^2 adjusted


library(tidyverse)

# ──────────────────────────────────────────────
# SETTINGS
# ──────────────────────────────────────────────
input_dir   <- "."          # cartella con le sottocartelle MAXXX.X/
p_thresh    <- 0.05         # Soglia per P-value grezzo
r2_thresh   <- 0.25


# ──────────────────────────────────────────────
# 1. LEGGI TUTTI I FILE
# ──────────────────────────────────────────────
files <- list.files(
  path       = input_dir,
  pattern    = "pgls5_castes\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(files) == 0) stop("Nessun file pgls5_castes.csv trovato.")

message("File trovati: ", length(files))

df_all <- map_dfr(files, function(f) {
  motivo <- basename(dirname(f))   # es. "MA0001.1"
  dat    <- read.csv(f, header = TRUE)
  dat$motivo <- motivo
  dat
})

message("OG totali caricati: ", nrow(df_all))

# ──────────────────────────────────────────────
# 2. RIMUOVI FAIL E NA
# ──────────────────────────────────────────────
df_clean <- df_all %>%
  filter(Status == "success") %>%
  filter(!is.na(P_value), !is.na(R_squared_adj))

message("OG con status success e p-value non-NA: ", nrow(df_clean))

# ──────────────────────────────────────────────
# 3. FILTRA: P-value grezzo < 0.05 e R² adj > 0.25
# ──────────────────────────────────────────────
df_sig <- df_clean %>%
  filter(P_value       < p_thresh,
         R_squared_adj > r2_thresh) %>%
  arrange(motivo, P_value)

message("OG significativi dopo filtro: ", nrow(df_sig))
message("Motivi con almeno un OG significativo: ",
        n_distinct(df_sig$motivo))

# ──────────────────────────────────────────────
# 4. OUTPUT
# ──────────────────────────────────────────────

output_dir  <- "OG_pvalue_Rsquaredadj_signif"

# Crea la cartella principale se non esiste
if (!dir.exists(output_dir)) dir.create(output_dir)

# 4a. Tabella globale completa dei significativi
write.table(df_sig,
          file      = file.path(output_dir, "pgls5_significant_OGs.tsv"),
          row.names = FALSE)

# 4b. Sommario globale per motivo
summary_motivo <- df_sig %>%
  group_by(motivo) %>%
  summarise(
    n_OG_sig          = n(),
    median_R2_adj     = round(median(R_squared_adj), 4),
    median_pval       = round(median(P_value), 4),
    n_positive_coef   = sum(Regression_Coefficient > 0),
    n_negative_coef   = sum(Regression_Coefficient < 0),
    direction_dominant = ifelse(n_positive_coef >= n_negative_coef,
                                "polymorphic", "monomorphic"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_OG_sig))

write.table(summary_motivo,
          file      = file.path(output_dir, "pgls5_summary_by_motif.tsv"),
          row.names = FALSE)

# 4c. NUOVO: Salva un singolo file CSV per ogni motivo estratto
motif_folder <- file.path(output_dir, "single_motif")
if (!dir.exists(motif_folder)) dir.create(motif_folder)

message("Scrittura dei singoli file per motivo...")

df_sig %>%
  group_split(motivo) %>%
  walk(function(sub_df) {
    # Prende il nome del motivo dalla prima riga del sotto-gruppo
    motif_name <- sub_df$motivo[1]

    # Costruisce il percorso finale (es. "00_OG_significant/single_motif/MA0001.1.tsv")
    file_out <- file.path(motif_folder, paste0(motif_name, ".tsv"))

    # Salva il file per questo specifico motivo
    write.table(sub_df, file = file_out, row.names = FALSE)
  })

message("Processo completato con successo!")
message("→ Controlla la cartella: ", motif_folder)
