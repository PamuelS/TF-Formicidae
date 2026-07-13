#!/usr/bin/env Rscript
# =============================================================================
# motif_disco_vs_ortho_plot.R
#
# Analisi DISCO vs Orthofinder — versione R / ggplot2
#
# Utilizzo:
#   Rscript motif_disco_vs_ortho_plot.R \
#       --results-dir risultati_ortho_vs_disco \
#       --outdir     output_R
#
# Input atteso in --results-dir:
#   <results-dir>/
#     batch_summary.tsv          (prodotto da Python: concatenazione summary_stats per motivo)
#     <motif_name>/
#       per_og_comparison.tsv    (uno per motivo)
#       summary_stats.tsv        (uno per motivo)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(stringr)
  library(forcats)
  library(gridExtra)
  library(cowplot)
  library(grid)
})

# ─────────────────────────────────────────────────────────────────────────────
# 0. Argomenti da riga di comando
# ─────────────────────────────────────────────────────────────────────────────

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    results_dir = "risultati_ortho_vs_disco",
    outdir      = "output_R"
  )
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--results-dir" = { opts$results_dir <- args[i+1]; i <- i+2 },
      "--outdir"      = { opts$outdir      <- args[i+1]; i <- i+2 },
      { cat("Argomento sconosciuto:", args[i], "\n"); i <- i+1 }
    )
  }
  opts
}

opts <- parse_args()
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
cat("=== DISCO vs Orthofinder — analisi R ===\n")
cat(sprintf("  results-dir : %s\n", opts$results_dir))
cat(sprintf("  outdir      : %s\n\n", opts$outdir))

# ─────────────────────────────────────────────────────────────────────────────
# 1. Caricamento batch_summary.tsv
# ─────────────────────────────────────────────────────────────────────────────

batch_path <- file.path(opts$results_dir, "batch_summary.tsv")
if (!file.exists(batch_path))
  stop("batch_summary.tsv non trovato in: ", opts$results_dir)

batch <- read_tsv(batch_path, show_col_types = FALSE)
cat(sprintf("batch_summary caricato: %d righe, motivi = %d\n",
            nrow(batch), n_distinct(batch$motif)))

# Sottotabelle per categoria
sc <- batch %>% filter(category == "single_complete")
sp <- batch %>% filter(category == "split")
oo <- batch %>% filter(category == "ortho_only")

all_motifs <- sort(unique(batch$motif))
n_motifs   <- length(all_motifs)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Aggregazione dei per_og_comparison.tsv per grafici per-motivo
# ─────────────────────────────────────────────────────────────────────────────

cat("Lettura per_og_comparison.tsv per motivo...\n")

per_og_list <- list()
motif_dirs <- list.dirs(opts$results_dir, full.names = TRUE, recursive = FALSE)

for (md in motif_dirs) {
  motif_name <- basename(md)
  per_og_path <- file.path(md, "per_og_comparison.tsv")
  if (!file.exists(per_og_path)) next
  df <- tryCatch(
    read_tsv(per_og_path, show_col_types = FALSE) %>% mutate(motif = motif_name),
    error = function(e) NULL
  )
  if (!is.null(df)) per_og_list[[motif_name]] <- df
}

cat(sprintf("  Letti %d file per_og_comparison.tsv\n", length(per_og_list)))

if (length(per_og_list) == 0)
  stop("Nessun per_og_comparison.tsv trovato. Controlla --results-dir.")

per_og_all <- bind_rows(per_og_list)

if ("og_base" %in% names(per_og_all)) {
  per_og_all <- per_og_all %>% rename(og = og_base)
} else {
  colnames(per_og_all)[1] <- "og"
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Calcolo statistiche per pannello [A] — riepilogo globale
# ─────────────────────────────────────────────────────────────────────────────

uniform_or_median <- function(x, label = "") {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("n.d.")
  uv <- unique(x)
  if (length(uv) == 1) return(as.character(as.integer(uv)))
  sprintf("%d (mediana)", as.integer(median(x)))
}

n_raw_val  <- uniform_or_median(sc$n_og_raw_ortho)
n_sc_val   <- uniform_or_median(sc$n_og)
n_sp_val   <- uniform_or_median(sp$n_og)
n_oo_val   <- uniform_or_median(oo$n_og)

mw_col <- "mw_pval_sc_vs_sp"
mw_sig_motifs <- character(0)
if (mw_col %in% names(sc)) {
  mw_sig_motifs <- sc %>%
    filter(!is.na(.data[[mw_col]]), .data[[mw_col]] < 0.05) %>%
    pull(motif) %>% sort()
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. Palette colori e tema ggplot condiviso
# ─────────────────────────────────────────────────────────────────────────────

COL_SC   <- "#4C8FBF"
COL_SP   <- "#E8703A"
COL_OO   <- "#6AAF6A"
COL_LOSS <- "#C0392B"
COL_OK   <- "#E8703A"

theme_disco <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 1,
                                     margin = margin(b = 6)),
      plot.subtitle   = element_text(size = base_size - 1, color = "#555555",
                                     margin = margin(b = 4)),
      axis.title      = element_text(size = base_size),
      axis.text       = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = base_size - 1)
    )
}

# =============================================================================
# 5. Pannello [A] — Boxplot ratio DISCO/ORTHO
# =============================================================================

# =============================================================================
# 5. Pannello [A] — Boxplot ratio DISCO/ORTHO
# =============================================================================

make_panel_A <- function() {

  ratio_long <- bind_rows(
    sc %>% select(motif, median_ratio_sum, mw_pval_sc_vs_sp) %>%
      mutate(cat_label = "Single-complete"),
    sp %>% select(motif, median_ratio_sum, mw_pval_sc_vs_sp) %>%
      mutate(cat_label = "Split (\u22652 forme)")
  ) %>% filter(!is.na(median_ratio_sum))

  med_sc <- median(sc$median_ratio_sum, na.rm = TRUE)
  med_sp <- median(sp$median_ratio_sum, na.rm = TRUE)
  n_sig  <- if ("mw_pval_sc_vs_sp" %in% names(sc))
              sum(sc$mw_pval_sc_vs_sp < 0.05, na.rm = TRUE) else 0
  n_tot  <- n_distinct(ratio_long$motif)

  # --- MODIFICA APPLICATA QUI ---
  ann_df <- data.frame(
    cat_label = c("Single-complete", "Split (\u22652 forme)"),
    label     = c(
      sprintf("mediana = %.3f\n%d motivi MW p<0.05", med_sc, n_sig),
      sprintf("mediana = %.3f\n ", med_sp) # Aggiunto "\n " per simulare la seconda riga
    ),
    y         = c(med_sc, med_sc) # Usa med_sc per allineare entrambe alla stessa altezza
  )

  ggplot(ratio_long, aes(x = cat_label, y = median_ratio_sum, fill = cat_label)) +
    geom_violin(alpha = 0.18, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.35, alpha = 0.85, outlier.shape = 21,
                 outlier.size = 1.6, outlier.alpha = 0.5,
                 outlier.fill = "white", linewidth = 0.5) +
    stat_summary(fun = median, geom = "point", shape = 23,
                 size = 4, fill = "white", color = "#333333", stroke = 1.1) +
    geom_hline(yintercept = 1.0, linetype = "dashed",
               color = "#555555", linewidth = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dotted",
               color = COL_LOSS, linewidth = 0.8) +
    geom_text(data = ann_df,
              aes(x = cat_label, y = y, label = label),
              vjust = -1.2, hjust = 0.5, size = 3.5,
              color = "#333333", fontface = "italic",
              inherit.aes = FALSE) +
    annotate("text", x = 2.45, y = 1.01, label = "ratio = 1\n(nessuna perdita)",
             hjust = 1, vjust = 0, size = 3.2, color = "#555555") +
    annotate("text", x = 2.45, y = 0.51, label = "ratio = 0.5",
             hjust = 1, vjust = 0, size = 3.2, color = COL_LOSS) +
    scale_fill_manual(
      values = c("Single-complete" = COL_SC, "Split (\u22652 forme)" = COL_SP)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.18)),
      breaks = seq(0, 2, by = 0.25)
    ) +
    scale_x_discrete(expand = expansion(add = 0.6)) +
    labs(
      title    = "[A]  Ratio DISCO/ORTHO — tutti i motivi",
      subtitle = sprintf(
        "Ogni punto = mediana di un motivo (n = %d)  |  Losanga bianca = mediana globale",
        n_tot
      ),
      x = NULL,
      y = "Ratio mediano (DISCO sum / ORTHO sum)"
    ) +
    theme_disco(base_size = 11) +
    theme(
      legend.position    = "none",
      panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.4),
      plot.subtitle      = element_text(size = 9.5, color = "#555555")
    )
}

# =============================================================================
# 6. Pannello [B] — Stacked bar proporzionale: distribuzione OG budget
# =============================================================================

make_panel_B <- function() {

  n_raw      <- as.integer(str_extract(n_raw_val,  "^[0-9]+"))
  n_sc       <- as.integer(str_extract(n_sc_val,   "^[0-9]+"))
  n_sp       <- as.integer(str_extract(n_sp_val,   "^[0-9]+"))
  n_oo       <- as.integer(str_extract(n_oo_val,   "^[0-9]+"))

  total <- n_sc + n_sp + n_oo

  seg_df <- tibble(
    gruppo = factor(c(
      "Ortho-only",
      "Split (>=2 forme)",
      "Single-complete"
    ), levels = c(
      "Ortho-only",
      "Split (>=2 forme)",
      "Single-complete"
    )),
    n      = c(n_oo, n_sp, n_sc),
    colore = c(COL_OO, COL_SP, COL_SC)
  ) %>%
    mutate(
      pct       = n / total,
      pct_lab   = sprintf("%.1f%%", pct * 100),
      n_lab     = format(n, big.mark = ".", scientific = FALSE),
      leg_label = sprintf("%s — n = %s (%s)", gruppo, n_lab, pct_lab),
      order_idx = row_number()
    )

  n_side      <- 10
  n_cells_tot <- n_side * n_side

  cells_alloc <- seg_df %>%
    mutate(
      raw_cells  = pct * n_cells_tot,
      base_cells = floor(raw_cells),
      remainder  = raw_cells - base_cells
    )

  deficit <- n_cells_tot - sum(cells_alloc$base_cells)
  if (deficit > 0) {
    idx_top <- order(cells_alloc$remainder, decreasing = TRUE)[seq_len(deficit)]
    cells_alloc$base_cells[idx_top] <- cells_alloc$base_cells[idx_top] + 1
  }

  cells_alloc <- cells_alloc %>% arrange(order_idx)

  gruppo_seq <- rep(as.character(cells_alloc$gruppo), times = cells_alloc$base_cells)
  length(gruppo_seq) <- n_cells_tot
  if (anyNA(gruppo_seq)) {
    gruppo_seq[is.na(gruppo_seq)] <- as.character(
      cells_alloc$gruppo[which.max(cells_alloc$base_cells)]
    )
  }

  waffle_df <- tibble(
    idx    = seq_len(n_cells_tot),
    col    = ((idx - 1) %% n_side) + 1,
    row    = ((idx - 1) %/% n_side) + 1,
    y      = n_side - row + 1,
    gruppo = factor(gruppo_seq, levels = levels(seg_df$gruppo))
  )

  cell_value <- total / n_cells_tot

  ggplot(waffle_df, aes(x = col, y = y, fill = gruppo)) +
    geom_tile(color = "white", linewidth = 1.1, width = 0.88, height = 0.88) +
    scale_fill_manual(
      values = setNames(as.character(seg_df$colore), as.character(seg_df$gruppo)),
      labels = setNames(seg_df$leg_label, as.character(seg_df$gruppo)),
      breaks = levels(seg_df$gruppo),
      guide  = guide_legend(ncol = 1, override.aes = list(color = NA))
    ) +
    coord_equal(expand = FALSE, clip = "off") +
    scale_x_continuous(limits = c(0.4, n_side + 0.6)) +
    scale_y_continuous(limits = c(0.4, n_side + 1.3)) +
    annotate(
      "text",
      x = (n_side + 1) / 2, y = n_side + 1.05,
      label = sprintf("Totale OG ORTHO: %s",
                      format(total, big.mark = ".", scientific = FALSE)),
      size = 3.8, fontface = "bold", color = "#333333", hjust = 0.5
    ) +
    labs(
      title    = "[B]  Distribuzione degli OG ORTHO",
      subtitle = sprintf(
        "La strutttura rimane invariata per 296 motivi",
        n_motifs, 100 / n_cells_tot,
        format(round(cell_value), big.mark = ".", scientific = FALSE)
      ),
      fill = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0,
                                     margin = margin(b = 6)),
      plot.subtitle   = element_text(size = 9.5, color = "#555555",
                                     margin = margin(b = 8)),
      legend.position = "right",
      legend.text     = element_text(size = 9),
      legend.key.size = unit(0.55, "cm"),
      plot.margin     = margin(t = 6, r = 6, b = 6, l = 6)
    )
}

# =============================================================================
# 7. Pannello [C] — Perdita copertura specie
# =============================================================================

make_panel_C <- function() {

  delta_df <- sp %>%
    mutate(
      delta      = median_ortho_frac_nz - median_disco_frac_nz,
      perdita    = delta > 0,
      motif_ord  = fct_reorder(motif, delta, .desc = FALSE)
    ) %>%
    filter(!is.na(delta))

  n_loss  <- sum(delta_df$perdita)
  n_ok    <- sum(!delta_df$perdita)
  n_valid <- nrow(delta_df)

  subtitle_str <- sprintf(
    "Perdita (DISCO < ORTHO): %d/%d motivi",
    n_loss, n_valid, n_ok, n_valid
  )

  ggplot(delta_df, aes(x = motif_ord, y = delta, fill = perdita)) +
    geom_col(width = 0.85, alpha = 0.88) +
    geom_hline(yintercept = 0, linewidth = 0.7, color = "#333333") +
    scale_fill_manual(
      values = c("TRUE" = COL_LOSS, "FALSE" = COL_OK),
      labels = c("TRUE" = "Perdita copertura (DISCO < ORTHO)",
                 "FALSE" = "Nessuna perdita (DISCO \u2265 ORTHO)")
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0.05, 0.10))
    ) +
    annotate(
      "text",
      x = 1, y = max(delta_df$delta, na.rm = TRUE) * 0.93,
      label = subtitle_str,
      hjust = 0, vjust = 1, size = 3.8, color = COL_LOSS, fontface = "italic"
    ) +
    labs(
      title    = "[C]  Split \u22652 forme \u2014 perdita di copertura specie per motivo",
      subtitle = paste0(
        "Delta = frac specie con punteggio > 0 in ORTHO  \u2212  frac in DISCO  ",
        "(barre rosse: DISCO copre meno specie di ORTHO)"
      ),
      x = sprintf("Motivi (n = %d, ordinati per delta crescente)", n_valid),
      y = "Delta frac specie > 0  (ORTHO \u2212 DISCO)"
    ) +
    theme_disco(base_size = 12) +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.4),
      legend.position    = "bottom",
      legend.key.size    = unit(0.55, "cm"),
      legend.text        = element_text(size = 11),
      
      plot.title         = element_text(margin = margin(b = 2)),
      
      plot.subtitle      = element_text(margin = margin(b = 2))
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# 8. Composizione finale e salvataggio PDF + PNG
# ─────────────────────────────────────────────────────────────────────────────

cat("Generazione grafici...\n")

panel_A <- make_panel_A()
panel_B <- make_panel_B()
panel_C <- make_panel_C()

fig_w <- max(22, n_motifs * 0.055 + 6)
fig_h <- 22

out_pdf <- file.path(opts$outdir, "batch_overview.pdf")
out_png <- file.path(opts$outdir, "batch_overview.png")

save_figure <- function(device, path, ...) {
  tryCatch({
    device(path, ...)

    top_row <- plot_grid(
      panel_A, panel_B,
      ncol = 2, align = "h", axis = "tb",
      rel_widths = c(1, 1)
    )

    # rel_heights: pannello [C] cresce proporzionalmente con il numero di motivi
    # (più motivi → barre più strette → serve più spazio verticale relativo).
    # Formula: top row fissa al 35–45%, bottom row prende il resto.
    h_top <- max(0.30, min(0.45, 0.50 - n_motifs * 0.0005))
    h_bot <- 1 - h_top

    full <- plot_grid(
      top_row, panel_C,
      ncol = 1,
      rel_heights = c(h_top, h_bot)
    )

    title_grob <- ggdraw() +
      draw_label(
        "Confronto DISCO vs Orthofinder — riepilogo globale tutti i motivi",
        fontface = "bold", size = 14, x = 0.5, hjust = 0.5
      )

    final <- plot_grid(title_grob, full, ncol = 1, rel_heights = c(0.04, 0.96))
    print(final)

    dev.off()
    cat(sprintf("  Salvato: %s\n", path))
  }, error = function(e) {
    cat(sprintf("  ERRORE salvataggio %s: %s\n", path, e$message))
  })
}

save_figure(
  function(p, ...) pdf(p, ...),
  out_pdf,
  width = fig_w, height = fig_h, paper = "special"
)

save_figure(
  function(p, ...) png(p, ...),
  out_png,
  width = fig_w, height = fig_h, units = "in", res = 180
)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Grafici aggiuntivi per motivo (opzionale)
# ─────────────────────────────────────────────────────────────────────────────

cat("\nGenerazione grafici per-motivo (ratio + scatter)...\n")

motif_pdf <- file.path(opts$outdir, "per_motif_plots.pdf")
pdf(motif_pdf, width = 14, height = 6, paper = "special")

for (mot in all_motifs) {
  df_m <- per_og_all %>% filter(motif == mot)
  if (nrow(df_m) == 0) next

  df_plot <- df_m %>%
    filter(category %in% c("single_complete", "split")) %>%
    mutate(cat_label = recode(category,
      single_complete = "Single-complete",
      split           = "Split (\u22652 forme)"
    ))

  p_ratio <- ggplot(df_plot, aes(x = cat_label, y = ratio_sum, fill = cat_label)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.4,
                 linewidth = 0.5) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "#666666") +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = COL_LOSS) +
    scale_fill_manual(values = c("Single-complete" = COL_SC,
                                 "Split (\u22652 forme)" = COL_SP)) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = sprintf("%s — Ratio DISCO/ORTHO", mot),
      subtitle = "ratio = 1: nessuna perdita  |  linea tratteggiata rossa = 0.5",
      x = NULL, y = "ratio (DISCO sum / ORTHO sum)"
    ) +
    theme_disco(base_size = 10) +
    theme(legend.position = "none")

  p_scatter <- ggplot(df_plot,
      aes(x = ortho_sum, y = disco_sum, color = cat_label)) +
    geom_point(alpha = 0.35, size = 1.5, stroke = 0) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "#444444", linewidth = 0.8) +
    scale_color_manual(values = c("Single-complete" = COL_SC,
                                  "Split (\u22652 forme)" = COL_SP)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      title    = sprintf("%s — DISCO vs ORTHO (sum)", mot),
      subtitle = "Diagonale tratteggiata = parit\u00e0",
      x = "ORTHO sum punteggio",
      y = "DISCO sum punteggio",
      color = NULL
    ) +
    theme_disco(base_size = 10) +
    theme(legend.position = "bottom")

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_ratio, newpage = FALSE)
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p_scatter, newpage = FALSE)
  popViewport()
}

dev.off()
cat(sprintf("  Salvato: %s\n", motif_pdf))

cat(sprintf("\n=== Done. Output in: %s ===\n", opts$outdir))