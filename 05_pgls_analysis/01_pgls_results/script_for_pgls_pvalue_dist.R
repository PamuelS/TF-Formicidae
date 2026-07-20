library(tidyverse)

# ──────────────────────────────────────────────
# 1. Lettura del file
# ──────────────────────────────────────────────
data_filtered <- read.table("pgls5_significant_OGs_pvalue_only.tsv", 
                             header = TRUE, sep = " ")

# ──────────────────────────────────────────────
# 2. Conta OG per motivo e direzione
# ──────────────────────────────────────────────
df_plot <- data_filtered %>%
  group_by(motivo) %>%
  summarise(
    polymorphic_bias = sum(Regression_Coefficient > 0),
    monomorphic_bias = -sum(Regression_Coefficient < 0),
    .groups = "drop"
  ) %>%
  arrange(desc(polymorphic_bias)) %>%
  mutate(motivo_order = factor(motivo, levels = motivo))

# ──────────────────────────────────────────────
# 3. Formato lungo per ggplot
# ──────────────────────────────────────────────
df_long <- df_plot %>%
  pivot_longer(
    cols      = c(polymorphic_bias, monomorphic_bias),
    names_to  = "direction",
    values_to = "n_OG"
  ) %>%
  mutate(direction = factor(direction,
                            levels = c("polymorphic_bias", "monomorphic_bias")))

# ──────────────────────────────────────────────
# 4. Valori per annotazioni
# ──────────────────────────────────────────────
max_val <- max(df_plot$polymorphic_bias)
n_motifs <- nrow(df_plot)

# Punto di crossover: primo motivo dove |monomorphic| supera polymorphic
crossover <- which(abs(df_plot$monomorphic_bias) > df_plot$polymorphic_bias)[1]
if (is.na(crossover)) crossover <- nrow(df_plot)

# ──────────────────────────────────────────────
# 5. Plot
# ──────────────────────────────────────────────
ggplot(df_long, aes(x = motivo_order, y = n_OG, fill = direction)) +
  geom_col(width = 0.8) +
  
  # Linea dello zero
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +

  scale_fill_manual(
    values = c("polymorphic_bias" = "#47218d",
               "monomorphic_bias" = "#dfc21d"),
    labels = c("Polymorphic Bias", "Monomorphic Bias")
  ) +
  
  # Margine extra in alto per non tagliare l'annotazione
  scale_y_continuous(
    breaks = pretty(c(min(df_long$n_OG), max(df_long$n_OG)), n = 6),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  labs(x = "Motifs", y = "Number of OGs", fill = paste0("Total Motifs: ", n_motifs)) +
  
  # clip = "off" permette al testo di uscire dai bordi del pannello
  coord_cartesian(clip = "off") +
  
  theme_classic() +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    legend.position    = c(0.85, 0.85),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.title       = element_text(face = "bold", size = 10),
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.margin        = margin(t = 30, r = 10, b = 10, l = 10)
  )

ggsave("motif_pvalue_bias_plot.pdf", width = 12, height = 5, dpi = 300)
ggsave("motif_pvalue_bias_plot.png", width = 12, height = 5, dpi = 300)

