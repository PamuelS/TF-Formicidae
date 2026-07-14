library(tidyverse)

# File reading
data_filtered <- read.table("pgls5_significant_OGs_pvalue_only.tsv", 
                             header = TRUE, sep = " ")

# Counting OGs and their direction
df_plot <- data_filtered %>%
  group_by(motivo) %>%
  summarise(
    polymorphic_bias = sum(Regression_Coefficient > 0),
    monomorphic_bias = -sum(Regression_Coefficient < 0),
    .groups = "drop"
  ) %>%
  arrange(desc(polymorphic_bias)) %>%
  mutate(motivo_order = factor(motivo, levels = motivo))

# Creation of the dataframe for the plot
df_long <- df_plot %>%
  pivot_longer(
    cols      = c(polymorphic_bias, monomorphic_bias),
    names_to  = "direction",
    values_to = "n_OG"
  ) %>%
  mutate(direction = factor(direction,
                            levels = c("polymorphic_bias", "monomorphic_bias")))

# Part of the plot
ggplot(df_long, aes(x = motivo_order, y = n_OG, fill = direction)) +
  geom_col(width = 0.8) +
  
  # Line 0
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  
  # Coloring of the bars
  scale_fill_manual(
    values = c("polymorphic_bias" = "#47218d",
               "monomorphic_bias" = "#dfc21d"),
    labels = c("Polymorphic Bias", "Monomorphic Bias")
  ) +
  
  # Extention of the Y-axe
  scale_y_continuous(
    breaks = pretty(c(min(df_long$n_OG), max(df_long$n_OG)), n = 6),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  # Instruction about the axes
  labs(x = "Motifs", y = "Number of OGs", fill = NULL) +
  
  # Words can go over the edge of the plot
  coord_cartesian(clip = "off") +
  
  theme_classic() +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    legend.position    = c(0.85, 0.85),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.margin        = margin(t = 30, r = 10, b = 10, l = 10)
  )

ggsave("motif_bias_plot.pdf", width = 12, height = 5, dpi = 300)
ggsave("motif_bias_plot.png", width = 12, height = 5, dpi = 300)
