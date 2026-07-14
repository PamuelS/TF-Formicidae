library(ggplot2)
library(ggrepel)
library(viridis) # for the color scale of the dots

dataset_filtered <- read.table("pgls5_significant_OGs_pvalue_only.tsv", header = TRUE, sep = " ")

# letting know the program that the column are numeric
dataset_filtered$P_value <- as.numeric(as.character(dataset_filtered$P_value))
dataset_filtered$R_squared_adj <- as.numeric(as.character(dataset_filtered$R_squared_adj))

dataset_filtered <- dataset_filtered[!is.na(as.character(dataset_filtered$P_value)) & !is.na(as.character(dataset_filtered$R_squared_adj)), ]

# count of the OGs with both parameter setted
num_OGs <- sum(dataset_filtered$R_squared_adj > 0.25)
cat("Found ", num_OGs, "OGs with R^2 > 0.25\n")

# Actual rappresentation of the data
ggplot(dataset_filtered,
       aes(R_squared_adj, P_value)) +

    geom_point(
        color = "grey85",
        alpha = 0.05,
        size = 0.5
    ) +

    geom_point(
        data = subset(dataset_filtered,
                      R_squared_adj > 0.25),
        aes(color = R_squared_adj),
        alpha = 0.8,
        size = 1.5
    ) +

    scale_color_viridis_c(
    option = "plasma",
    name = expression("Adjusted " * R^2),
    breaks = c(0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95),
    limits = c(0.25, max(dataset_filtered$R_squared_adj))
    ) +

    geom_vline(
        xintercept = 0.25,
        linetype = "dashed"
    ) +
    labs(
        title = "PGLS-castes signnificant OG",
        x = expression("Adjusted " * R^2),
        y = expression(P[value])
    )

    ggsave("pgls5_significant_OGs_pvalue_only.png", 
        width = 10, 
        height = 8, 
        dpi = 300
    )
    theme_minimal()
