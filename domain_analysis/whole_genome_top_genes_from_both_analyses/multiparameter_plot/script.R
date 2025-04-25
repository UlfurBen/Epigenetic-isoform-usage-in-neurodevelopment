library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)

# Load the ranked data
ranked <- read_csv("intersected_genes_ranked_by_composite_score.csv", show_col_types = FALSE)

# Compute -log10(min FDR) for plotting
ranked <- ranked %>%
  mutate(
    min_fdr = pmin(isoform_fdr, exon_fdr),
    neg_log10_min_fdr = -log10(min_fdr + 1e-10),
    is_infinite = is.infinite(odds_ratio)  # Flag infinite odds ratios
  )

# Get top 10 genes by composite score
top10 <- ranked %>% slice_max(order_by = composite_score, n = 10)

# Plot
ggplot(ranked, aes(x = odds_ratio, y = neg_log10_min_fdr, color = is_infinite)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_manual(
    values = c("FALSE" = "purple3", "TRUE" = "red3"),
    labels = c("FALSE" = "Finite Odds Ratio", "TRUE" = "Infinite Odds Ratio")
  ) +
  labs(
    title = "Composite scoring of exon + isoform significant genes",
    subtitle = "Odds ratio (x-axis), -log10(min FDR) (y-axis), color shows infinite vs finite OR",
    x = "Odds Ratio (exon level)",
    y = "-log10(FDR) (min of isoform & exon)",
    color = "Odds Ratio Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_text_repel(data = top10, aes(label = gene_name), size = 4.5, color = "black", box.padding = 0.5)

# Save the plot
ggsave("multiparameter_gene_scoring_plot_inf_highlighted.pdf", width = 10, height = 6, dpi = 300)
