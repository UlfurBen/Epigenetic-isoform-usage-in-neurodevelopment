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
    neg_log10_min_fdr = -log10(min_fdr + 1e-10)
  )

# Get top 10 by composite score
top10 <- ranked %>% slice_max(order_by = composite_score, n = 10)

# Plot
ggplot(ranked, aes(x = odds_ratio, y = neg_log10_min_fdr, color = composite_score)) +
  geom_point(size = 3.5, alpha = 0.8) +
  scale_color_gradient(low = "gray80", high = "purple3") +
  labs(
    title = "Composite scoring of exon + isoform significant genes",
    subtitle = "Odds ratio (x-axis), -log10(min FDR) (y-axis), composite score (color)",
    x = "Odds Ratio (exon level)",
    y = "-log10(FDR) (min of isoform & exon)",
    color = "Composite Score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),  # <- updated here
    axis.ticks = element_line(color = "black", linewidth = 0.8), # <- updated here
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_text_repel(data = top10, aes(label = gene_name), size = 4.5, color = "black", box.padding = 0.5)

  ggsave("multiparameter_gene_scoring_plot.pdf", width = 10, height = 6, dpi = 300)
