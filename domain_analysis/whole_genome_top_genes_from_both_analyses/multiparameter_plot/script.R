library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)

# 1. Load the ranked data
ranked <- read_csv("intersected_non_overlapping_genes_ranked_by_composite_score.csv", show_col_types = FALSE)

# 2. Compute -log10(min FDR), log2(odds ratio), and flag infinite ORs
ranked <- ranked %>%
  mutate(
    is_infinite = is.infinite(odds_ratio),
    odds_ratio_capped = ifelse(is_infinite, 1024, odds_ratio),
    log2_odds_ratio = pmin(log2(odds_ratio_capped + 1e-6), 10),
    min_fdr = pmin(isoform_fdr, exon_fdr),
    neg_log10_min_fdr = -log10(min_fdr + 1e-10)
  )

# 3. Get top 10 by composite score
top10 <- ranked %>% slice_max(order_by = composite_score, n = 10)

# 4. Create base plot (color = composite_score, shape = infinite flag)
whole_plot <- ggplot(ranked, aes(x = log2_odds_ratio, y = neg_log10_min_fdr,
                                 color = composite_score, shape = is_infinite)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_gradient(low = "gray80", high = "purple3") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Finite", "TRUE" = "Infinite")
  ) +
  labs(
    title = "Composite scoring of exon + isoform significant genes",
    subtitle = "Color = composite score, Shape = infinite vs finite odds ratio",
    x = "log2(Odds Ratio)",
    y = "-log10(FDR) (min of isoform & exon)",
    color = "Composite Score",
    shape = "Odds Ratio Type"
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
  geom_text_repel(data = top10, aes(label = gene_name),
                  size = 4.5, color = "black", box.padding = 0.5)

# 5. Remove shape legend if all ORs are finite
if (!any(ranked$is_infinite)) {
  whole_plot <- whole_plot + guides(shape = "none")
}

# 6. Save and show plot
ggsave("multiparameter_non_overlapping_gene_scoring_plot_with_shapes.pdf", plot = whole_plot, width = 10, height = 6, dpi = 300)
print(whole_plot)
