library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)

# 1. Load ranked EM gene data
ranked <- read_csv("EM_genes_non_overlapping_intersected_genes_ranked_by_composite_score.csv", show_col_types = FALSE)

# 2. Compute plotting metrics
ranked <- ranked %>%
  mutate(
    is_infinite = is.infinite(odds_ratio),
    odds_ratio_capped = ifelse(is_infinite, 1024, odds_ratio),
    log2_odds_ratio = pmin(log2(odds_ratio_capped + 1e-6), 10),
    neg_log10_min_fdr = -log10(pmin(isoform_fdr, exon_fdr) + 1e-10)
  )

# 3. Get top 10 genes by composite score for labeling
top10 <- ranked %>% slice_max(order_by = composite_score, n = 10)

# 4. Base plot (color = composite_score, shape = is_infinite)
em_plot <- ggplot(ranked, aes(x = log2_odds_ratio, y = neg_log10_min_fdr,
                              color = composite_score, shape = is_infinite)) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_gradient(low = "gray80", high = "darkred") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Finite", "TRUE" = "Infinite")
  ) +
  labs(
    title = "Composite scoring of EM genes",
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

# 5. Hide shape legend if all entries are finite
if (!any(ranked$is_infinite)) {
  em_plot <- em_plot + guides(shape = "none")
}

# 6. Save and show plot
ggsave("em_genes_non_overlapping_composite_plot_shapes_conditional_legend.pdf", plot = em_plot, width = 10, height = 6, dpi = 300)
print(em_plot)
