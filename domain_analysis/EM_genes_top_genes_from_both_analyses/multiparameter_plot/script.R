library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)

# 1. Load ranked EM gene data
ranked <- read_csv("EM_genes_intersected_genes_ranked_by_composite_score.csv", show_col_types = FALSE)

# 2. Handle infinite odds ratios and compute log-transformed metrics
ranked <- ranked %>%
  mutate(
    is_infinite = is.infinite(odds_ratio),  # Flag infinite values
    odds_ratio_capped = ifelse(is_infinite, 1024, odds_ratio),  # Cap for log scale
    log2_odds_ratio = pmin(log2(odds_ratio_capped + 1e-6), 10),  # Cap at log2(1024)
    neg_log10_min_fdr = -log10(pmin(isoform_fdr, exon_fdr) + 1e-10)
  )

# 3. Select top 10 by composite score for labeling
top10 <- ranked %>% slice_max(order_by = composite_score, n = 10)

# 4. Plot
em_plot <- ggplot(ranked, aes(x = log2_odds_ratio, y = neg_log10_min_fdr, color = is_infinite)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(
    values = c("FALSE" = "darkred", "TRUE" = "dodgerblue3"),
    labels = c("FALSE" = "Finite Odds Ratio", "TRUE" = "Infinite Odds Ratio")
  ) +
  labs(
    title = "Composite scoring of EM genes",
    subtitle = "log2(odds ratio), -log10(min FDR), color highlights infinite OR",
    x = "log2(Odds Ratio)",
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

# 5. Save plot
ggsave("em_genes_composite_plot_inf_highlighted.pdf", plot = em_plot, width = 10, height = 6, dpi = 300)

# Optional: Show it
print(em_plot)
