library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# 1. Load EM gene isoform expression data
expr_data <- read.csv("isoform_expression_levels_with_npc.csv")  # Ensure this is the EM genes expression file

# 2. Convert to long format
long_df <- expr_data %>%
  pivot_longer(cols = -c(isoform_id, gene_name), names_to = "sample", values_to = "expression") %>%
  mutate(
    day = case_when(
      grepl("NPC", sample) ~ "NPC",
      grepl("Day3", sample) ~ "Day3",
      grepl("Day6", sample) ~ "Day6",
      grepl("Day12", sample) ~ "Day12",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(day), !is.na(expression))

# 3. Function to create volcano data for a pairwise comparison
pairwise_volcano_data <- function(df, dayA, dayB, comparison_name) {
  df %>%
    filter(day %in% c(dayA, dayB)) %>%
    group_by(isoform_id, gene_name) %>%
    summarise(
      meanA = mean(expression[day == dayA]),
      meanB = mean(expression[day == dayB]),
      log2FC = log2((meanA + 1e-6) / (meanB + 1e-6)),
      p_value = tryCatch(t.test(expression[day == dayA], expression[day == dayB])$p.value, error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(
      fdr = p.adjust(p_value, method = "fdr"),
      comparison = comparison_name
    )
}

# 4. Generate volcano data for all four comparisons
volcano_npc_vs_day3  <- pairwise_volcano_data(long_df, "NPC", "Day3", "NPC vs Day3")
volcano_day3_vs_day6 <- pairwise_volcano_data(long_df, "Day3", "Day6", "Day3 vs Day6")
volcano_day6_vs_day12 <- pairwise_volcano_data(long_df, "Day6", "Day12", "Day6 vs Day12")
volcano_day3_vs_day12 <- pairwise_volcano_data(long_df, "Day3", "Day12", "Day3 vs Day12")

# 5. Combine all into one dataframe
combined_volcano_df <- bind_rows(
  volcano_npc_vs_day3,
  volcano_day3_vs_day6,
  volcano_day6_vs_day12,
  volcano_day3_vs_day12
)

# 6. Highlight top 10 hits per comparison
top_labels <- combined_volcano_df %>%
  filter(!is.na(fdr)) %>%
  group_by(comparison) %>%
  slice_min(fdr, n = 10) %>%
  ungroup()

# 7. Plot
volcano_plot <- ggplot(combined_volcano_df, aes(x = log2FC, y = -log10(fdr))) +
  geom_point(alpha = 0.5, color = "#4169E1") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.25,
    point.padding = 0.2
  ) +
  facet_wrap(~comparison, scales = "free") +
  theme_minimal(base_size = 13) +
  labs(
    title = "EM Genes Faceted Volcano Plot of Isoform Expression",
    x = "log2(Fold Change)",
    y = "-log10(FDR)"
  )

# 8. Save plot
ggsave("em_genes_faceted_isoform_volcano_plot_4comparisons.png", volcano_plot, width = 14, height = 8, dpi = 300)
print(volcano_plot)
