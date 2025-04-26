###############################################################
# Script 2: Plotting Isoform Expression (Top Significant, EM Genes, Correct NPC Split)
###############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Load expression and significant isoform data
expr_df <- read.csv("isoform_expression_levels_with_npc.csv")
sig_isoforms <- read.csv("isoform_significant_with_npc_fdr_fc.csv")

# Reshape expression data
expr_df_long <- expr_df %>%
  pivot_longer(cols = -c(isoform_id, gene_name), names_to = "sample", values_to = "expression") %>%
  mutate(Day = case_when(
    grepl("NPC", sample) ~ "NPC",
    grepl("Day3", sample) ~ "Day3",
    grepl("Day6", sample) ~ "Day6",
    grepl("Day12", sample) ~ "Day12"
  ))

# Calculate maximum expression per isoform per day
max_expr_per_isoform <- expr_df_long %>%
  group_by(isoform_id, Day) %>%
  summarise(max_expr = max(expression, na.rm = TRUE), .groups = "drop")

# Find isoforms where NPC expression is highest
npc_high_isoforms <- max_expr_per_isoform %>%
  group_by(isoform_id) %>%
  filter(max_expr == max(max_expr)) %>%
  filter(Day == "NPC") %>%
  pull(isoform_id)

# Find isoforms where NPC is NOT the highest
non_npc_high_isoforms <- max_expr_per_isoform %>%
  group_by(isoform_id) %>%
  filter(max_expr == max(max_expr)) %>%
  filter(Day != "NPC") %>%
  pull(isoform_id)

# Top 10 by FDR where NPC is highest
top_isoforms_npc_high <- sig_isoforms %>%
  arrange(fdr) %>%
  filter(isoform_id %in% npc_high_isoforms) %>%
  slice_head(n = 10) %>%
  pull(isoform_id)

# Top 10 by FDR where NPC is NOT highest
top_isoforms_not_npc_high <- sig_isoforms %>%
  arrange(fdr) %>%
  filter(isoform_id %in% non_npc_high_isoforms) %>%
  slice_head(n = 10) %>%
  pull(isoform_id)

# Function to prepare plot data
prepare_plot_data <- function(selected_isoforms) {
  plot_data <- expr_df_long %>%
    filter(isoform_id %in% selected_isoforms) %>%
    group_by(isoform_id, gene_name, Day) %>%
    summarise(Mean_Expression = mean(expression), .groups = "drop") %>%
    left_join(sig_isoforms %>% dplyr::select(isoform_id, fdr), by = "isoform_id")
  
  plot_data$Day <- factor(plot_data$Day, levels = c("NPC", "Day3", "Day6", "Day12"))
  return(plot_data)
}

# Prepare data for both plots
plot_data_npc_high <- prepare_plot_data(top_isoforms_npc_high)
plot_data_not_npc_high <- prepare_plot_data(top_isoforms_not_npc_high)

# Function to create plot with limited height per isoform
create_plot <- function(plot_data, title_text) {
  
  facet_y_limits <- plot_data %>%
    group_by(isoform_id) %>%
    summarise(max_y = max(Mean_Expression, na.rm = TRUE)) %>%
    mutate(y_limit = max_y * 1.5)  # Leave 50% headroom
  
  plot_data <- plot_data %>%
    left_join(facet_y_limits, by = "isoform_id")
  
  label_data <- plot_data %>%
    group_by(isoform_id, gene_name, fdr, y_limit) %>%
    summarise(Mean_Expression = max(Mean_Expression), .groups = "drop") %>%
    mutate(label = sprintf("FDR = %.1e", fdr))
  
  ggplot(plot_data, aes(x = Day, y = Mean_Expression, fill = Day)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = title_text,
         x = "Condition", y = "Mean Expression Level") +
    theme_minimal() +
    scale_fill_manual(values = c("NPC" = "#999999", "Day3" = "#E69F00", "Day6" = "#56B4E9", "Day12" = "#009E73")) +
    facet_wrap(~ paste(gene_name, isoform_id, sep = " - "), ncol = 2, scales = "free_y") +
    geom_text(data = label_data, aes(x = 2.5, y = y_limit * 0.95, label = label),
              inherit.aes = FALSE, size = 2.5, vjust = 0, hjust = 0.5) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create both plots
p1 <- create_plot(plot_data_npc_high, "Top 10 Isoforms (NPC Highest Expression, EM Genes)")
p2 <- create_plot(plot_data_not_npc_high, "Top 10 Isoforms (Day3/6/12 Highest Expression, EM Genes)")

# Save to PDF side-by-side
pdf("top_significant_isoform_expression_side_by_side_em_genes_correct_npc.pdf", width = 16, height = 10)
(p1 | p2)  # patchwork side-by-side
dev.off()
