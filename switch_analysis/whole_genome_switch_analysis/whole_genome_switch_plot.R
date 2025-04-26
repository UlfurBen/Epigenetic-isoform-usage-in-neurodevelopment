library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # for side-by-side plots!

# Load data
expr_df <- read.csv("whole_genome_isoform_expression_levels_with_npc.csv")
sig_isoforms <- read.csv("whole_genome_isoform_significant_with_npc_fdr_fc.csv")

# Reshape expression data and extract Day information
expr_df_long <- expr_df %>%
  pivot_longer(cols = starts_with("NPC") | starts_with("Day"), 
               names_to = "sample", values_to = "expression") %>%
  mutate(Day = case_when(
    grepl("NPC", sample) ~ "NPC",
    grepl("Day3", sample) ~ "Day3",
    grepl("Day6", sample) ~ "Day6",
    grepl("Day12", sample) ~ "Day12"
  ))

# Find isoforms with expression at NPC
npc_expression <- expr_df_long %>%
  filter(Day == "NPC") %>%
  group_by(isoform_id) %>%
  summarise(max_expr = max(expression), .groups = "drop") %>%
  filter(max_expr > 0) %>%
  pull(isoform_id)

# Top 10 isoforms including NPC-high
top_isoforms_incl_npc <- sig_isoforms %>%
  arrange(fdr) %>%
  slice_head(n = 10) %>%
  pull(isoform_id)

# Top 10 isoforms excluding NPC-high
top_isoforms_excl_npc <- sig_isoforms %>%
  arrange(fdr) %>%
  filter(!(isoform_id %in% npc_expression)) %>%
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
plot_data_incl_npc <- prepare_plot_data(top_isoforms_incl_npc)
plot_data_excl_npc <- prepare_plot_data(top_isoforms_excl_npc)

# Function to create plot with limited height per facet
create_plot <- function(plot_data, title_text) {
  
  # Calculate maximum y for each isoform to limit bar heights
  facet_y_limits <- plot_data %>%
    group_by(isoform_id) %>%
    summarise(max_y = max(Mean_Expression, na.rm = TRUE)) %>%
    mutate(y_limit = max_y * 1.5) # leave some headroom
  
  # Merge limits into the main plot_data
  plot_data <- plot_data %>%
    left_join(facet_y_limits, by = "isoform_id")
  
  # Only one row per isoform for FDR labeling
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
p1 <- create_plot(plot_data_incl_npc, "Top 10 Significant Isoforms (Including NPC-high Isoforms)")
p2 <- create_plot(plot_data_excl_npc, "Top 10 Significant Isoforms (Excluding NPC-high Isoforms)")

# Save to PDF side-by-side
pdf("top_significant_isoform_expression_side_by_side.pdf", width = 16, height = 10)
(p1 | p2)  # patchwork magic: side-by-side
dev.off()
