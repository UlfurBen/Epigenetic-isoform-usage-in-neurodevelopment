library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Load data
expr_df <- read.csv("whole_genome_isoform_expression_levels_with_npc.csv")
sig_isoforms <- read.csv("whole_genome_isoform_significant_with_npc_fdr_fc.csv")

# Pick top 10 by FDR
top_isoforms <- sig_isoforms %>% arrange(fdr) %>% head(10) %>% pull(isoform_id)

# Prepare data for plotting
plot_data <- expr_df %>%
  filter(isoform_id %in% top_isoforms) %>%
  pivot_longer(cols = -c(isoform_id, gene_name), names_to = "sample", values_to = "expression") %>%
  mutate(Day = case_when(
    grepl("NPC", sample) ~ "NPC",
    grepl("Day3", sample) ~ "Day3",
    grepl("Day6", sample) ~ "Day6",
    grepl("Day12", sample) ~ "Day12"
  )) %>%
  group_by(isoform_id, gene_name, Day) %>%
  summarise(Mean_Expression = mean(expression), .groups = "drop")

plot_data$Day <- factor(plot_data$Day, levels = c("NPC", "Day3", "Day6", "Day12"))

# Plot to PDF
pdf("top_significant_isoform_expression_plots_with_npc_whole_genome.pdf", width = 12, height = 8)
p <- ggplot(plot_data, aes(x = Day, y = Mean_Expression, fill = Day)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Top Significantly Upregulated Isoforms (Whole Genome, incl. NPC)",
       x = "Condition", y = "Mean Expression Level") +
  theme_minimal() +
  scale_fill_manual(values = c("NPC" = "#999999", "Day3" = "#E69F00", "Day6" = "#56B4E9", "Day12" = "#009E73")) +
  facet_wrap(~ paste(gene_name, isoform_id, sep = " - "), ncol = 3, scales = "free_y") +
  theme(strip.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
