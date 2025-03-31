# 1) Load Required Libraries -------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 2) Read the Isoform Expression Data -----------------------------------
isoform_data <- read.csv("~/downloads/isoform_expression_levels.csv", header = TRUE)

# 3) Aggregate Data by Averaging Replicates -----------------------------
plot_data <- isoform_data %>%
  pivot_longer(cols = -c(isoform_id, gene_name), 
               names_to = "Replicate", 
               values_to = "Expression") %>%
  mutate(Day = case_when(
    grepl("Day3_", Replicate) ~ "Day3",
    grepl("Day6_", Replicate) ~ "Day6",
    grepl("Day12_", Replicate) ~ "Day12"
  )) %>%
  group_by(isoform_id, gene_name, Day) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Ensure correct factor ordering
plot_data$Day <- factor(plot_data$Day, levels = c("Day3", "Day6", "Day12"))

# 4) Split Data into Pages of 9 Plots -----------------------------------
unique_isoforms <- unique(plot_data$isoform_id)
isof_pages <- split(unique_isoforms, ceiling(seq_along(unique_isoforms) / 9))

pdf("isoform_expression_plots.pdf", width = 12, height = 8)

for (page in isof_pages) {
  subset_data <- plot_data %>% filter(isoform_id %in% page)
  
  p <- ggplot(subset_data, aes(x = Day, y = Mean_Expression, fill = Day)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Isoform Expression Across Conditions",
         x = "Condition",
         y = "Mean Expression Level") +
    theme_minimal() +
    scale_fill_manual(values = c("Day3" = "#E69F00", "Day6" = "#56B4E9", "Day12" = "#009E73")) +
    facet_wrap(~ paste(gene_name, isoform_id, sep = " - "), ncol = 3, nrow = 3, scales = "free_y") +
    theme(strip.text = element_text(size = 12, face = "bold"), 
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

dev.off()

# Print Completion Message ----------------------------------------------
print("Multi-panel plot saved with 9 plots per page and averaged replicates.")
