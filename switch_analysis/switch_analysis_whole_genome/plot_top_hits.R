# 1) Load Required Libraries -------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 2) Read the Isoform Expression Data -----------------------------------
isoform_data <- read.csv("~/downloads/top_hits_whole_genome_significant_isoform_upregulation_by_day.csv", header = TRUE)

# 3) Define Top Significant Isoforms -----------------------------------
top_isoforms <- c("ENSMUST00000037206", "ENSMUST00000043739", "ENSMUST00000097952", "ENSMUST00000121983", "ENSMUST00000123283",
                 "ENSMUST00000125105", "ENSMUST00000127104", "ENSMUST00000129371", "ENSMUST00000130899", "ENSMUST00000136971")

# 4) Filter and Aggregate Data -----------------------------------------
plot_data <- isoform_data %>%
  filter(isoform_id %in% top_isoforms) %>%
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

# 5) Plot All in One Page (Max 10 Isoforms) ----------------------------
pdf("whole_genome_top_significant_isoform_expression_plots.pdf", width = 12, height = 8)

p <- ggplot(plot_data, aes(x = Day, y = Mean_Expression, fill = Day)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Top Significantly Upregulated Isoforms",
       x = "Condition",
       y = "Mean Expression Level") +
  theme_minimal() +
  scale_fill_manual(values = c("Day3" = "#E69F00", "Day6" = "#56B4E9", "Day12" = "#009E73")) +
  facet_wrap(~ paste(gene_name, isoform_id, sep = " - "), ncol = 3, nrow = 4, scales = "free_y") +
  theme(strip.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()

# Print Completion Message ----------------------------------------------
print("Top significant isoform plots saved to whole_genome_top_significant_isoform_expression_plots.pdf")
