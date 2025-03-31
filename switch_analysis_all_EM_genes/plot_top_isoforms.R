# 1) Load Required Libraries -------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# 2) Read the Isoform Expression Data -----------------------------------
isoform_data <- read.csv("~/downloads/isoform_expression_levels.csv", header = TRUE)

# 3) Define Top Significant Isoforms -----------------------------------
top_isoforms <- c("ENSMUST00000092971", "ENSMUST00000108661", "ENSMUST00000140327",
                  "ENSMUST00000216540", "ENSMUST00000128981", "ENSMUST00000209474",
                  "ENSMUST00000114446", "ENSMUST00000189562", "ENSMUST00000135909",
                  "ENSMUST00000122992")

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
pdf("top_significant_isoform_expression_plots.pdf", width = 12, height = 8)

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
print("Top significant isoform plots saved.")
