# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the file containing isoform usage results
isoform_data <- read.csv("~/downloads/isoform_usage_results.csv", header = TRUE)

# Define isoforms of interest and their corresponding gene names
gene_mapping <- data.frame(
  isoform_id = c("ENSMUST00000033770", "ENSMUST00000023165", "ENSMUST00000176030",
                 "ENSMUST00000087916", "ENSMUST00000113573", "ENSMUST00000099490"),
  gene_name = c("MECP2", "CREBBP", "SMARCA2", "HDAC8", "ATRX", "NSD1")
)

# Filter data for the selected isoforms
isoform_subset <- isoform_data %>%
  filter(isoform_id %in% gene_mapping$isoform_id)

# Merge gene names into the dataset
isoform_subset <- isoform_subset %>%
  left_join(gene_mapping, by = "isoform_id")

# Reshape the data to combine all comparisons into a single plot
plot_data <- isoform_subset %>%
  pivot_longer(cols = c(IF1, IF2), names_to = "IF_type", values_to = "Isoform_Fraction") %>%
  mutate(Day = ifelse(IF_type == "IF1", condition_1, condition_2)) %>%
  select(isoform_id, gene_name, Day, Isoform_Fraction)

# Ensure Days are ordered correctly
plot_data$Day <- factor(plot_data$Day, levels = c("Day3", "Day6", "Day12"))

# Create a multi-panel bar plot for all isoforms, including gene names
p <- ggplot(plot_data, aes(x = Day, y = Isoform_Fraction, fill = Day)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Isoform Usage Across Conditions",
       x = "Day",
       y = "Isoform Fraction (IF)") +
  theme_minimal() +
  scale_fill_manual(values = c("Day3" = "#E69F00", "Day6" = "#56B4E9", "Day12" = "#009E73")) +
  facet_wrap(~ paste(gene_name, isoform_id, sep = " - "), ncol = 2, scales = "free_y") +  # Facet by gene name & isoform
  theme(strip.text = element_text(size = 12, face = "bold"))

# Save the plot as a PDF
ggsave(filename = "~/downloads/isoform_usage_multiple_isoforms.pdf", plot = p, width = 10, height = 6)

# Display message
print("Multi-panel plot saved for all selected isoforms with gene names.")
