# Load required libraries
library(ggplot2)
library(dplyr)

# Read the file you transferred
isoform_data <- read.csv("~/downloads/isoform_usage_results.csv", header = TRUE)

# Check column structure
str(isoform_data)

# Filter for the isoform of interest
iso_of_interest <- "ENSMUST00000023741"
isoform_subset <- isoform_data %>%
  filter(isoform_id == iso_of_interest)

# Check if isoform data exists
if (nrow(isoform_subset) == 0) {
  stop("No data found for the specified isoform!")
}

# Reshape Data for ggplot2
isoform_long <- isoform_subset %>%
  select(condition_1, condition_2, IF1, IF2) %>%
  tidyr::pivot_longer(cols = c(IF1, IF2), names_to = "IF_type", values_to = "IF_value")

# Modify labels for clarity
isoform_long$Condition_Comparison <- paste(isoform_long$condition_1, "vs", isoform_long$condition_2)

# Plot isoform fraction across conditions
ggplot(isoform_long, aes(x = Condition_Comparison, y = IF_value, fill = IF_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = paste("Isoform Usage of", iso_of_interest),
       x = "Condition Comparison",
       y = "Isoform Fraction (IF)",
       fill = "Isoform Usage") +
  theme_minimal() +
  scale_fill_manual(values = c("IF1" = "blue", "IF2" = "red"))

# Save as PDF
ggsave("isoform_usage_plot.pdf", width = 8, height = 6)
