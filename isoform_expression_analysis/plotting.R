# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the file containing isoform usage results
isoform_data <- read.csv("~/downloads/isoform_usage_results.csv", header = TRUE)

# Define isoforms of interest
isoforms_of_interest <- c("ENSMUST00000033770", "ENSMUST00000169922")

# Filter data for the selected isoforms
isoform_subset <- isoform_data %>%
  filter(isoform_id %in% isoforms_of_interest)

# Reshape the data to combine all comparisons into a single plot
plot_data <- isoform_subset %>%
  pivot_longer(cols = c(IF1, IF2), names_to = "IF_type", values_to = "Isoform_Fraction") %>%
  mutate(Day = ifelse(IF_type == "IF1", condition_1, condition_2)) %>%
  select(isoform_id, Day, Isoform_Fraction)

# Ensure Days are ordered correctly
plot_data$Day <- factor(plot_data$Day, levels = c("Day3", "Day6", "Day12"))

# Create the bar plot for each isoform
for (isoform in isoforms_of_interest) {
  p <- ggplot(plot_data %>% filter(isoform_id == isoform), aes(x = Day, y = Isoform_Fraction, fill = Day)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("Isoform Usage for", isoform),
         x = "Day",
         y = "Isoform Fraction (IF)") +
    theme_minimal() +
    scale_fill_manual(values = c("Day3" = "#E69F00", "Day6" = "#56B4E9", "Day12" = "#009E73"))
  
  # Save the plot as a PDF
  ggsave(filename = paste0("isoform_usage_", isoform, ".pdf"), plot = p, width = 8, height = 6)
}

# Display message
print("Plots saved for each isoform with all days together.")

