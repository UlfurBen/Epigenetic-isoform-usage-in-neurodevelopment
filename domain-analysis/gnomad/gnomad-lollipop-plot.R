# This code creates a lollipop graph using the data from the output file from gnomad-plot-preprocessing.R

###############################################################################
# LOLLIPOP PLOT FOR gnomAD Missense Variants (Domains Below Lollipops)
###############################################################################
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ✅ Load domain data (with Start, End, Name, and variant counts)
chd3_domains <- read_csv("chd3_domain_variant_counts.csv", show_col_types = FALSE)

# ✅ Load gnomAD variant data (from the modified script)
gnomad_for_plot <- read_csv("chd3_variants.csv", show_col_types = FALSE)

# ✅ Extract protein positions from `hgvsp`
gnomad_for_plot <- gnomad_for_plot %>%
  mutate(PROTEIN_POS = as.numeric(str_extract(hgvsp, "[0-9]+"))) %>%
  filter(!is.na(PROTEIN_POS))

# ✅ Count variants per amino acid position
gnomad_counts <- gnomad_for_plot %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n()) %>%
  ungroup()

# ✅ Determine height of domain boxes
domain_height <- max(gnomad_counts$Count, na.rm = TRUE) * 0.2

# ✅ Custom theme for plot
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# ✅ Plot gnomAD variants with domains below
plot_gnomad <- ggplot() +
  
  # Colored boxes for domains
  geom_rect(
    data = chd3_domains,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0, fill = Name),
    alpha = 0.4,
    color = "black"
  ) +
  
  # Lollipop stems
  geom_segment(
    data = gnomad_counts,
    aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
    color = "black"
  ) +
  
  # Lollipop heads
  geom_point(
    data = gnomad_counts,
    aes(x = PROTEIN_POS, y = Count),
    color = "black",
    size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous("Number of Variants", limits = c(-domain_height, max(gnomad_counts$Count, na.rm = TRUE) + 1)) +
  ggtitle("gnomAD Missense Variants in CHD3") +
  custom_theme

# ✅ Save the figure
ggsave("lollipop_gnomAD_chd3.png", plot_gnomad, width = 10, height = 4, dpi = 300)
message("✅ gnomAD lollipop plot saved successfully!")

