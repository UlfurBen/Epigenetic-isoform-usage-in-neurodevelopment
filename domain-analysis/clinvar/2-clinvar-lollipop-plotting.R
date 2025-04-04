# This code inputs the clinvar plot preprocessed file and plots it in a lollipop plot. 

###############################################################################
# LOLLIPOP PLOT FOR ClinVar Missense Variants in MECP2 (Domains Below)
###############################################################################
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ✅ 1) Load domain data
mecp2_domains <- read_tsv("mecp2_domain_variant_counts.tsv", show_col_types = FALSE)

# ✅ 2) Load ClinVar variant data
clinvar_raw <- read_tsv("clinvar_result_MECP2.txt", show_col_types = FALSE)

# ✅ 3) Extract protein positions
clinvar_for_plot <- clinvar_raw %>%
  filter(`Molecular consequence` == "missense variant") %>%
  filter(str_detect(`Protein change`, "p\\.")) %>%
  mutate(PROTEIN_POS = as.numeric(str_extract(`Protein change`, "[0-9]+"))) %>%
  filter(!is.na(PROTEIN_POS))

# ✅ 4) Count variants per amino acid position
clinvar_counts <- clinvar_for_plot %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n()) %>%
  ungroup()

# ✅ 5) Determine height for domain boxes
domain_height <- max(clinvar_counts$Count, na.rm = TRUE) * 0.2

# ✅ 6) Custom ggplot theme
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# ✅ 7) Create the lollipop plot
plot_clinvar_mecp2 <- ggplot() +
  
  # Add domain boxes below the lollipops
  geom_rect(
    data = mecp2_domains,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0, fill = Name),
    alpha = 0.4,
    color = "black"
  ) +
  
  # Add lollipop stems
  geom_segment(
    data = clinvar_counts,
    aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
    color = "firebrick"
  ) +
  
  # Add lollipop heads
  geom_point(
    data = clinvar_counts,
    aes(x = PROTEIN_POS, y = Count),
    color = "firebrick",
    size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous("Number of Variants", limits = c(-domain_height, max(clinvar_counts$Count, na.rm = TRUE) + 1)) +
  ggtitle("ClinVar Missense Variants in MECP2") +
  custom_theme

# ✅ 8) Save the plot
ggsave("lollipop_ClinVar_mecp2.png", plot_clinvar_mecp2, width = 10, height = 4, dpi = 300)
message("✅ ClinVar lollipop plot for MECP2 saved successfully!")
