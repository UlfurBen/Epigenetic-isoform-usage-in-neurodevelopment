###############################################################################
# COMBINED SCRIPT: Lollipop Plot of Missense Variants with InterPro Domains
###############################################################################

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

###############################################################################
# 1) Load and parse ClinVar file
###############################################################################

clinvar_raw <- read_tsv("clinvar_result_KMT2A.txt", show_col_types = FALSE)

clinvar_for_plot <- clinvar_raw %>%
  filter(`Molecular consequence` == "missense variant") %>%
  mutate(
    protein_annotation = str_extract(Name, "\\(p\\.[^\\)]+\\)"),
    PROTEIN_POS = str_extract_all(protein_annotation, "[0-9]+")
  ) %>%
  unnest(PROTEIN_POS) %>%
  mutate(PROTEIN_POS = as.numeric(PROTEIN_POS)) %>%
  filter(!is.na(PROTEIN_POS))

# Count variants per amino acid position
clinvar_counts <- clinvar_for_plot %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n(), .groups = "drop")

###############################################################################
# 2) Load and process InterPro domain annotation
###############################################################################

domains_raw <- read_tsv("entry-matching-Q03164.tsv", show_col_types = FALSE)

domain_data <- domains_raw %>%
  filter(!is.na(Matches)) %>%
  separate_rows(Matches, sep = ",") %>%
  mutate(
    Match_Pos = str_extract(Matches, "[0-9]+\\.\\.[0-9]+"),
    Start = as.numeric(str_extract(Match_Pos, "^[0-9]+")),
    End = as.numeric(str_extract(Match_Pos, "[0-9]+$"))
  ) %>%
  select(Name, Type, Start, End) %>%
  arrange(Start)

###############################################################################
# 3) Plot lollipop plot
###############################################################################

# Set height for domain rectangles below x-axis
domain_height <- max(clinvar_counts$Count, na.rm = TRUE) * 0.2

# Custom theme
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Combine everything in one plot
plot_lollipop <- ggplot() +
  
  # Domain boxes
  geom_rect(
    data = domain_data,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0),
    fill = "gray70",
    color = "black",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  
  # Lollipop stems
  geom_segment(
    data = clinvar_counts,
    aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
    color = "firebrick"
  ) +
  
  # Lollipop heads
  geom_point(
    data = clinvar_counts,
    aes(x = PROTEIN_POS, y = Count),
    color = "firebrick",
    size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous("Number of Variants", limits = c(-domain_height, max(clinvar_counts$Count) + 1), expand = c(0, 0)) +
  ggtitle("ClinVar Missense Variants in KMT2A") +
  custom_theme

###############################################################################
# 4) Save plot
###############################################################################

ggsave("lollipop_ClinVar_kmt2a.png", plot_lollipop, width = 10, height = 4, dpi = 300)
message("✅ Lollipop plot with domains saved as: lollipop_ClinVar_kmt2a.png")
