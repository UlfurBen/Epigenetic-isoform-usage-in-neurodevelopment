###############################################################################
# LOCAL PLOTTING CODE
###############################################################################
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# 1) READ CSV FILES
mecp2_domains      <- read.csv("mecp2_domains.csv", stringsAsFactors = FALSE)
clinvar_for_plot   <- read.csv("plotA_clinvar_mecp2.csv", stringsAsFactors = FALSE)
gnomad_for_plot    <- read.csv("plotB_gnomad_mecp2.csv", stringsAsFactors = FALSE)

# ✅ Extract protein positions from `hgvsp` in gnomAD data
gnomad_for_plot <- gnomad_for_plot %>%
  mutate(PROTEIN_POS = as.numeric(str_extract(hgvsp, "[0-9]+")))

# 2) EXAMPLE: PLOT A (ClinVar)
plotA <- ggplot() +
  geom_rect(
    data = mecp2_domains,
    aes(
      xmin = Start,
      xmax = End,
      ymin = 0,
      ymax = 1,
      fill = Name  
    ),
    alpha = 0.5,  
    color = "black"
  ) +
  geom_point(
    data = clinvar_for_plot,
    aes(
      x = PROTEIN_POS,   
      y = 0.5
    ),
    size = 2
  ) +
  scale_x_continuous("Protein position (amino acid)") +
  scale_y_continuous("", breaks = NULL) +
  ggtitle("Plot A: ClinVar Missense Variants in MECP2") +
  theme_bw()

# Save to file
ggsave("plotA_ClinVar_MECP2.png", plotA, width = 8, height = 3, dpi = 300)

# 3) EXAMPLE: PLOT B (gnomAD)
plotB <- ggplot() +
  geom_rect(
    data = mecp2_domains,
    aes(
      xmin = Start,
      xmax = End,
      ymin = 0,
      ymax = 1,
      fill = Name
    ),
    alpha = 0.5,
    color = "black"
  ) +
  geom_point(
    data = gnomad_for_plot,
    aes(
      x = PROTEIN_POS,  
      y = 0.5
    ),
    size = 2
  ) +
  scale_x_continuous("Protein position (amino acid)") +
  scale_y_continuous("", breaks = NULL) +
  ggtitle("Plot B: gnomAD Missense Variants in MECP2") +
  theme_bw()

ggsave("plotB_gnomAD_MECP2.png", plotB, width = 8, height = 3, dpi = 300)

message("✅ Plots created locally!")
