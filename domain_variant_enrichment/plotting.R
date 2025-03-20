###############################################################################
# LOCAL PLOTTING CODE
###############################################################################
library(ggplot2)
library(dplyr)
library(readr)

# 1) READ CSV FILES
mecp2_domains      <- read.csv("mecp2_domains.csv", stringsAsFactors = FALSE)
clinvar_for_plot   <- read.csv("plotA_clinvar_mecp2.csv", stringsAsFactors = FALSE)
gnomad_for_plot    <- read.csv("plotB_gnomad_mecp2.csv", stringsAsFactors = FALSE)

# 2) EXAMPLE: PLOT A (ClinVar)
#    We assume you want an x-axis from 1..~486 (length of MECP2 or so).
#    We'll draw each domain as a rectangle and place points for each variant
#    by protein_position. If you do NOT have actual protein coordinates,
#    you'd need to use another approach or omit the points.

plotA <- ggplot() +
  # Draw domain rectangles (each row in mecp2_domains)
  geom_rect(
    data = mecp2_domains,
    aes(
      xmin = Start,
      xmax = End,
      ymin = 0,
      ymax = 1,   # a simple height for the domain bars
      fill = Name
    ),
    alpha = 0.5,  # So they do not fully cover each other if they overlap
    color = "black"
  ) +
  # Plot each ClinVar variant as a point on top
  geom_point(
    data = clinvar_for_plot,
    aes(
      x = protein_position,   # numeric position along x-axis
      y = 0.5                 # put them in the middle of the domain bar
    ),
    # color, shape, etc. can be adjusted as you wish
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
      x = protein_position,
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
