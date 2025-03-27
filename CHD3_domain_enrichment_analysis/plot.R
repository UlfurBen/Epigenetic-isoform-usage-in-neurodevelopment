# Ensembl VEP was used to get protein position of ClinVar variants

###############################################################################
# FIXED LOLLIPOP PLOTS FOR ClinVar & gnomAD (Domains Below Lollipops)
###############################################################################
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# 1) FUNCTION TO EXTRACT DOMAIN START & END POSITIONS
get_chd3_domains_from_tsv <- function(file_path = "entry-matching-Q12873-chd3.tsv") {
  # Read the InterPro TSV file
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  # Select relevant columns
  domain_df <- df %>%
    dplyr::select(Accession, Name, `Source Database`, Type, Matches)
  
  # Extract start and end positions from "Matches"
  domain_df <- domain_df %>%
    separate(Matches, into = c("Start", "End"), sep = "\\.\\.", convert = TRUE)
  
  # Convert to numeric
  domain_df <- domain_df %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
    arrange(Start)  # Sort by position
  
  # Save to CSV
  write.csv(domain_df, "chd3_domains.csv", row.names = FALSE)
  message("✅ chd3 domain data saved to: chd3_domains.csv")
  
  return(domain_df)
}

# ✅ GET DOMAIN DATA
chd3_domains <- get_chd3_domains_from_tsv("entry-matching-Q12873-chd3.tsv")

# ✅ READ CLINVAR & GNOMAD VARIANT DATA
clinvar_for_plot   <- read_csv("plotA_clinvar_chd3.csv", show_col_types = FALSE)
gnomad_for_plot    <- read_csv("plotB_gnomad_chd3.csv", show_col_types = FALSE)

# ✅ Extract protein positions from `hgvsp` in gnomAD data
gnomad_for_plot <- gnomad_for_plot %>%
  mutate(PROTEIN_POS = as.numeric(str_extract(hgvsp, "[0-9]+")))

# ✅ Count variants per amino acid position
clinvar_counts <- clinvar_for_plot %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n()) %>%
  ungroup()

gnomad_counts <- gnomad_for_plot %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n()) %>%
  ungroup()

# ✅ Determine domain box height below the x-axis
domain_height <- max(max(clinvar_counts$Count, na.rm = TRUE), max(gnomad_counts$Count, na.rm = TRUE)) * 0.2

# ✅ COMMON THEME FUNCTION FOR CONSISTENCY
custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# 2) PLOT A: ClinVar Lollipop Plot
plot_clinvar <- ggplot() +
  
  # ✅ Add colored boxes for domains BELOW the lollipops
  geom_rect(
    data = chd3_domains,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0, fill = Name),
    alpha = 0.4,
    color = "black"
  ) +
  
  # ✅ Add lollipops for ClinVar variants
  geom_segment(
    data = clinvar_counts,
    aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
    color = "red"
  ) +
  
  geom_point(
    data = clinvar_counts,
    aes(x = PROTEIN_POS, y = Count),
    color = "red",
    size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous("Number of Variants", limits = c(-domain_height, max(clinvar_counts$Count, na.rm = TRUE) + 1)) +
  ggtitle("ClinVar Missense Variants in chd3") +
  custom_theme

ggsave("lollipop_ClinVar_chd3.png", plot_clinvar, width = 10, height = 4, dpi = 300)

# 3) PLOT B: gnomAD Lollipop Plot
plot_gnomad <- ggplot() +
  
  # ✅ Add colored boxes for domains BELOW the lollipops
  geom_rect(
    data = chd3_domains,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0, fill = Name),
    alpha = 0.4,
    color = "black"
  ) +
  
  # ✅ Add lollipops for gnomAD variants
  geom_segment(
    data = gnomad_counts,
    aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
    color = "black"
  ) +
  
  geom_point(
    data = gnomad_counts,
    aes(x = PROTEIN_POS, y = Count),
    color = "black",
    size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous("Number of Variants", limits = c(-domain_height, max(gnomad_counts$Count, na.rm = TRUE) + 1)) +
  ggtitle("gnomAD Missense Variants in chd3") +
  custom_theme

ggsave("lollipop_gnomAD_chd3.png", plot_gnomad, width = 10, height = 4, dpi = 300)

message("✅ Lollipop plots saved successfully!")
