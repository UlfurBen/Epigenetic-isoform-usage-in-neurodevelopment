###############################################################################
# Lollipop Plot of Missense Variants with Pfam Domains and Canonical Exons
###############################################################################

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(biomaRt)

###############################################################################
# 1) Load and filter ClinVar missense variants for CHD3
###############################################################################

clinvar_raw <- read_tsv("clinvar_result_CHD3.txt", show_col_types = FALSE)

clinvar_for_plot <- clinvar_raw %>%
  filter(`Molecular consequence` == "missense variant") %>%
  mutate(
    protein_annotation = str_extract(Name, "\\(p\\.[^\\)]+\\)"),
    PROTEIN_POS = str_extract_all(protein_annotation, "[0-9]+")
  ) %>%
  unnest(PROTEIN_POS) %>%
  mutate(PROTEIN_POS = as.numeric(PROTEIN_POS)) %>%
  filter(!is.na(PROTEIN_POS))

clinvar_counts <- clinvar_for_plot %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n(), .groups = "drop")

###############################################################################
# 2) Load Pfam domain annotations for CHD3
###############################################################################

domains_raw <- read_tsv("entry-matching-Q12873.tsv", show_col_types = FALSE)

domain_data <- domains_raw %>%
  filter(`Source Database` == "pfam", !is.na(Matches)) %>%
  separate_rows(Matches, sep = ",") %>%
  mutate(
    Match_Pos = str_extract(Matches, "[0-9]+\\.\\.[0-9]+"),
    Start = as.numeric(str_extract(Match_Pos, "^[0-9]+")),
    End = as.numeric(str_extract(Match_Pos, "[0-9]+$"))
  ) %>%
  dplyr::select(Name, Type, Start, End) %>%
  arrange(Start)

###############################################################################
# 3) Get Unique Non-Canonical Exons using biomaRt
###############################################################################

ensembl_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = "CHD3",
  mart = ensembl_mart
)

transcript_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

canonical_transcripts <- transcript_data %>%
  group_by(ensembl_gene_id) %>%
  slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(ensembl_transcript_id)

exon_data_full <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "cds_start", "cds_end"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
) %>%
  filter(!is.na(cds_start), !is.na(cds_end)) %>%
  mutate(
    Start = as.numeric(cds_start) %/% 3,
    End = as.numeric(cds_end) %/% 3
  )

# Identify non-canonical exons unique to a single isoform
noncanonical_unique_exons <- exon_data_full %>%
  group_by(ensembl_exon_id) %>%
  summarise(
    Start = min(Start),
    End = max(End),
    n_transcripts = n_distinct(ensembl_transcript_id),
    is_canonical = any(ensembl_transcript_id %in% canonical_transcripts$ensembl_transcript_id),
    .groups = "drop"
  ) %>%
  filter(n_transcripts == 1, !is_canonical) %>%
  arrange(Start) %>%
  mutate(overlap = Start < lag(End, default = -Inf)) %>%
  filter(!overlap) %>%
  mutate(Exon_Label = ensembl_exon_id) %>%
  dplyr::select(Start, End, Exon_Label)

###############################################################################
# 4) Plot lollipop chart with colored domains, exon labels, and legend
###############################################################################

domain_height <- max(clinvar_counts$Count, na.rm = TRUE) * 0.2
label_offset <- domain_height * 0.6
exon_bar_height <- domain_height * 0.5

custom_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )

# Create a named vector of colors for each domain name
domain_colors <- setNames(rainbow(length(unique(domain_data$Name))), unique(domain_data$Name))

plot_lollipop <- ggplot() +
  
  # Canonical exons as low gray bars
  geom_rect(
    data = exon_data,
    aes(xmin = Start, xmax = End, ymin = -domain_height - exon_bar_height, ymax = -domain_height),
    fill = "gray90", color = "black", alpha = 0.8, inherit.aes = FALSE
  ) +
  
  # Exon labels
  geom_text(
    data = exon_data,
    aes(x = (Start + End) / 2, y = -domain_height - exon_bar_height - label_offset * 0.5, label = Exon_Label),
    size = 2.8, angle = 45, hjust = 1, vjust = 1, inherit.aes = FALSE
  ) +
  
  # Pfam domains, colored by Name with legend
  geom_rect(
    data = domain_data,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0, fill = Name),
    color = "black", alpha = 0.5
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
    color = "firebrick", size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous(
    "Number of Variants",
    limits = c(-domain_height - exon_bar_height - label_offset * 4.5, max(clinvar_counts$Count) + 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(name = "Pfam Domains", values = domain_colors) +
  ggtitle("Pathogenic Missense Variants in CHD3 (Pfam Domains + Unique Non-Canonical Exons)") +
  custom_theme

###############################################################################
# 5) Save output
###############################################################################

ggsave("lollipop_ClinVar_CHD3_Pfam_and_Exons.png", plot_lollipop, width = 12, height = 5, dpi = 300)
message("✅ Saved as: lollipop_ClinVar_CHD3_Pfam_and_Exons.png")
