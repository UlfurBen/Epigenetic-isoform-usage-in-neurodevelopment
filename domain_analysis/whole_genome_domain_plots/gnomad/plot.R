###############################################################################
# gnomAD Missense Variant Lollipop Plot with Pfam Domains and Unique Non-Canonical Exons
###############################################################################

# Load required libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(httr)
library(jsonlite)
library(ggplot2)
library(biomaRt)

###############################################################################
# 1) Load Pfam Domain TSV
###############################################################################

domain_data <- read_tsv("entry-matching-P18507.tsv", show_col_types = FALSE) %>%
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
# 2) Fetch gnomAD Missense Variants for GABRG2 via GraphQL
###############################################################################

get_gnomad_variants <- function(gene = "GABRG2") {
  base_url <- "https://gnomad.broadinstitute.org/api"
  query <- list(query = paste0('{
    gene(gene_symbol: "', gene, '", reference_genome: GRCh38) {
      variants(dataset: gnomad_r3) {
        variant_id
        consequence
        hgvsp
      }
    }
  }'))
  
  res <- POST(base_url, body = query, encode = "json")
  dat <- fromJSON(rawToChar(res$content))
  
  df <- dat$data$gene$variants %>%
    as.data.frame() %>%
    filter(consequence == "missense_variant") %>%
    mutate(
      PROTEIN_POS = as.numeric(str_extract(hgvsp, "[0-9]+"))
    ) %>%
    filter(!is.na(PROTEIN_POS))
  
  return(df)
}

variants_df <- get_gnomad_variants("GABRG2")
variant_counts <- variants_df %>%
  group_by(PROTEIN_POS) %>%
  summarise(Count = n(), .groups = "drop")

###############################################################################
# 3) Get Unique Non-Canonical Exons using biomaRt
###############################################################################

ensembl_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = "GABRG2",
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
# 4) Plot Lollipop with Domains + Unique Non-Canonical Exons
###############################################################################

domain_height <- max(variant_counts$Count, na.rm = TRUE) * 0.2
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

# Unique color per Pfam domain
domain_colors <- setNames(rainbow(length(unique(domain_data$Name))), unique(domain_data$Name))

plot_gnomad <- ggplot() +
  
  # Non-canonical unique exons
  geom_rect(
    data = noncanonical_unique_exons,
    aes(xmin = Start, xmax = End, ymin = -domain_height - exon_bar_height, ymax = -domain_height),
    fill = "gray90", color = "black", alpha = 0.8
  ) +
  
  # Exon labels angled
  geom_text(
    data = noncanonical_unique_exons,
    aes(x = (Start + End) / 2, y = -domain_height - exon_bar_height - label_offset * 0.5, label = Exon_Label),
    size = 2.8, angle = 45, hjust = 1, vjust = 1
  ) +
  
  # Pfam domains colored by name
  geom_rect(
    data = domain_data,
    aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0, fill = Name),
    color = "black", alpha = 0.5
  ) +
  
  # Lollipop stems
  geom_segment(
    data = variant_counts,
    aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
    color = "firebrick"
  ) +
  
  # Lollipop heads
  geom_point(
    data = variant_counts,
    aes(x = PROTEIN_POS, y = Count),
    color = "firebrick", size = 3
  ) +
  
  scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
  scale_y_continuous(
    "Number of Variants",
    limits = c(-domain_height - exon_bar_height - label_offset * 4.5, max(variant_counts$Count) + 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(name = "Pfam Domains", values = domain_colors) +
  ggtitle("gnomAD Missense Variants in GABRG2 (Pfam Domains + Unique Non-Canonical Exons)") +
  custom_theme

###############################################################################
# 5) Save Output
###############################################################################

ggsave("lollipop_gnomAD_GABRG2_noncanonical_exons.png", plot_gnomad, width = 12, height = 5, dpi = 300)
message("✅ Saved as: lollipop_gnomAD_GABRG2_noncanonical_exons.png")
