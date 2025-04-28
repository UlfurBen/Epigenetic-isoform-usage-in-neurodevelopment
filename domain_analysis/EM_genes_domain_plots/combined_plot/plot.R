# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(biomaRt)

# -------------------------------
# 1. Load and clean ClinVar data
# -------------------------------
clinvar_raw <- read_tsv("clinvar_result_CHD3.txt", show_col_types = FALSE)

clinvar_variants <- clinvar_raw %>%
  filter(`Molecular consequence` == "missense variant") %>%
  mutate(
    protein_annotation = str_extract(Name, "\\(p\\.[^\\)]+\\)"),
    PROTEIN_POS = str_extract_all(protein_annotation, "[0-9]+")
  ) %>%
  unnest(PROTEIN_POS) %>%
  mutate(PROTEIN_POS = as.numeric(PROTEIN_POS)) %>%
  filter(!is.na(PROTEIN_POS))

# -------------------------------
# 2. Load and clean gnomAD data
# -------------------------------
gnomad_raw <- read_csv("CHD3_gnomad_variants.csv", show_col_types = FALSE)

gnomad_variants <- gnomad_raw %>%
  mutate(PROTEIN_POS = str_extract(hgvsp, "[0-9]+")) %>%
  mutate(PROTEIN_POS = as.numeric(PROTEIN_POS)) %>%
  filter(!is.na(PROTEIN_POS))

# -------------------------------
# 3. Load Pfam domain annotations
# -------------------------------
domains_raw <- read_tsv("entry-matching-Q12873.tsv", show_col_types = FALSE)

domain_data <- domains_raw %>%
  filter(`Source Database` == "pfam", !is.na(Matches)) %>%
  separate_rows(Matches, sep = ",") %>%
  mutate(
    Match_Pos = str_extract(Matches, "[0-9]+\\.\\.[0-9]+"),
    Start = as.numeric(str_extract(Match_Pos, "^[0-9]+")),
    End = as.numeric(str_extract(Match_Pos, "[0-9]+$"))
  ) %>%
  dplyr::select(Name, Start, End)

# -------------------------------
# 4. Load Non-canonical unique exon regions (biomaRt)
# -------------------------------
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
  dplyr::select(Start, End)

# -------------------------------
# 5. Count variants per domain/non-canonical exon
# -------------------------------
count_variants_in_domains_and_exons <- function(variant_df, source_name) {
  variant_with_annotations <- variant_df %>%
    rowwise() %>%
    mutate(
      Domains = paste(domain_data$Name[PROTEIN_POS >= domain_data$Start & PROTEIN_POS <= domain_data$End], collapse = ","),
      In_noncanonical_exon = any(PROTEIN_POS >= noncanonical_unique_exons$Start & PROTEIN_POS <= noncanonical_unique_exons$End)
    ) %>%
    ungroup()
  
  # Variants in domains
  in_domain <- variant_with_annotations %>%
    filter(Domains != "") %>%
    separate_rows(Domains, sep = ",") %>%
    group_by(Domains) %>%
    summarise(Count = n(), .groups = "drop") %>%
    rename(Domain = Domains)
  
  # Variants in non-canonical exons but *not* inside domains
  in_noncanonical_exon <- variant_with_annotations %>%
    filter(Domains == "", In_noncanonical_exon) %>%
    summarise(Count = n()) %>%
    mutate(Domain = "Non-canonical Unique Exon")
  
  # Remaining variants (neither domain nor non-canonical exon)
  non_domain <- variant_with_annotations %>%
    filter(Domains == "", !In_noncanonical_exon) %>%
    summarise(Count = n()) %>%
    mutate(Domain = "Non-domain")
  
  bind_rows(in_domain, in_noncanonical_exon, non_domain) %>%
    mutate(Source = source_name)
}

clinvar_counts <- count_variants_in_domains_and_exons(clinvar_variants, "ClinVar")
gnomad_counts  <- count_variants_in_domains_and_exons(gnomad_variants,  "gnomAD")

# -------------------------------
# 6. Calculate percentages
# -------------------------------
clinvar_total <- nrow(clinvar_variants)
gnomad_total  <- nrow(gnomad_variants)

clinvar_counts <- clinvar_counts %>%
  mutate(Percent = (Count / clinvar_total) * 100)

gnomad_counts <- gnomad_counts %>%
  mutate(Percent = (Count / gnomad_total) * 100)

# Merge for plotting
combined_df <- bind_rows(clinvar_counts, gnomad_counts)

# -------------------------------
# 7. Plot
# -------------------------------
ggplot(combined_df, aes(x = reorder(Domain, -Percent), y = Percent, color = Source)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("ClinVar" = "firebrick", "gnomAD" = "gray40"),
    labels = c("ClinVar" = "ClinVar Pathogenic And/Or Likely Pathogenic Missense Variants", "gnomAD" = "gnomAD")
  ) +
  labs(
    title = "Percent of Missense Variants in Pfam Domains and Non-Canonical Exons (CHD3)",
    x = "Domain / Exon Category",
    y = "Percent of Variants",
    color = "Variant Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Save
ggsave("EM_genes_variant_percent_by_domain_and_noncanoexon.png", width = 11, height = 5, dpi = 300)
