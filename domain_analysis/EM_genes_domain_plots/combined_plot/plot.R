# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

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
# 4. Count variants per domain (with "Non-domain")
# -------------------------------
count_variants_in_domains <- function(variant_df, source_name) {
  variant_with_domain <- variant_df %>%
    rowwise() %>%
    mutate(
      Domains = paste(domain_data$Name[PROTEIN_POS >= domain_data$Start & PROTEIN_POS <= domain_data$End], collapse = ",")
    ) %>%
    ungroup()
  
  # Split into domain and non-domain
  in_domain <- variant_with_domain %>%
    filter(Domains != "") %>%
    separate_rows(Domains, sep = ",") %>%
    group_by(Domains) %>%
    summarise(Count = n(), .groups = "drop") %>%
    rename(Domain = Domains)
  
  non_domain <- variant_with_domain %>%
    filter(Domains == "" | is.na(Domains)) %>%
    summarise(Count = n()) %>%
    mutate(Domain = "Non-domain")
  
  bind_rows(in_domain, non_domain) %>%
    mutate(Source = source_name)
}

clinvar_counts <- count_variants_in_domains(clinvar_variants, "ClinVar")
gnomad_counts  <- count_variants_in_domains(gnomad_variants,  "gnomAD")

# -------------------------------
# 5. Calculate percentages
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
# 6. Plot
# -------------------------------
ggplot(combined_df, aes(x = reorder(Domain, -Percent), y = Percent, color = Source)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = c("ClinVar" = "firebrick", "gnomAD" = "gray40"),
    labels = c("ClinVar" = "ClinVar Pathogenic And/Or Likely Pathogenic Missense Variants", "gnomAD" = "gnomAD")
  ) +
  labs(
    title = "Percent of Missense Variants in Pfam Domains for CHD3",
    x = "Domain (including Non-domain)",
    y = "Percent of Variants",
    color = "Variant Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# Save
ggsave("EM_genes_variant_percent_by_domain_with_nondomain.png", width = 11, height = 5, dpi = 300)
