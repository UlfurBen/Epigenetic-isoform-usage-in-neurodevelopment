# This code inputs:
# clinvar variants file for a gene containing protein position which was downloaded manually from the clinvar website,
# domains for the gene downloaded from interpro.
# The code counts the number of variants in each domain and outputs domain variants count with domains labeled.

# Load required libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

###############################################################
# 1) Load ClinVar file and extract protein positions
###############################################################
clinvar_data <- read_tsv("clinvar_result_KMT2A.txt", show_col_types = FALSE)

# Only keep missense variants with valid protein changes
clinvar_variants <- clinvar_data %>%
  filter(`Molecular consequence` == "missense variant") %>%
  mutate(
    protein_annotation = str_extract(Name, "\\(p\\.[^\\)]+\\)"),
    PROTEIN_POS = str_extract_all(protein_annotation, "[0-9]+")
  ) %>%
  unnest(PROTEIN_POS) %>%
  mutate(PROTEIN_POS = as.numeric(PROTEIN_POS)) %>%
  filter(!is.na(PROTEIN_POS))



###############################################################
# 2) Load InterPro domain annotation
###############################################################
domains_raw <- read_tsv("entry-matching-Q03164.tsv", show_col_types = FALSE)

# Extract start and end positions from Match column
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

###############################################################
# 3) Count number of variants that fall in each domain
###############################################################
domain_counts <- domain_data %>%
  rowwise() %>%
  mutate(Variant_Count = sum(clinvar_variants$PROTEIN_POS >= Start &
                             clinvar_variants$PROTEIN_POS <= End, na.rm = TRUE)) %>%
  ungroup()

###############################################################
# 4) Save output
###############################################################
write_tsv(domain_counts, "kmt2a_domain_variant_counts.tsv")
message("✅ Domain variant counts saved to: kmt2a_domain_variant_counts.tsv")
