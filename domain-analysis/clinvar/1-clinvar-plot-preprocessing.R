# This code inputs:
# clinvar variants file for a gene containing protein position,
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
clinvar_data <- read_tsv("clinvar_result_MECP2.txt", show_col_types = FALSE)

# Only keep missense variants with valid protein changes
clinvar_variants <- clinvar_data %>%
  filter(`Molecular consequence` == "missense variant") %>%
  filter(str_detect(`Protein change`, "p\\.")) %>%
  mutate(
    # Extract number from first protein position (e.g. p.Val493Met → 493)
    PROTEIN_POS = as.numeric(str_extract(`Protein change`, "[0-9]+"))
  ) %>%
  filter(!is.na(PROTEIN_POS))

###############################################################
# 2) Load InterPro domain annotation
###############################################################
domains_raw <- read_tsv("entry-matching-P51608.tsv", show_col_types = FALSE)

# Extract start and end positions from Match column
domain_data <- domains_raw %>%
  filter(!is.na(Matches)) %>%
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
write_tsv(domain_counts, "mecp2_domain_variant_counts.tsv")
message("✅ Domain variant counts saved to: mecp2_domain_variant_counts.tsv")
