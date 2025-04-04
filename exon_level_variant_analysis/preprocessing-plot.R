# This code requires as input files from the exon-level-variant-analysis
# It calculates the ratio of clinvar/gnomad missense variants for each exon

# Script 1: Generate CSV for ClinVar/gnomAD ratios with constraint scores

library(dplyr)
library(readr)

# Load inputs
clinvar <- read_csv("exon_clinvar_variant_counts.csv", show_col_types = FALSE)
gnomad  <- read_csv("exon_gnomad_variant_counts.csv", show_col_types = FALSE)
constraints <- read_tsv("constraint.txt", show_col_types = FALSE)  # or constraint.txt.bgz if unzipped

# Merge variant counts
merged <- full_join(clinvar, gnomad, 
                    by = c("gene", "ensembl_exon_id", "exon_chr", "exon_start", "exon_end", "isoform_type"),
                    suffix = c("_clinvar", "_gnomad")) %>%
  mutate(across(starts_with("variant_count"), ~replace_na(., 0)))

# Add transcript constraint data (isoform-level)
constraint_cols <- c("transcript", "gene", "mis.oe", "mis.z_score", "constraint_flags", "canonical")
constraints_clean <- constraints %>%
  select(any_of(constraint_cols)) %>%
  rename_with(~gsub("\\.", "_", .))  # to match R variable naming

# (Optional) Add a representative transcript ID per exon if you have mapping
# merged <- merged %>% left_join(exon_transcript_map, by = "ensembl_exon_id")

# Merge on gene level only (or use transcript if available)
merged_with_constraints <- merged %>%
  left_join(constraints_clean, by = "gene")

# Calculate ratio
merged_with_constraints <- merged_with_constraints %>%
  mutate(
    clinvar_gnomad_ratio = variant_count_clinvar / (variant_count_gnomad + 1e-6)
  )

# Save to file
write_csv(merged_with_constraints, "exon_clinvar_gnomad_ratio_output.csv")
cat("✔ Saved ClinVar/gnomAD ratio data with constraint scores to 'exon_clinvar_gnomad_ratio_output.csv'\n")
