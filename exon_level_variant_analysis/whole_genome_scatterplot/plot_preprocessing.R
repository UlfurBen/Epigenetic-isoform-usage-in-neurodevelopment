# This code requires input files from the exon-level-variant-analysis
# It calculates the ratio of clinvar/gnomad missense variants for each exon

# Script 1: Generate CSV for ClinVar/gnomAD ratios with constraint scores

library(dplyr)
library(readr)

# Load inputs
clinvar <- read_csv("whole_genome_exon_clinvar_variant_counts.csv", show_col_types = FALSE)
gnomad  <- read_csv("whole_genome_exon_gnomad_variant_counts.csv", show_col_types = FALSE)

# Merge variant counts on exon coordinates
merged <- full_join(clinvar, gnomad, 
                    by = c("gene", "ensembl_exon_id", "exon_chr", "exon_start", "exon_end", "isoform_type"),
                    suffix = c("_clinvar", "_gnomad")) %>%
  mutate(across(starts_with("variant_count"), ~replace_na(., 0)))

# Calculate ClinVar/gnomAD ratio
merged <- merged %>%
  mutate(
    clinvar_gnomad_ratio = variant_count_clinvar / (variant_count_gnomad + 1e-6)  # avoid division by zero
  )

# Save to file
write_csv(merged, "whole_genome_exon_clinvar_gnomad_ratio_output.csv")
cat("✔ Saved ClinVar/gnomAD ratio data with constraint scores to 'whole_genome_exon_clinvar_gnomad_ratio_output.csv'\n")
