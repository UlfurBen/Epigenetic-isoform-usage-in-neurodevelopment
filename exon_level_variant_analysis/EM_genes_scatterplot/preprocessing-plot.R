# This code requires input files from the exon-level-variant-analysis
# It calculates the ratio of clinvar/gnomad missense variants for each exon

# Script 1: Generate CSV for ClinVar/gnomAD ratios with constraint scores

library(dplyr)
library(readr)

# Load inputs
clinvar <- read_csv("EM_genes_non_overlapping_exon_clinvar_variant_counts.csv", show_col_types = FALSE)
gnomad  <- read_csv("EM_genes_non_overlapping_exon_gnomad_variant_counts.csv", show_col_types = FALSE)

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

# Aggregate ClinVar and gnomAD variant counts by isoform type
agg_counts <- merged %>%
  group_by(isoform_type) %>%
  summarise(
    clinvar_total = sum(variant_count_clinvar, na.rm = TRUE),
    gnomad_total = sum(variant_count_gnomad, na.rm = TRUE)
  ) %>%
  mutate(total_variants = clinvar_total + gnomad_total)

# View aggregated counts
print(agg_counts)

# Extract counts
canonical <- filter(agg_counts, isoform_type == "canonical")
noncanonical <- filter(agg_counts, isoform_type == "non-canonical")

# Run proportion test
prop.test(
  x = c(noncanonical$clinvar_total, canonical$clinvar_total),  # successes
  n = c(noncanonical$total_variants, canonical$total_variants),  # trials
  alternative = "greater"  # test if non-canonical has a greater ClinVar proportion
)

# Save to file
write_csv(merged, "EM_genes_non_overlapping_exon_clinvar_gnomad_ratio_output.csv")
cat("✔ Saved ClinVar/gnomAD ratio data with constraint scores to 'EM_genes_non_overlapping_exon_clinvar_gnomad_ratio_output.csv'\n")
