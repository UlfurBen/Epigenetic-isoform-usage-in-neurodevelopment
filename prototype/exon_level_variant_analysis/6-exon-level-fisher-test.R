library(dplyr)
library(readr)

# Load variant counts and exon metadata
clinvar_df <- read_csv("exon_clinvar_variant_counts.csv", show_col_types = FALSE)
gnomad_df <- read_csv("exon_gnomad_variant_counts.csv", show_col_types = FALSE)
all_exons <- read_csv("all_exons_canonical_and_noncanonical.csv", show_col_types = FALSE)

# Merge ClinVar and gnomAD counts on exon ID
merged_variants <- full_join(clinvar_df, gnomad_df,
                             by = c("gene", "ensembl_exon_id", "exon_chr", "exon_start", "exon_end", "isoform_type"),
                             suffix = c("_clinvar", "_gnomad")) %>%
  mutate(across(c(variant_count_clinvar, variant_count_gnomad), ~replace_na(., 0)))

# Get gene-level totals for all exons (per source)
total_variant_counts <- merged_variants %>%
  group_by(gene) %>%
  summarize(
    total_clinvar_gene = sum(variant_count_clinvar, na.rm = TRUE),
    total_gnomad_gene = sum(variant_count_gnomad, na.rm = TRUE),
    .groups = "drop"
  )

# Merge gene-level totals
merged_variants <- merged_variants %>%
  left_join(total_variant_counts, by = "gene") %>%
  mutate(
    a = variant_count_clinvar,
    b = total_clinvar_gene - variant_count_clinvar,
    c = variant_count_gnomad,
    d = total_gnomad_gene - variant_count_gnomad
  )

# Run Fisher's Exact Test for each exon
fisher_results <- merged_variants %>%
  rowwise() %>%
  mutate(
    fisher_p = tryCatch({
      mat <- matrix(c(a, b, c, d), nrow = 2)
      fisher.test(mat)$p.value
    }, error = function(e) NA_real_)
  ) %>%
  ungroup()

# Add FDR correction
fisher_results <- fisher_results %>%
  mutate(fdr = p.adjust(fisher_p, method = "fdr"))

# Select only desired columns
fisher_results_limited <- fisher_results %>%
  dplyr::select(
    gene,
    ensembl_exon_id,
    isoform_type,
    variant_count_clinvar = a,
    variant_count_gnomad = c,
    fisher_p,
    fdr
  )

# Save full results (limited columns)
write_csv(fisher_results_limited, "exon_fisher_enrichment_results.csv")
cat("✔ Fisher's exact test results saved to 'exon_fisher_enrichment_results.csv'\n")

# Save top hits (FDR < 0.05), sorted by enrichment
#top_hits <- fisher_results_limited %>%
#  filter(fdr < 0.05) %>%
#  arrange(fdr)
#
#write_csv(top_hits, "fisher_exon_top_clinvar_enriched.csv")
#cat("✔ Fisher's exact top exon results saved to 'exon_top_clinvar_enriched.csv'\n")
