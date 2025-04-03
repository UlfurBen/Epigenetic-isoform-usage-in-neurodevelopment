library(dplyr)
library(readr)
library(tidyr)

# Load the exon-level Fisher results (which already include a/c and isoform type)
merged <- read_csv("exon_fisher_enrichment_results.csv", show_col_types = FALSE)

# Check if needed columns are present
stopifnot(all(c("gene", "isoform_type", "a", "c") %in% colnames(merged)))

# Summarize ClinVar and gnomAD counts by gene and isoform type
summary_by_gene <- merged %>%
  group_by(gene, isoform_type) %>%
  summarize(
    clinvar = sum(a, na.rm = TRUE),
    gnomad = sum(c, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to wide format to separate canonical and non-canonical counts per gene
wide_by_gene <- summary_by_gene %>%
  pivot_wider(
    names_from = isoform_type,
    values_from = c(clinvar, gnomad),
    values_fill = 0  # fill with 0 if a gene lacks one exon type
  )

# Optional: Keep only genes with both isoform types
wide_by_gene <- wide_by_gene %>%
  filter(!is.na(clinvar_canonical), !is.na(clinvar_non_canonical))

# Run Fisher's exact test per gene
fisher_results_per_gene <- wide_by_gene %>%
  rowwise() %>%
  mutate(
    fisher_p = tryCatch({
      mat <- matrix(c(
        clinvar_non_canonical, clinvar_canonical,
        gnomad_non_canonical, gnomad_canonical
      ), nrow = 2)
      fisher.test(mat)$p.value
    }, error = function(e) NA_real_)
  ) %>%
  ungroup()

# Add FDR correction for multiple testing
fisher_results_per_gene <- fisher_results_per_gene %>%
  mutate(
    fdr = p.adjust(fisher_p, method = "fdr")
  )

# Add enrichment ratio
fisher_results_per_gene <- fisher_results_per_gene %>%
  mutate(
    clinvar_ratio = clinvar_non_canonical / (clinvar_non_canonical + clinvar_canonical + 1e-6),
    gnomad_ratio  = gnomad_non_canonical / (gnomad_non_canonical + gnomad_canonical + 1e-6),
    enrichment_ratio = clinvar_ratio / gnomad_ratio
  )

# Save output
write_csv(fisher_results_per_gene, "per_gene_fisher_results.csv")
cat("✔ Per-gene Fisher's test results saved to 'per_gene_fisher_results.csv'\n")
