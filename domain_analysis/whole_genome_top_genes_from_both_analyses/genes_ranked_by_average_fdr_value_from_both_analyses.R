library(readr)
library(dplyr)

# 1. Load exon-level data (prioritized by odds ratio)
exon_data <- read_csv("whole_genome_non_overlapping_exons_prioritized_by_fdr_and_ratio.csv", show_col_types = FALSE) %>%
  mutate(gene = trimws(gene))

# 2. Get the exon with highest odds ratio (or lowest FDR if tied)
exon_top_hits <- exon_data %>%
  filter(!is.na(fdr), !is.na(odds_ratio)) %>%
  group_by(gene) %>%
  slice_max(order_by = odds_ratio, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(gene, exon_fdr = fdr, odds_ratio)

# 3. Load isoform-level expression significance results
isoform_data <- read_csv("whole_genome_isoform_significant_with_npc_fdr_fc.csv", show_col_types = FALSE) %>%
  mutate(gene_name = trimws(gene_name)) %>%
  filter(!is.na(fdr)) %>%
  group_by(gene_name) %>%
  slice_min(order_by = fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(gene_name, isoform_fdr = fdr)

# 4. Join by gene name (case-insensitive)
merged <- inner_join(
  isoform_data %>% mutate(gene_name_lower = tolower(gene_name)),
  exon_top_hits %>% mutate(gene_lower = tolower(gene)),
  by = c("gene_name_lower" = "gene_lower")
)

# 5. Compute composite score with capped log2 odds ratio
ranked <- merged %>%
  mutate(
    # Add small constants to avoid log(0)
    log10_isoform_fdr = -log10(isoform_fdr + 1e-10),
    log10_exon_fdr    = -log10(exon_fdr + 1e-10),
    log2_odds_ratio   = pmin(log2(odds_ratio + 1e-6), 10),  # Cap at log2(1024)
    composite_score   = log10_isoform_fdr + log10_exon_fdr + log2_odds_ratio,
    min_fdr           = pmin(isoform_fdr, exon_fdr)
  ) %>%
  arrange(desc(composite_score))

# 6. Save top results
write_csv(ranked, "intersected_non_overlapping_genes_ranked_by_composite_score.csv")

# 7. Show summary
cat("✅ Intersected genes saved to 'intersected_non_overlapping_genes_ranked_by_composite_score.csv'\n")
print(head(ranked %>% dplyr::select(gene_name, isoform_fdr, exon_fdr, odds_ratio, composite_score), 10))
