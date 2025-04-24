library(readr)
library(dplyr)

# 1. Load exon-level data prioritized by FDR and odds ratio
exon_data <- read_csv("EM_genes_exons_prioritized_by_fdr_and_oddsratio.csv", show_col_types = FALSE) %>%
  mutate(gene = trimws(gene))

# 2. Filter to top exon hit per gene (highest odds ratio, tie-breaker by lowest FDR)
exon_top_hits <- exon_data %>%
  filter(!is.na(fdr), !is.na(odds_ratio)) %>%
  group_by(gene) %>%
  slice_max(order_by = odds_ratio, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(gene, exon_fdr = fdr, odds_ratio)

# 3. Load isoform-level expression data (significant DE isoforms only)
isoform_data <- read_csv("em_significant_isoforms_fdr_below_0.05_fc_above_2.csv", show_col_types = FALSE) %>%
  mutate(gene_name = trimws(gene_name)) %>%
  filter(!is.na(fdr)) %>%
  group_by(gene_name) %>%
  slice_min(order_by = fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(gene_name, isoform_fdr = fdr)

# 4. Join exon and isoform data by gene (case-insensitive)
merged <- inner_join(
  isoform_data %>% mutate(gene_name_lower = tolower(gene_name)),
  exon_top_hits %>% mutate(gene_lower = tolower(gene)),
  by = c("gene_name_lower" = "gene_lower")
)

# 5. Compute composite score with log protections
ranked <- merged %>%
  mutate(
    composite_score = -log10(isoform_fdr + 1e-10) + -log10(exon_fdr + 1e-10) + log2(odds_ratio + 1e-6),
    min_fdr = pmin(isoform_fdr, exon_fdr)
  ) %>%
  arrange(desc(composite_score))

# 6. Save top results
write_csv(head(ranked, 10), "EM_genes_intersected_genes_ranked_by_composite_score.csv")

# 7. Print summary
cat("✅ Intersected and ranked EM genes saved to 'EM_genes_intersected_genes_ranked_by_composite_score.csv'\n")
print(head(ranked %>% dplyr::select(gene_name, isoform_fdr, exon_fdr, odds_ratio, composite_score), 10))
