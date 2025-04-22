library(readr)
library(dplyr)

# 1. Load exon-level data (tab-delimited!)
exon_data <- read_csv("whole_genome_exons_prioritized_by_fdr_and_ratio.csv", show_col_types = FALSE) %>%
  mutate(gene = trimws(gene))

# 2. Ensure min FDR per gene exists
exon_min_fdr <- exon_data %>%
  group_by(gene) %>%
  summarise(exon_fdr = min(fdr, na.rm = TRUE), .groups = "drop")

# 3. Load isoform-level FDR (this one is comma-delimited)
isoform_data <- read_csv("whole_genome_significant_isoforms_fdr_below_0.05_fc_above_2.csv", show_col_types = FALSE) %>%
  mutate(gene_name = trimws(gene_name))

# 4. Join by gene name (case-insensitive just in case)
merged <- inner_join(
  isoform_data %>% mutate(gene_name_lower = tolower(gene_name)),
  exon_min_fdr %>% mutate(gene_lower = tolower(gene)),
  by = c("gene_name_lower" = "gene_lower")
)

# 5. Compute average FDR and sort
ranked <- merged %>%
  mutate(avg_fdr = (fdr + exon_fdr) / 2) %>%
  arrange(avg_fdr)

# 6. Save result
write_csv(head(ranked, 10), "intersected_genes_ranked_by_avg_fdr.csv")

# 7. Print summary
cat("✅ Intersected and ranked genes saved to 'intersected_genes_ranked_by_avg_fdr.csv'\n")
print(head(ranked %>% dplyr::select(gene_name, fdr, exon_fdr, avg_fdr), 10))
