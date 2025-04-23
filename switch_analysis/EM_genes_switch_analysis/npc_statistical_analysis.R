
###############################################################
# Script 1: Statistical Analysis of Isoform Expression (with NPC)
###############################################################

library(dplyr)
library(tidyr)
library(broom)

# Define gene list
gene_list <- c("CDY2A", "CDYL", "CDYL2", "CECR2", "CHD1", "CHD2", "CHD3")
gene_list_lower <- tolower(gene_list)

# Load expression data (with gene column lowercased)
day3 <- read.delim("unfiltered_transcript_counts_with_genes_day3.tsv", header = TRUE)
day6 <- read.delim("unfiltered_transcript_counts_with_genes_day6.tsv", header = TRUE)
day12 <- read.delim("unfiltered_transcript_counts_with_genes_day12.tsv", header = TRUE)

# Filter by gene list (case-insensitive)
day3 <- day3 %>% filter(tolower(gene_name) %in% gene_list_lower)
day6 <- day6 %>% filter(tolower(gene_name) %in% gene_list_lower)
day12 <- day12 %>% filter(tolower(gene_name) %in% gene_list_lower)

# Subset required columns with full annotation of dplyr::select
day3_sub <- day3 %>% dplyr::select(isoform_id = feature_id, gene_name, NPC_A, NPC_B, NPC_C, Day3_A, Day3_B, Day3_C)
day6_sub <- day6 %>% dplyr::select(isoform_id = feature_id, gene_name, Day6_A, Day6_B, Day6_C)
day12_sub <- day12 %>% dplyr::select(isoform_id = feature_id, gene_name, Day12_A, Day12_B, Day12_C)

# Merge data across all time points
merged <- day3_sub %>%
  full_join(day6_sub, by = c("isoform_id", "gene_name")) %>%
  full_join(day12_sub, by = c("isoform_id", "gene_name"))

# Convert to numeric and drop missing
expr_cols <- setdiff(colnames(merged), c("isoform_id", "gene_name"))
merged[expr_cols] <- lapply(merged[expr_cols], as.numeric)
merged <- merged %>% drop_na()

# Save merged counts
write.csv(merged, "isoform_expression_levels_with_npc.csv", row.names = FALSE)

# Convert to long format
long_df <- merged %>%
  pivot_longer(cols = -c(isoform_id, gene_name), names_to = "sample", values_to = "expression") %>%
  mutate(day = case_when(
    grepl("NPC", sample) ~ "NPC",
    grepl("Day3", sample) ~ "Day3",
    grepl("Day6", sample) ~ "Day6",
    grepl("Day12", sample) ~ "Day12",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(expression))

# Define test function
get_significance <- function(df) {
  df <- df %>% filter(!is.na(expression))
  if (length(unique(df$day)) < 2) return(data.frame(best_day = NA, p_value = NA, fold_change = NA))
  model <- tryCatch(aov(expression ~ day, data = df), error = function(e) return(NULL))
  if (is.null(model)) return(data.frame(best_day = NA, p_value = NA, fold_change = NA))
  tukey <- tryCatch(TukeyHSD(model), error = function(e) return(NULL))
  if (is.null(tukey)) return(data.frame(best_day = NA, p_value = NA, fold_change = NA))
  results <- as.data.frame(tukey$day)
  sig_results <- results %>% filter(`p adj` < 0.05)
  if (nrow(sig_results) == 0) return(data.frame(best_day = NA, p_value = NA, fold_change = NA))
  avg_expr <- df %>% group_by(day) %>% summarise(mean_expr = mean(expression), .groups = "drop")
  top <- avg_expr[which.max(avg_expr$mean_expr), ]
  bottom <- avg_expr[which.min(avg_expr$mean_expr), ]
  fc <- top$mean_expr / (bottom$mean_expr + 1e-6)
  return(data.frame(best_day = top$day, p_value = min(sig_results$`p adj`), fold_change = fc))
}

# Apply per isoform
anova_results <- long_df %>%
  group_by(isoform_id, gene_name) %>%
  group_modify(~get_significance(.x)) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_value, method = "fdr"))

# Write results
write.csv(anova_results, "isoform_statistical_results_with_npc.csv", row.names = FALSE)
significant <- anova_results %>% filter(!is.na(best_day), fdr < 0.05, fold_change > 2)
write.csv(significant, "isoform_significant_with_npc_fdr_fc.csv", row.names = FALSE)
