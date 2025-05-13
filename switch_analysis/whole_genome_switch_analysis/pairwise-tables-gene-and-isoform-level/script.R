# Combined Script – Generate 9 Output Files for Differential Expression (Whole Genome)

library(readr)
library(dplyr)
library(tidyr)

# Load preprocessed files
expr_df <- read.csv("whole_genome_isoform_expression_levels_with_npc.csv")
stats_df <- read.csv("whole_genome_isoform_statistical_results_with_npc.csv")

# --- Global Isoform-Level Top 50 ---
top_isoforms <- stats_df %>%
  filter(!is.na(fdr)) %>%
  arrange(fdr) %>%
  slice_head(n = 50)

write.csv(top_isoforms, "top50_differentially_expressed_isoforms.csv", row.names = FALSE)

# --- Long Format ---
long_df <- expr_df %>%
  pivot_longer(cols = -c(isoform_id, gene_name), names_to = "sample", values_to = "expression") %>%
  mutate(day = case_when(
    grepl("NPC", sample) ~ "NPC",
    grepl("Day3", sample) ~ "Day3",
    grepl("Day6", sample) ~ "Day6",
    grepl("Day12", sample) ~ "Day12"
  )) %>%
  filter(!is.na(expression))

# --- Gene-Level Summarized Data ---
gene_expr <- long_df %>%
  group_by(gene_name, sample, day) %>%
  summarise(expression = sum(expression), .groups = "drop")

# --- Comparison Setup ---
pairwise_contrasts <- list(c("NPC", "Day3"), c("Day3", "Day6"), c("Day6", "Day12"), c("Day3", "Day12"))

# --- Gene-Level Top 50 by T-test ---
get_top_genes <- function(df, g1, g2) {
  df_filtered <- df %>% filter(day %in% c(g1, g2))
  df_filtered %>%
    group_by(gene_name) %>%
    summarise(
      expr1 = list(expression[day == g1]),
      expr2 = list(expression[day == g2]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      p_value = tryCatch(t.test(unlist(expr1), unlist(expr2))$p.value, error = function(e) NA),
      log2fc = log2(mean(unlist(expr2)) + 1e-6) - log2(mean(unlist(expr1)) + 1e-6)
    ) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
    arrange(fdr) %>%
    slice_head(n = 50) %>%
    select(-expr1, -expr2)
}

# --- Isoform-Level Top 50 by T-test ---
get_top_isoforms <- function(df, g1, g2) {
  df_filtered <- df %>% filter(day %in% c(g1, g2))
  df_filtered %>%
    group_by(isoform_id, gene_name) %>%
    summarise(
      expr1 = list(expression[day == g1]),
      expr2 = list(expression[day == g2]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      p_value = tryCatch(t.test(unlist(expr1), unlist(expr2))$p.value, error = function(e) NA),
      log2fc = log2(mean(unlist(expr2)) + 1e-6) - log2(mean(unlist(expr1)) + 1e-6)
    ) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
    arrange(fdr) %>%
    slice_head(n = 50) %>%
    select(-expr1, -expr2)
}

# --- Output Pairwise Gene and Isoform Comparisons ---
for (contrast in pairwise_contrasts) {
  g1 <- contrast[1]; g2 <- contrast[2]

  # Gene-level
  top_genes <- get_top_genes(gene_expr, g1, g2)
  write.csv(top_genes, paste0("top50_genes_", g1, "_vs_", g2, ".csv"), row.names = FALSE)

  # Isoform-level
  top_isoforms_pair <- get_top_isoforms(long_df, g1, g2)
  write.csv(top_isoforms_pair, paste0("top50_isoforms_", g1, "_vs_", g2, ".csv"), row.names = FALSE)
}

cat("\n✅ All 9 output tables successfully created.\n")
