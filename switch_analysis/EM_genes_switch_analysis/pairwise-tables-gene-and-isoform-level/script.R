# EM_GENES_DE_ANALYSIS.R

library(dplyr)
library(tidyr)
library(readr)

# Load EM expression data (from previous merge)
merged <- read.csv("isoform_expression_levels_with_npc.csv")
expr_cols <- setdiff(names(merged), c("isoform_id", "gene_name"))
merged[expr_cols] <- lapply(merged[expr_cols], as.numeric)
merged <- merged %>% drop_na()

# Long format
long_df <- merged %>%
  pivot_longer(cols = -c(isoform_id, gene_name), names_to = "sample", values_to = "expression") %>%
  mutate(day = case_when(
    grepl("NPC", sample) ~ "NPC",
    grepl("Day3", sample) ~ "Day3",
    grepl("Day6", sample) ~ "Day6",
    grepl("Day12", sample) ~ "Day12"
  )) %>%
  filter(!is.na(expression))

# Statistical test for isoforms
get_significance <- function(df) {
  df <- df %>% filter(!is.na(expression))
  if (length(unique(df$day)) < 2) return(data.frame(p_value = NA, fold_change = NA))
  model <- tryCatch(aov(expression ~ day, data = df), error = function(e) return(NULL))
  if (is.null(model)) return(data.frame(p_value = NA, fold_change = NA))
  tukey <- tryCatch(TukeyHSD(model), error = function(e) return(NULL))
  if (is.null(tukey)) return(data.frame(p_value = NA, fold_change = NA))
  results <- as.data.frame(tukey$day)
  sig <- results %>% filter(`p adj` < 0.05)
  if (nrow(sig) == 0) return(data.frame(p_value = NA, fold_change = NA))
  avg <- df %>% group_by(day) %>% summarise(mean_expr = mean(expression), .groups = "drop")
  fc <- max(avg$mean_expr) / (min(avg$mean_expr) + 1e-6)
  return(data.frame(p_value = min(sig$`p adj`), fold_change = fc))
}

# Isoform-level DE
isoform_stats <- long_df %>%
  group_by(isoform_id, gene_name) %>%
  group_modify(~get_significance(.x)) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_value, method = "fdr"))

# Top 50 isoforms
top_isoforms <- isoform_stats %>%
  filter(!is.na(fdr)) %>%
  arrange(fdr) %>%
  slice_head(n = 50)
write.csv(top_isoforms, "top50_differentially_expressed_isoforms.csv", row.names = FALSE)

# Gene-level DE function
gene_expr <- long_df %>%
  group_by(gene_name, sample, day) %>%
  summarise(expression = sum(expression), .groups = "drop")

pairwise_contrasts <- list(c("NPC", "Day3"), c("Day3", "Day6"), c("Day6", "Day12"), c("Day3", "Day12"))

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
    slice_head(n = 50)
}

# Gene-level tables
for (contrast in pairwise_contrasts) {
  g1 <- contrast[1]; g2 <- contrast[2]
  res <- get_top_genes(gene_expr, g1, g2)
  write.csv(res, paste0("top50_genes_", g1, "_vs_", g2, ".csv"), row.names = FALSE)
}
