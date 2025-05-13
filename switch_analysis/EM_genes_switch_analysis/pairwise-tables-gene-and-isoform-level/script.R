# Combined Script – Generate 9 Output Files for Differential Expression (EM Genes)

library(readr)
library(dplyr)
library(tidyr)

# Define gene list directly
gene_list <- c("AIRE", "AKAP1", "ALG13", "ASH1L", "ASXL1", "ASXL2", "ASXL3", "ATAD2", "ATAD2B", "ATRX",
              "BAHCC1", "BAHD1", "BAZ1A", "BAZ1B", "BAZ2A", "BAZ2B", "BPTF", "BRD1", "BRD2", "BRD3",
              "BRD4", "BRD7", "BRD8", "BRD9", "BRDT", "BRPF1", "BRPF3", "BRWD1", "BRWD3", "C14orf169",
              "CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8", "CDY1", "CDY2A", "CDYL",
              "CDYL2", "CECR2", "CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6", "CHD7", "CHD8", "CHD9",
              "CREBBP", "CXXC1", "CXXC4", "CXXC5", "DIDO1", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DPF1",
              "DPF2", "DPF3", "EED", "EHMT1", "EHMT2", "EP300", "EP400", "EZH1", "EZH2", "FBXL19", "G2E3",
              "GLYR1", "HDAC1", "HDAC10", "HDAC11", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7",
              "HDAC8", "HDAC9", "HDGF", "HDGFL1", "HDGFRP2", "HDGFRP3", "HIF1AN", "HR", "HSPBAP1", "ING1",
              "ING2", "ING3", "ING4", "ING5", "INO80", "INTS12", "JADE1", "JADE2", "JADE3", "JARID2", "JMJD1C",
              "JMJD4", "JMJD6", "JMJD7", "JMJD8", "KAT2A", "KAT2B", "KAT5", "KAT6A", "KAT6B", "KAT7", "KAT8",
              "KDM1A", "KDM1B", "KDM2A", "KDM2B", "KDM3A", "KDM3B", "KDM4A", "KDM4B", "KDM4C", "KDM4D",
              "KDM4E", "KDM5A", "KDM5B", "KDM5C", "KDM5D", "KDM6A", "KDM6B", "KDM7A", "KDM8", "KMT2A",
              "KMT2B", "KMT2C", "KMT2D", "KMT2E", "KMT5A", "KMT5B", "KMT5C", "L3MBTL1", "L3MBTL2", "L3MBTL3",
              "L3MBTL4", "LBR", "MBD1", "MBD2", "MBD3", "MBD4", "MBD5", "MBD6", "MBTD1", "MECP2", "MINA",
              "MLLT10", "MLLT6", "MORC1", "MORC2", "MORC3", "MORC4", "MPHOSPH8", "MSH6", "MSL3", "MTA1",
              "MTA2", "MTA3", "MTF2", "MUM1", "MUM1L1", "NSD1", "ORC1", "PBRM1", "PHF1", "PHF10", "PHF11",
              "PHF12", "PHF13", "PHF14", "PHF19", "PHF2", "PHF20", "PHF20L1", "PHF21A", "PHF21B", "PHF23",
              "PHF3", "PHF6", "PHF7", "PHF8", "PHIP", "PHRF1", "PRDM1", "PRDM10", "PRDM11", "PRDM12",
              "PRDM13", "PRDM14", "PRDM15", "PRDM16", "PRDM2", "PRDM4", "PRDM5", "PRDM6", "PRDM7", "PRDM8",
              "PRDM9", "PSIP1", "PWWP2A", "PWWP2B", "PYGO1", "PYGO2", "RAG2", "RAI1", "RERE", "RNF17",
              "RSF1", "SCMH1", "SCML2", "SETD1A", "SETD1B", "SETD2", "SETD3", "SETD4", "SETD5", "SETD6",
              "SETD7", "SETD9", "SETDB1", "SETDB2", "SETMAR", "SFMBT1", "SFMBT2", "SHPRH", "SIRT1", "SIRT2",
              "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCA5", "SMN1",
              "SMNDC1", "SMYD1", "SMYD2", "SMYD3", "SMYD4", "SMYD5", "SND1", "SP110", "SP140", "SP140L",
              "SRCAP", "STK31", "SUV39H1", "SUV39H2", "TAF1", "TAF1L", "TAF3", "TCF19", "TCF20", "TDRD1",
              "TDRD10", "TDRD12", "TDRD15", "TDRD3", "TDRD5", "TDRD6", "TDRD7", "TDRD9", "TDRKH", "TET1",
              "TET2", "TET3", "TNRC18", "TRIM24", "TRIM28", "TRIM33", "TRIM66", "TYW5", "UBR7", "UHRF1",
              "UHRF2", "UTY", "WHSC1", "WHSC1L1", "ZCWPW1", "ZCWPW2", "ZMYND11", "ZMYND8")
gene_list_lower <- tolower(gene_list)

# Load preprocessed files
expr_df <- read.csv("isoform_expression_levels_with_npc.csv")
stats_df <- read.csv("isoform_statistical_results_with_npc.csv")

# Filter to EM genes
expr_df <- expr_df %>% filter(tolower(gene_name) %in% gene_list_lower)
stats_df <- stats_df %>% filter(tolower(gene_name) %in% gene_list_lower)

# --- Global Isoform-Level Top 50 ---
top_isoforms <- stats_df %>%
  filter(!is.na(fdr)) %>%
  arrange(fdr) %>%
  slice_head(n = 50)

write.csv(top_isoforms, "top50_differentially_expressed_isoforms_em.csv", row.names = FALSE)

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
  write.csv(top_genes, paste0("top50_genes_em_", g1, "_vs_", g2, ".csv"), row.names = FALSE)

  # Isoform-level
  top_isoforms_pair <- get_top_isoforms(long_df, g1, g2)
  write.csv(top_isoforms_pair, paste0("top50_isoforms_em_", g1, "_vs_", g2, ".csv"), row.names = FALSE)
}

cat("\n✅ All 9 EM gene output tables successfully created.\n")
