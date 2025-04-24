
###############################################################
# Script 1: Statistical Analysis of Isoform Expression (with NPC)
###############################################################

library(dplyr)
library(tidyr)
library(broom)

# Define gene list
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
