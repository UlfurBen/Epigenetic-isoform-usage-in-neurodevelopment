###############################################################
# R SCRIPT: EM Genes Isoform Expression Analysis with FDR
###############################################################

# 1) Load or Install Packages -------------------------------------------
if(!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
if(!requireNamespace("IsoformSwitchAnalyzeR", quietly=TRUE)) {
  BiocManager::install("IsoformSwitchAnalyzeR")
}
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(tidyr)
library(broom)

# 2) Define EM Gene List ------------------------------------------------
gene_list <- tolower(c("AIRE", "AKAP1", "ALG13", "ASH1L", "ASXL1", "ASXL2", "ASXL3", "ATAD2", "ATAD2B", "ATRX",
                       "BAHCC1", "BAHD1", "BAZ1A", "BAZ1B", "BAZ2A", "BAZ2B", "BPTF", "BRD1", "BRD2", "BRD3",
                       "BRD4", "BRD7", "BRD8", "BRD9", "BRDT", "BRPF1", "BRPF3", "BRWD1", "BRWD3", "C14orf169",
                       "CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8", "CDY1", "CDY2A", "CDYL",
                       "CDYL2", "CECR2", "CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6", "CHD7", "CHD8", "CHD9",
                       "CREBBP", "CXXC1", "CXXC4", "CXXC5", "DIDO1", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DPF1",
                       "DPF2", "DPF3", "EED", "EHMT1", "EHMT2", "EP300", "EP400", "EZH1", "EZH2", "FBXL19", "G2E3",
                       "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDGF"))

# 3) File Paths to HPC Data ---------------------------------------------
gtf_path   <- "/proj/hpcdata/Mimir/shared/From_Katrin_For_Kaan/km125_Perturbseq/genomes/Mus_musculus.GRCm38.98.gtf.gz"
cdna_path  <- "/proj/hpcdata/Mimir/shared/kmoller/km127_RNAseq/Mus_musculus.GRCm38.cdna.all.fa.gz"

day3_file  <- "unfiltered_transcript_counts_with_genes_day3.tsv"
day6_file  <- "unfiltered_transcript_counts_with_genes_day6.tsv"
day12_file <- "unfiltered_transcript_counts_with_genes_day12.tsv"

# 4) Read Counts --------------------------------------------------------
day3_counts  <- read.delim(day3_file, header=TRUE)
day6_counts  <- read.delim(day6_file, header=TRUE)
day12_counts <- read.delim(day12_file, header=TRUE)

# 5) Filter for EM Genes ------------------------------------------------
day3_counts  <- day3_counts %>% filter(tolower(gene_name) %in% gene_list)
day6_counts  <- day6_counts %>% filter(tolower(gene_name) %in% gene_list)
day12_counts <- day12_counts %>% filter(tolower(gene_name) %in% gene_list)

# 6) Rename Columns if Needed -------------------------------------------
for(df_name in c("day3_counts", "day6_counts", "day12_counts")) {
  tmp <- get(df_name)
  if ("feature_id" %in% colnames(tmp)) {
    names(tmp)[names(tmp) == "feature_id"] <- "isoform_id"
  }
  assign(df_name, tmp)
}

# 7) Merge All Time Points ----------------------------------------------
day3_sub <- day3_counts %>% dplyr::select(isoform_id, gene_name, matches("Day3_"))
day6_sub <- day6_counts %>% dplyr::select(isoform_id, gene_name, matches("Day6_"))
day12_sub <- day12_counts %>% dplyr::select(isoform_id, gene_name, matches("Day12_"))

merged_3_6  <- merge(day3_sub, day6_sub,  by=c("isoform_id", "gene_name"), all=TRUE)
all_3_6_12  <- merge(merged_3_6, day12_sub, by=c("isoform_id", "gene_name"), all=TRUE)

# 8) Clean Numeric Values -----------------------------------------------
rep_cols <- setdiff(colnames(all_3_6_12), c("isoform_id", "gene_name"))
for (col_name in rep_cols) {
  all_3_6_12[[col_name]] <- as.numeric(all_3_6_12[[col_name]])
}
all_3_6_12 <- all_3_6_12 %>% drop_na()

# 9) Save Intermediate Table --------------------------------------------
write.csv(all_3_6_12, "em_isoform_expression_levels.csv", row.names = FALSE)
print("✅ EM isoform expression levels saved.")

# 10) Reshape to Long Format --------------------------------------------
long_df <- all_3_6_12 %>%
  pivot_longer(cols = -c(isoform_id, gene_name),
               names_to = "sample",
               values_to = "expression") %>%
  mutate(
    day = case_when(
      grepl("Day3", sample) ~ "Day3",
      grepl("Day6", sample) ~ "Day6",
      grepl("Day12", sample) ~ "Day12",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(day), !is.na(expression))

# 11) ANOVA + Tukey Function --------------------------------------------
get_dominant_day <- function(df) {
  df <- df %>% mutate(day = as.character(day), expression = as.numeric(expression)) %>%
    filter(!is.na(expression))

  if (length(unique(df$day)) < 2) {
    return(data.frame(best_day = NA, p_value = NA))
  }

  model <- tryCatch(aov(expression ~ day, data = df), error = function(e) return(NULL))
  if (is.null(model)) return(data.frame(best_day = NA, p_value = NA))

  tukey <- tryCatch(TukeyHSD(model), error = function(e) return(NULL))
  if (is.null(tukey)) return(data.frame(best_day = NA, p_value = NA))

  results <- as.data.frame(tukey$day)
  results$comparison <- rownames(results)
  sig_results <- results %>% filter(`p adj` < 0.05)
  if (nrow(sig_results) == 0) return(data.frame(best_day = NA, p_value = NA))

  avg_expr <- df %>% group_by(day) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE)) %>%
    arrange(desc(mean_expr)) %>%
    slice_head(n = 1)

  return(data.frame(best_day = avg_expr$day, p_value = min(sig_results$`p adj`)))
}

# 12) Apply Per Isoform + FDR -------------------------------------------
anova_results <- long_df %>%
  group_by(isoform_id, gene_name) %>%
  group_modify(~get_dominant_day(.x)) %>%
  ungroup()

anova_results <- anova_results %>%
  mutate(fdr = p.adjust(p_value, method = "fdr"))

# 13) Save All and Significant Hits -------------------------------------
write.csv(anova_results, "em_significant_isoform_upregulation_by_day_fdr.csv", row.names = FALSE)

significant_isoforms <- anova_results %>%
  filter(!is.na(best_day), fdr < 0.05)

write.csv(significant_isoforms, "em_significant_isoforms_fdr_below_0.05.csv", row.names = FALSE)

# 14) Show Top Hits -----------------------------------------------------
top_hits <- significant_isoforms %>%
  arrange(fdr) %>%
  head(10)

print("✅ Top significantly upregulated EM isoforms (FDR < 0.05):")
print(top_hits)
