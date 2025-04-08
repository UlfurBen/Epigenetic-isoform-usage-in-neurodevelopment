###############################################################
# R SCRIPT EXAMPLE: Merge Day3, Day6, Day12 with 
# multiple replicates and run IsoformSwitchAnalyzeR 
# using MOUSE GTF & FASTA from HPC cluster, filtering for specific genes
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

# 2) Define Gene List --------------------------------------------------
gene_list <- read_csv("gene_symbols.csv", col_names = FALSE, show_col_types = FALSE) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  pull(gene) %>%
  tolower()

# 3) File Paths to HPC Data (Adapt as needed) ---------------------------
gtf_path   <- "/proj/hpcdata/Mimir/shared/From_Katrin_For_Kaan/km125_Perturbseq/genomes/Mus_musculus.GRCm38.98.gtf.gz"
cdna_path  <- "/proj/hpcdata/Mimir/shared/kmoller/km127_RNAseq/Mus_musculus.GRCm38.cdna.all.fa.gz"

day3_file  <- "unfiltered_transcript_counts_with_genes_day3.tsv"
day6_file  <- "unfiltered_transcript_counts_with_genes_day6.tsv"
day12_file <- "unfiltered_transcript_counts_with_genes_day12.tsv"

# 4) Read Each Day's File Without Summation -----------------------------
day3_counts  <- read.delim(day3_file, header=TRUE)
day6_counts  <- read.delim(day6_file, header=TRUE)
day12_counts <- read.delim(day12_file, header=TRUE)

# 5) Filter Data to Include Only Selected Genes (Case Insensitive) ------
day3_counts  <- day3_counts %>% filter(tolower(gene_name) %in% gene_list)
day6_counts  <- day6_counts %>% filter(tolower(gene_name) %in% gene_list)
day12_counts <- day12_counts %>% filter(tolower(gene_name) %in% gene_list)

# 6) Rename feature_id -> isoform_id if present -------------------------
for(df_name in c("day3_counts","day6_counts","day12_counts")) {
  tmp <- get(df_name)
  if("feature_id" %in% colnames(tmp)) {
    names(tmp)[names(tmp) == "feature_id"] <- "isoform_id"
  }
  assign(df_name, tmp)
}

# 7) Merge All Days by isoform_id ---------------------------------------
day3_sub <- day3_counts %>% dplyr::select(isoform_id, gene_name, matches("Day3_"))
day6_sub <- day6_counts %>% dplyr::select(isoform_id, gene_name, matches("Day6_"))
day12_sub <- day12_counts %>% dplyr::select(isoform_id, gene_name, matches("Day12_"))

merged_3_6  <- merge(day3_sub, day6_sub,  by=c("isoform_id", "gene_name"), all=TRUE)
all_3_6_12  <- merge(merged_3_6, day12_sub, by=c("isoform_id", "gene_name"), all=TRUE)

# 8) Convert replicate columns to numeric and Remove NA Rows -------------
rep_cols <- setdiff(colnames(all_3_6_12), c("isoform_id", "gene_name"))
for(col_name in rep_cols) {
  all_3_6_12[[col_name]] <- as.numeric(all_3_6_12[[col_name]])
}
all_3_6_12 <- all_3_6_12 %>% drop_na()

# 9) Save Processed Data for Plotting -----------------------------------
write.csv(all_3_6_12, "whole_genome_isoform_expression_levels.csv", row.names = FALSE)

# Print Summary ----------------------------------------------------------
print("Isoform expression levels saved for plotting.")

###############################################################
# ADDITIONAL SECTION:
# Run ANOVA + Tukey's HSD per isoform to detect significantly 
# upregulated day (Day3, Day6, or Day12) based on expression
###############################################################

# 10) Reshape to long format ---------------------------------------------
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

# 11) Define a function to run ANOVA and TukeyHSD -----------------------
get_dominant_day <- function(df) {
  # Coerce day and expression to basic types
  df <- df %>%
    mutate(day = as.character(day),
           expression = as.numeric(expression)) %>%
    filter(!is.na(expression))

  if (length(unique(df$day)) < 2) {
    return(data.frame(best_day = NA, p_value = NA))
  }

  # Run ANOVA safely
  model <- tryCatch(aov(expression ~ day, data = df), error = function(e) return(NULL))
  if (is.null(model)) return(data.frame(best_day = NA, p_value = NA))

  tukey <- tryCatch(TukeyHSD(model), error = function(e) return(NULL))
  if (is.null(tukey)) return(data.frame(best_day = NA, p_value = NA))

  results <- as.data.frame(tukey$day)
  results$comparison <- rownames(results)

  sig_results <- results %>% filter(`p adj` < 0.05)
  if (nrow(sig_results) == 0) {
    return(data.frame(best_day = NA, p_value = NA))
  }

  avg_expr <- df %>% group_by(day) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE)) %>%
    arrange(desc(mean_expr)) %>%
    slice_head(n = 1)

  return(data.frame(best_day = avg_expr$day, p_value = min(sig_results$`p adj`)))
}


# 12) Apply to each isoform ---------------------------------------------
anova_results <- long_df %>%
  group_by(isoform_id, gene_name) %>%
  group_modify(~get_dominant_day(.x)) %>%
  ungroup()

# 13) Save the results --------------------------------------------------
write.csv(anova_results, "whole_genome_significant_isoform_upregulation_by_day.csv", row.names = FALSE)

# Print top hits
top_hits <- anova_results %>%
  filter(!is.na(best_day)) %>%
  arrange(p_value) %>%
  head(10)

print("Top significantly upregulated isoforms:")
print(top_hits)

write.csv(top_hits, "top_hits_whole_genome_significant_isoform_upregulation_by_day.csv", row.names = FALSE)

                    
