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

# 2) Define Gene List --------------------------------------------------
gene_list <- tolower(c("AIRE", "AKAP1", "ALG13", "ASH1L", "ASXL1", "ASXL2", "ASXL3", "ATAD2", "ATAD2B", "ATRX",
                       "BAHCC1", "BAHD1", "BAZ1A", "BAZ1B", "BAZ2A", "BAZ2B", "BPTF", "BRD1", "BRD2", "BRD3",
                       "BRD4", "BRD7", "BRD8", "BRD9", "BRDT", "BRPF1", "BRPF3", "BRWD1", "BRWD3", "C14orf169",
                       "CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8", "CDY1", "CDY2A", "CDYL",
                       "CDYL2", "CECR2", "CHD1", "CHD2", "CHD3", "CHD4", "CHD5", "CHD6", "CHD7", "CHD8", "CHD9",
                       "CREBBP", "CXXC1", "CXXC4", "CXXC5", "DIDO1", "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L", "DPF1",
                       "DPF2", "DPF3", "EED", "EHMT1", "EHMT2", "EP300", "EP400", "EZH1", "EZH2", "FBXL19", "G2E3",
                       "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDGF"))

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
write.csv(all_3_6_12, "isoform_expression_levels.csv", row.names = FALSE)

# Print Summary ----------------------------------------------------------
print("Isoform expression levels saved for plotting.")

###############################################################
# Done! You have:
#  - Filtered isoforms by gene list case insensitively
#  - Merged and processed expression data with gene names
#  - Removed rows containing NA values
#  - Saved expression levels for every isoform of all selected genes
###############################################################
