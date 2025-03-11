###############################################################
# R SCRIPT EXAMPLE: Merge Day3, Day6, Day12 with 
# multiple replicates and run IsoformSwitchAnalyzeR 
# using MOUSE GTF & FASTA from HPC cluster
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

# 2) File Paths to HPC Data (Adapt as needed) ---------------------------
# HPC references for mouse:
gtf_path   <- "/proj/hpcdata/Mimir/shared/From_Katrin_For_Kaan/km125_Perturbseq/genomes/Mus_musculus.GRCm38.98.gtf.gz"
cdna_path  <- "/proj/hpcdata/Mimir/shared/kmoller/km127_RNAseq/Mus_musculus.GRCm38.cdna.all.fa.gz"
# If your day3, day6, day12 .tsv files are also in HPC directory:
day3_file  <- "unfiltered_transcript_counts_with_genes_day3.tsv"
day6_file  <- "unfiltered_transcript_counts_with_genes_day6.tsv"
day12_file <- "unfiltered_transcript_counts_with_genes_day12.tsv"

# 3) Read Each Day's File Without Summation -----------------------------
# Suppose columns are: gene_id, gene_name, feature_id, Day3_A, Day3_B, etc.

day3_counts  <- read.delim(day3_file, header=TRUE)
day6_counts  <- read.delim(day6_file, header=TRUE)
day12_counts <- read.delim(day12_file, header=TRUE)

# 4) Rename feature_id -> isoform_id if present -------------------------
for(df_name in c("day3_counts","day6_counts","day12_counts")) {
  tmp <- get(df_name)
  if("feature_id" %in% colnames(tmp)) {
    names(tmp)[names(tmp) == "feature_id"] <- "isoform_id"
  }
  assign(df_name, tmp)
}

# 5) Subset to isoform_id + replicate columns for each day --------------
# For instance, if day3_counts has columns:
#   isoform_id, gene_id, gene_name, Day3_A, Day3_B, Day3_C, NPC_A, NPC_B, ...
#   (Adjust as needed)

day3_sub <- day3_counts %>%
  dplyr::select(isoform_id, matches("Day3_"))  # Keep replicate columns with "Day3_" in them

day6_sub <- day6_counts %>%
  dplyr::select(isoform_id, matches("Day6_"))

day12_sub <- day12_counts %>%
  dplyr::select(isoform_id, matches("Day12_"))

# You can also keep gene_id or gene_name if you want, 
# but typically we remove them so that each replicate col is numeric.

# 6) Merge All Days by isoform_id ---------------------------------------
# We'll do it step by step
merged_3_6  <- merge(day3_sub, day6_sub,  by="isoform_id", all=TRUE)
all_3_6_12  <- merge(merged_3_6, day12_sub, by="isoform_id", all=TRUE)

# 7) Convert replicate columns to numeric -------------------------------
rep_cols <- setdiff(colnames(all_3_6_12), "isoform_id")
for(col_name in rep_cols) {
  all_3_6_12[[col_name]] <- as.numeric(all_3_6_12[[col_name]])
}
all_3_6_12$isoform_id <- as.character(all_3_6_12$isoform_id)

# 8) Build Design Matrix for All Replicates / All Days ------------------
# Example: Day3 has columns Day3_A, Day3_B, Day3_C 
#          Day6 has Day6_A, Day6_B, Day6_C
#          Day12 has Day12_A, Day12_B, Day12_C
# Adjust if you have different replicate naming.

myDesign <- data.frame(
  sampleID = c("Day3_A","Day3_B","Day3_C",
               "Day6_A","Day6_B","Day6_C",
               "Day12_A","Day12_B","Day12_C"),
  condition = c(rep("Day3",3),
                rep("Day6",3),
                rep("Day12",3))
)

# 9) Decide Pairwise Comparisons or Multi-Group -------------------------
# "2 and 2" might mean you want multiple pairs. 
# For example: Day3 vs Day6, Day3 vs Day12, Day6 vs Day12
# Define them all:

myComparisons <- data.frame(
  condition_1 = c("Day3","Day3","Day6"),
  condition_2 = c("Day6","Day12","Day12")
)

# 10) Import into IsoformSwitchAnalyzeR ---------------------------------
aSwitchList <- importRdata(
  isoformCountMatrix   = all_3_6_12,
  designMatrix         = myDesign,
  isoformExonAnnoation = gtf_path,
  isoformNtFasta       = cdna_path,
  comparisonsToMake    = myComparisons,
  ignoreAfterPeriod    = TRUE,
  ignoreAfterBar       = TRUE,
  ignoreAfterSpace     = TRUE,
  removeNonConvensionalChr = TRUE,
  showProgress         = TRUE
)

# 11) (Optional) Filter or Remove Isoforms Not in GTF -------------------
# If you want to remove isoforms not matching the GTF, you can do the 
# "problematic isoforms" code from the original script. 
# E.g., setdiff(quantified_isoforms, transcript_ids)

# 12) Run isoformSwitchTestDEXSeq (DexSeq Analysis) ---------------------
DiffSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList,
  alpha = 0.1,
  reduceToSwitchingGenes = FALSE  # keep all genes
)

# 13) Inspect Results ---------------------------------------------------
extractSwitchSummary(DiffSwitchListAnalyzed)

# 14) (Optional) Filter for multiple isoforms of interest or genes ---------------
iso_of_interest <- c(
"ENSMUST00000033770",
"ENSMUST00000023165",
"ENSMUST00000176030",
"ENSMUST00000087916",
"ENSMUST00000113573",
"ENSMUST00000099490"
)

# Filter the results to only include isoforms of interest
filtered_isoforms <- DiffSwitchListAnalyzed$isoformSwitchAnalysis %>%
  dplyr::filter(isoform_id %in% iso_of_interest)

# Save results to CSV
write.csv(filtered_isoforms, "isoform_usage_results.csv", row.names = FALSE)

# Print the filtered results
print(filtered_isoforms)


###############################################################
# Done! You have:
#  - Merged day3, day6, day12 raw counts w/o summing replicates
#  - Specified HPC mouse GTF/FASTA
#  - Set up pairwise comparisons among the 3 conditions
#  - Run the isoform usage test
###############################################################
