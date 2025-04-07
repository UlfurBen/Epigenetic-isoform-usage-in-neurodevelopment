# Get unique exons for these genes from list

# Load required libraries
library(biomaRt)
library(dplyr)

# Connect to Ensembl BioMart using the "useast" mirror
ensembl_mart <- useEnsembl(biomart = "ensembl", 
                           dataset = "hsapiens_gene_ensembl", 
                           mirror = "useast")

# List attributes containing "external_gene_name" for verification
# all_attrs <- listAttributes(ensembl_mart)
# gene_attrs <- all_attrs[grep("external_gene_name", all_attrs$name, ignore.case = TRUE), ]
# print(gene_attrs)

# Define the list of target genes (using their external gene names)
target_genes <- c("AIRE", "AKAP1", "ALG13", "ASH1L", "ASXL1", "ASXL2", "ASXL3", "ATAD2", "ATAD2B", "ATRX",
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

# Retrieve Ensembl Gene IDs for target genes using external_gene_name
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = target_genes,
  mart = ensembl_mart
)

# Check the retrieved gene_info column names
print(names(gene_info))

# Stop if no genes were found
if (nrow(gene_info) == 0) {
  stop("No matching genes found in Ensembl. Check spelling or dataset availability.")
}

# Retrieve transcript data (including transcript length)
transcript_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

# Identify canonical transcripts (longest transcript per gene)
canonical_transcripts <- transcript_data %>%
  group_by(ensembl_gene_id) %>%
  slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(ensembl_transcript_id)

# Retrieve exon data for all target genes
exon_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
                 "exon_chrom_start", "exon_chrom_end", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

# Merge gene info with exon data to include external gene names
exons_joined <- merge(exon_data, gene_info, by = "ensembl_gene_id")

# Print column names to check what we have
print(names(exons_joined))

# Rename columns using base R
colnames(exons_joined)[colnames(exons_joined) == "external_gene_name"] <- "gene"
colnames(exons_joined)[colnames(exons_joined) == "chromosome_name"] <- "exon_chr"
colnames(exons_joined)[colnames(exons_joined) == "exon_chrom_start"] <- "exon_start"
colnames(exons_joined)[colnames(exons_joined) == "exon_chrom_end"] <- "exon_end"

# (Optional) Verify the new names
print(names(exons_joined))

# Convert exon start and end positions to numeric
exons_joined$exon_start <- as.numeric(exons_joined$exon_start)
exons_joined$exon_end   <- as.numeric(exons_joined$exon_end)

# Group by exon and count how many distinct transcripts each exon appears in.
# Also, annotate isoform type based on whether any of the transcripts is canonical.
unique_exons <- exons_joined %>%
  group_by(ensembl_exon_id, gene, exon_chr, exon_start, exon_end) %>%
  summarize(
    n_transcripts = n_distinct(ensembl_transcript_id),
    isoform_type = ifelse(any(ensembl_transcript_id %in% canonical_transcripts$ensembl_transcript_id),
                          "canonical", "non_canonical"),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_transcripts == 1)  # Keep only exons unique to a single isoform

# Save all unique exons (both canonical and non-canonical)
write.csv(unique_exons, "all_exons_canonical_and_noncanonical.csv", row.names = FALSE)

# Save only canonical unique exons
canonical_exons_only <- unique_exons %>% filter(isoform_type == "canonical")
write.csv(canonical_exons_only, "canonical_unique_exons.csv", row.names = FALSE)

cat("✔ All unique exons saved to 'all_exons_canonical_and_noncanonical.csv'\n")
cat("✔ Canonical unique exons saved to 'canonical_unique_exons.csv'\n")






















# Get clinvar variants for genes from list

# Download clinvar vcf file from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# There download clinvar.vcf.gz and gunzip

# Download clinvar vcf file from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# There download clinvar.vcf.gz and gunzip

library(dplyr)
library(readr)
library(stringr)

# Function to extract missense variants for a given gene
extract_clinvar_missense <- function(vcf_file, gene) {
  # Read VCF file, skipping header lines (lines starting with '#')
  clinvar_data <- read_tsv(vcf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
  
  # Define column names (from VCF format)
  colnames(clinvar_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  
  # Filter rows containing the specific gene
  gene_variants <- clinvar_data %>%
    dplyr::filter(str_detect(INFO, paste0("GENEINFO=", gene, ":"))) %>%
    dplyr::filter(str_detect(INFO, "MC=SO:0001583\\|missense_variant"))  # Keep only missense variants
  
  # Extract Clinical Significance (CLNSIG) from INFO column
  gene_variants <- gene_variants %>%
    mutate(CLNSIG = str_extract(INFO, "CLNSIG=[^;]+")) %>%
    mutate(CLNSIG = str_replace(CLNSIG, "CLNSIG=", ""))  # Remove "CLNSIG=" prefix

    # Filter for pathogenic and likely pathogenic variants
  gene_variants <- gene_variants %>%
    dplyr::filter(str_detect(CLNSIG, "Pathogenic|Likely_pathogenic"))
  
  # Select relevant columns and rename them
  gene_variants <- gene_variants %>%
    dplyr::select(CHROM, POS, ID, REF, ALT, CLNSIG)
  
  # Save results to a gene-specific CSV file
  output_file <- paste0("clinvar_result_", gene, ".csv")
  write_csv(gene_variants, output_file)
  
  message("✅ Saved ", nrow(gene_variants), " variants to: ", output_file)
  return(gene_variants)
}

# Define input VCF file and gene list
vcf_file <- "clinvar.vcf"
genes <- c("AIRE", "AKAP1", "ALG13", "ASH1L", "ASXL1", "ASXL2", "ASXL3", "ATAD2", "ATAD2B", "ATRX",
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

# Process each gene
gene_results <- lapply(genes, function(gene) extract_clinvar_missense(vcf_file, gene))





















# Count clinvar variants in each exon

# Load required libraries
library(dplyr)
library(readr)  # For read_csv/read_delim

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read.csv("all_exons_canonical_and_noncanonical.csv", stringsAsFactors = FALSE)

# Convert exon coordinates to numeric
exons$exon_start <- as.numeric(exons$exon_start)
exons$exon_end   <- as.numeric(exons$exon_end)

# Define target genes (updated list)
 target_genes <- c("AIRE", "AKAP1", "ALG13", "ASH1L", "ASXL1", "ASXL2", "ASXL3", "ATAD2", "ATAD2B", "ATRX",
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

# ---- 2. Process ClinVar Variant Files for Each Gene ----
# Variant files are assumed to be named "clinvar_result_[GENE].csv" and contain:
#   - Chromosome column: "CHROM"
#   - Location column (variant coordinate): "POS"
#   - Classification column: "CLNSIG"

# Create an empty list to store exon variant counts
exon_variant_counts <- list()

for(g in target_genes) {
  variant_file <- paste0("clinvar_result_", g, ".csv")
  
  if(!file.exists(variant_file)) {
    message("File not found: ", variant_file, " ... skipping.")
    next
  }
  
  # Read the variant file
  variants <- read.csv(variant_file, stringsAsFactors = FALSE, header = TRUE)
  
  # Clean up column names: trim any leading/trailing whitespace
  colnames(variants) <- trimws(colnames(variants))
  
  # Check for the required columns
  required_cols <- c("CHROM", "POS", "CLNSIG")
  if(!all(required_cols %in% names(variants))) {
    warning("Variant file ", variant_file, " does not contain the required columns. Skipping.")
    next
  }
  
  # Convert CLNSIG column to a consistent format (lowercase, underscores instead of spaces)
  variants$CLNSIG <- tolower(gsub(" ", "_", trimws(variants$CLNSIG)))
  
  # List of pathogenic classifications
  pathogenic_labels <- c("pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic")
  
  # Filter for Pathogenic, Likely Pathogenic, and Pathogenic/Likely Pathogenic variants
  pathogenic_variants <- variants %>%
    dplyr::filter(CLNSIG %in% pathogenic_labels)
  
  # Debugging: Print message if no pathogenic variants are found
  if(nrow(pathogenic_variants) == 0) {
    message("No Pathogenic/Likely Pathogenic variants found for gene ", g, " in ", variant_file)
    next
  }
  
  # Convert the variant location ("POS") column to numeric
  pathogenic_variants[["POS"]] <- as.numeric(pathogenic_variants[["POS"]])
  
  # Ensure the chromosome column is character
  pathogenic_variants[["CHROM"]] <- as.character(pathogenic_variants[["CHROM"]])
  
  # Get the exon data for the current gene
  gene_exons <- exons %>% dplyr::filter(gene == g)
  
  if(nrow(gene_exons) == 0) {
    message("No exon data for gene ", g, " found in exon file.")
    next
  }
  
  # For each exon, count the number of pathogenic variants that fall within its coordinates
  gene_exon_counts <- gene_exons %>%
    rowwise() %>%
    mutate(
      variant_count = sum(
        pathogenic_variants[["CHROM"]] == as.character(exon_chr) &
          pathogenic_variants[["POS"]] >= exon_start &
          pathogenic_variants[["POS"]] <= exon_end,
        na.rm = TRUE
      ),
      # Calculate exon length (ensure length > 0 to avoid division by zero)
      exon_length = exon_end - exon_start,
      variant_density = ifelse(exon_length > 0, variant_count / exon_length, NA)
    ) %>%
    ungroup()
  
  # Keep only exons with at least one variant
  gene_exon_counts <- gene_exon_counts %>% dplyr::filter(variant_count > 0)
  
  if(nrow(gene_exon_counts) > 0) {
    exon_variant_counts[[g]] <- gene_exon_counts
  }
}

# Combine results from all genes
final_results <- bind_rows(exon_variant_counts)

# Write results to an output CSV file if any exons have variants
if(nrow(final_results) > 0) {
  output <- final_results %>%
    dplyr::select(gene, ensembl_exon_id, exon_chr, exon_start, exon_end, isoform_type, variant_count, variant_density)
  
  write.csv(output, "exon_clinvar_variant_counts.csv", row.names = FALSE)
  cat("Exon ClinVar variant counts and densities have been saved to 'exon_clinvar_variant_counts.csv'.\n")
} else {
  cat("No exons with Pathogenic or Likely Pathogenic variants were found.\n")
}
