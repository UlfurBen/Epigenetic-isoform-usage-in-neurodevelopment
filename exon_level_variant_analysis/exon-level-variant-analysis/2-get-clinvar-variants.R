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
genes <- c("MECP2")

# Process each gene
gene_results <- lapply(genes, function(gene) extract_clinvar_missense(vcf_file, gene))

# Combine results into one dataframe (optional)
# all_results <- bind_rows(gene_results)

# Print results
# print(all_results)
