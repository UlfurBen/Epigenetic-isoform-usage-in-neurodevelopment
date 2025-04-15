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
write.csv(unique_exons, "EM_genes_all_exons_canonical_and_noncanonical.csv", row.names = FALSE)

# Save only canonical unique exons
canonical_exons_only <- unique_exons %>% filter(isoform_type == "canonical")
write.csv(canonical_exons_only, "EM_genes_canonical_unique_exons.csv", row.names = FALSE)

cat("✔ All unique exons saved to 'EM_genes_all_exons_canonical_and_noncanonical.csv'\n")
cat("✔ Canonical unique exons saved to 'EM_genes_canonical_unique_exons.csv'\n")






















# Retrieve ClinVar variants from clinvar VCF file

# Download clinvar vcf file from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# There, download clinvar.vcf.gz and gunzip

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Function to extract missense variants for a given gene
extract_clinvar_missense <- function(vcf_file, gene) {
  # Define output file name for this gene
  output_file <- paste0("clinvar_result_", gene, ".csv")

  # Skip processing if file already exists
  if (file.exists(output_file)) {
    message("⚠️ Skipping ", gene, " — file already exists.")
    return(NULL)
  }

  # Read VCF file, skipping header lines (lines starting with '#')
  clinvar_data <- read_tsv(vcf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
  
  # Define column names (from VCF format)
  colnames(clinvar_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  
  # Filter rows containing the specific gene and missense variants
  gene_variants <- clinvar_data %>%
    filter(str_detect(INFO, paste0("GENEINFO=", gene, ":"))) %>%
    filter(str_detect(INFO, "MC=SO:0001583\\|missense_variant"))
  
  # Extract Clinical Significance (CLNSIG)
  gene_variants <- gene_variants %>%
    mutate(CLNSIG = str_extract(INFO, "CLNSIG=[^;]+")) %>%
    mutate(CLNSIG = str_replace(CLNSIG, "CLNSIG=", ""))

   # Keep only pathogenic and/or likely pathogenic variants
  gene_variants <- gene_variants %>%
  mutate(CLNSIG_clean = tolower(gsub("[\\s,]+", "_", CLNSIG))) %>%
  filter(CLNSIG_clean %in% c("pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"))
  
  # Select relevant columns
  gene_variants <- gene_variants %>%
    select(CHROM, POS, ID, REF, ALT, CLNSIG)
  
  # Save results to file
  write_csv(gene_variants, output_file)
  message("✅ Saved ", nrow(gene_variants), " variants to: ", output_file)
  return(gene_variants)
}

# Define input VCF file and specify gene list
vcf_file <- "clinvar.vcf"
genes <- target_genes

# Process each gene
gene_results <- lapply(genes, function(gene) extract_clinvar_missense(vcf_file, gene))




















# Count clinvar variants in each exon

# Load required libraries
library(dplyr)
library(readr)  # For read_csv/read_delim

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read_csv("EM_genes_all_exons_canonical_and_noncanonical.csv")

# Convert exon coordinates to numeric
exons$exon_start <- as.numeric(exons$exon_start)
exons$exon_end   <- as.numeric(exons$exon_end)

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
  variants <- read_csv(variant_file)
  
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
  
  write.csv(output, "EM_genes_exon_clinvar_variant_counts.csv", row.names = FALSE)
  cat("Exon ClinVar variant counts and densities have been saved to 'EM_genes_exon_clinvar_variant_counts.csv'.\n")
} else {
  cat("No exons with Pathogenic and/or Likely Pathogenic variants were found.\n")
}
























library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(readr)

# Function to get missense variants from gnomAD API for a given gene
get_missense_variants_gnomad <- function(gene, dataset = "gnomad_r3") {
  file_name <- paste0("gnomAD_", gene, ".csv")
  
  # Check if file already exists
  if (file.exists(file_name)) {
    message("⚠️ Skipping ", gene, " — file already exists.")
    return(NULL)
  }
  
  base_url <- "https://gnomad.broadinstitute.org/api"
  
  query <- list(query = paste0('
    {
      gene(gene_symbol: "', gene, '", reference_genome: GRCh38) {
        variants(dataset: ', dataset, ') {
          variant_id
          consequence
          transcript_id
          hgvsc
          hgvsp
          rsids
        }
      }
    }
  '))
  
  response <- POST(url = base_url, body = query, encode = "json")

  # Check for successful status code and correct content type
  if (http_type(response) != "application/json") {
    message("⚠️  Non-JSON response for gene: ", gene, " — skipping.")
    return(data.frame())
  }
  
  # Extract and parse JSON
  response_data <- rawToChar(response$content)
  parsed_data <- tryCatch(
    fromJSON(response_data, flatten = TRUE),
    error = function(e) {
      message("❌ JSON parse error for gene: ", gene)
      return(NULL)
    }
  )
  
  if (is.null(parsed_data) || is.null(parsed_data$data$gene$variants)) {
    message("No variants found for gene: ", gene)
    return(data.frame())
  }

  # Continue processing
  variants <- as.data.frame(parsed_data$data$gene$variants)
    
  if (!"consequence" %in% colnames(variants)) {
    message("Skipping gene: ", gene, " (no 'consequence' field found)")
    return(data.frame())
  }
    
  variants <- variants %>%
    filter(consequence == "missense_variant") %>%
    select(variant_id, transcript_id, hgvsc, hgvsp, rsids)
    
  # Convert list-type columns to comma-separated strings
  variants <- variants %>%
    mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))
    
  # Split variant ID
  variants <- variants %>%
    separate(variant_id, into = c("Chromosome", "Position", "Reference", "Alternate"), sep = "-", remove = FALSE) %>%
    select(Chromosome, Position, Reference, Alternate, transcript_id, hgvsc, hgvsp, rsids, variant_id)

  # Save to file
  write.csv(variants, file = file_name, row.names = FALSE)
  message("✅ Saved results to: ", file_name)
  return(variants)
}

# Specify gene list
genes <- target_genes

# Fetch missense variants for each gene, skipping existing files
missense_variants_list <- lapply(genes, get_missense_variants_gnomad)

# Combine results into a single dataframe (excluding NULLs from skipped genes)
missense_variants_df <- bind_rows(Filter(Negate(is.null), missense_variants_list))
























# GnomAd missense variant count in each unique exon 

# Load required libraries
library(dplyr)
library(readr)  # For read_csv

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read_csv("EM_genes_all_exons_canonical_and_noncanonical.csv", show_col_types = FALSE)

# Assumed columns: ensembl_exon_id, gene, exon_chr, exon_start, exon_end, isoform_type
# Convert exon coordinates to numeric
exons$exon_start <- as.numeric(exons$exon_start)
exons$exon_end   <- as.numeric(exons$exon_end)

# ---- 2. Process gnomAD Variant Files for Each Gene ----
# Variant files are assumed to be named "gnomAD_[GENE].csv" and contain:
#   - Chromosome column: "Chromosome"
#   - Position column: "Position"
#
# For each exon, we count the number of variants (rows) that fall within its coordinates.
# Then we compute the variant density = variant_count / (exon_end - exon_start).

# Create an empty list to store exon variant counts
exon_variant_counts <- list()

for (g in target_genes) {
  variant_file <- paste0("gnomAD_", g, ".csv")
  
  if (!file.exists(variant_file)) {
    message("File not found: ", variant_file, " ... skipping.")
    next
  }
  
  # Read the variant file (assuming csv format)
  variants <- read_csv(variant_file, show_col_types = FALSE)
  
  # Clean up column names: trim any leading/trailing whitespace
  colnames(variants) <- trimws(colnames(variants))
  
  # (Optional) Debug: print column names
  # cat("Columns in", variant_file, ":\n")
  # print(colnames(variants))
  
  # Check for the required columns
  required_cols <- c("Chromosome", "Position")
  if (!all(required_cols %in% names(variants))) {
    warning("Variant file ", variant_file, " does not contain the required columns. Skipping.")
    next
  }
  
  # Convert the Position column to numeric.
  variants[["Position"]] <- as.numeric(variants[["Position"]])
  
  # Ensure the Chromosome column is character.
  variants[["Chromosome"]] <- as.character(variants[["Chromosome"]])
  
  # Get the exon data for the current gene
  gene_exons <- dplyr::filter(exons, gene == g)
  
  if (nrow(gene_exons) == 0) {
    message("No exon data for gene ", g, " found in exon file.")
    next
  }
  
  # For each exon, count the number of variants (rows) that fall within its coordinates.
  gene_exon_counts <- gene_exons %>%
    rowwise() %>%
    mutate(
      variant_count = sum(
        variants[["Chromosome"]] == as.character(exon_chr) &
          variants[["Position"]] >= exon_start &
          variants[["Position"]] <= exon_end,
        na.rm = TRUE
      ),
      # Calculate variant density = variant_count divided by exon length.
      variant_density = variant_count / (exon_end - exon_start)
    ) %>%
    ungroup()
  
  # Keep only exons with at least one variant.
  gene_exon_counts <- gene_exon_counts %>% dplyr::filter(variant_count > 0)
  
  if (nrow(gene_exon_counts) > 0) {
    exon_variant_counts[[g]] <- gene_exon_counts
  }
}

# Combine results from all genes
final_results <- bind_rows(exon_variant_counts)

# Write results to an output CSV file if any exons have variants
if (nrow(final_results) > 0) {
  output <- final_results %>%
    dplyr::select(gene, ensembl_exon_id, exon_chr, exon_start, exon_end, isoform_type, variant_count, variant_density)
  
  write.csv(output, "EM_genes_exon_gnomad_variant_counts.csv", row.names = FALSE)
  cat("Exon gnomAD variant counts and densities have been saved to 'EM_genes_exon_gnomad_variant_counts.csv'.\n")
} else {
  cat("No exons with variants were found in the gnomAD files.\n")
}



















# Fisher test 

library(dplyr)
library(readr)

# Load variant counts and exon metadata
clinvar_df <- read_csv("EM_genes_exon_clinvar_variant_counts.csv", show_col_types = FALSE)
gnomad_df <- read_csv("EM_genes_exon_gnomad_variant_counts.csv", show_col_types = FALSE)
all_exons <- read_csv("EM_genes_all_exons_canonical_and_noncanonical.csv", show_col_types = FALSE)

# Merge ClinVar and gnomAD counts on exon ID
merged_variants <- full_join(clinvar_df, gnomad_df,
                             by = c("gene", "ensembl_exon_id", "exon_chr", "exon_start", "exon_end", "isoform_type"),
                             suffix = c("_clinvar", "_gnomad")) %>%
  mutate(across(c(variant_count_clinvar, variant_count_gnomad), ~replace_na(., 0)))

# Get gene-level totals for all exons (per source)
total_variant_counts <- merged_variants %>%
  group_by(gene) %>%
  summarize(
    total_clinvar_gene = sum(variant_count_clinvar, na.rm = TRUE),
    total_gnomad_gene = sum(variant_count_gnomad, na.rm = TRUE),
    .groups = "drop"
  )

# Merge gene-level totals
merged_variants <- merged_variants %>%
  left_join(total_variant_counts, by = "gene") %>%
  mutate(
    a = variant_count_clinvar,
    b = total_clinvar_gene - variant_count_clinvar,
    c = variant_count_gnomad,
    d = total_gnomad_gene - variant_count_gnomad
  )

# Run Fisher's Exact Test for each exon
fisher_results <- merged_variants %>%
  rowwise() %>%
  mutate(
    fisher_p = tryCatch({
      mat <- matrix(c(a, b, c, d), nrow = 2)
      fisher.test(mat)$p.value
    }, error = function(e) NA_real_)
  ) %>%
  ungroup()

# Add FDR correction
fisher_results <- fisher_results %>%
  mutate(fdr = p.adjust(fisher_p, method = "fdr"))

# Select only desired columns
fisher_results_limited <- fisher_results %>%
  dplyr::select(
    gene,
    ensembl_exon_id,
    isoform_type,
    variant_count_clinvar = a,
    variant_count_gnomad = c,
    fisher_p,
    fdr
  )

# Save full results (limited columns)
write_csv(fisher_results_limited, "EM_genes_exon_fisher_enrichment_results.csv")
cat("✔ Fisher's exact test results saved to 'EM_genes_exon_fisher_enrichment_results.csv'\n")

# Identify genes with at least one significant exon (FDR < 0.05),
# then select the exon with the lowest FDR per gene
significant_genes <- fisher_results_limited %>%
  filter(fdr < 0.05) %>%
  group_by(gene) %>%
  slice_min(order_by = fdr, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(min_fdr_for_gene = fdr)

# Count how many unique genes were retained
cat("✔ Number of genes with at least one exon significantly enriched (FDR < 0.05):", nrow(significant_genes), "\n")

# Save final filtered table
write_csv(significant_genes, "all_EM_genes_with_significant_exons.csv")
