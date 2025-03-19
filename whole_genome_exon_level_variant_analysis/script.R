# When running script again and again to remove all .csv files run rm *.csv

# Load required libraries
library(biomaRt)
library(dplyr)

# Connect to Ensembl BioMart using the "www" mirror
ensembl_mart <- useEnsembl(biomart = "ensembl", 
                           dataset = "hsapiens_gene_ensembl", 
                           mirror = "www")

# List attributes containing "external_gene_name" for verification
# all_attrs <- listAttributes(ensembl_mart)
# gene_attrs <- all_attrs[grep("external_gene_name", all_attrs$name, ignore.case = TRUE), ]
# print(gene_attrs)

# Define the list of target genes (using their external gene names)
# gene_symbols.csv was created by downloading all gene symbols from https://www.genenames.org/download/custom/?utm_source=chatgpt.com
# and then using R to save first column values to csv file with quotes around each gene symbol.
target_genes <- read_csv("gene_symbols.csv", col_names = FALSE, show_col_types = FALSE) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  pull(gene)

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
  dplyr::filter(n_transcripts == 1)  # Only exons unique to a single isoform

# Save the filtered exon information to a CSV file
write.csv(unique_exons, "all_exons_canonical_and_noncanonical.csv", row.names = FALSE)
cat("Unique exons (appearing in only one isoform) have been saved to 'all_exons_canonical_and_noncanonical.csv'.\n")

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
genes <- target_genes

# Process each gene
gene_results <- lapply(genes, function(gene) extract_clinvar_missense(vcf_file, gene))

# Combine results into one dataframe (optional)
# all_results <- bind_rows(gene_results)

# Print results
# print(all_results)












                       

# Load required libraries
library(dplyr)
library(readr)  # For read_csv/read_delim

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read.csv("all_exons_canonical_and_noncanonical.csv", stringsAsFactors = FALSE)

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
















                       
library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)

# Function to get missense variants from gnomAD API for a given gene
get_missense_variants_gnomad <- function(gene, dataset = "gnomad_r3") {
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
  
  # Extract JSON response
  response_data <- rawToChar(response$content)
  parsed_data <- fromJSON(response_data, flatten = TRUE)
  
  if (!is.null(parsed_data$data$gene$variants)) {
    variants <- as.data.frame(parsed_data$data$gene$variants)
    
    # Check if 'consequence' column exists
    if (!"consequence" %in% colnames(variants)) {
      message("Skipping gene: ", gene, " (no 'consequence' field found)")
      return(data.frame())  # Return an empty dataframe
    }
    
    variants <- variants %>%
      dplyr::filter(consequence == "missense_variant") %>%
      dplyr::select(variant_id, transcript_id, hgvsc, hgvsp, rsids)
    
    # Convert list-type columns to comma-separated strings
    variants <- variants %>% mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))
    
    # Split 'variant_id' into Chromosome, Position, Reference, and Alternate
    variants <- variants %>%
      separate(variant_id, into = c("Chromosome", "Position", "Reference", "Alternate"), sep = "-", remove = FALSE)
    
    # Reorder columns for better readability
    variants <- variants %>%
      dplyr::select(Chromosome, Position, Reference, Alternate, transcript_id, hgvsc, hgvsp, rsids, variant_id)
    
    file_name <- paste0("gnomAD_", gene, ".csv")
    write.csv(variants, file = file_name, row.names = FALSE)
    
    message("Saved results to: ", file_name)
    return(variants)
  } else {
    message("No variants found for gene: ", gene)
    return(data.frame())
  }
}

genes <- c("BRCA1")

# Fetch missense variants for each gene and save to individual files
missense_variants_list <- lapply(genes, get_missense_variants_gnomad)

# Combine results into a single dataframe
missense_variants_df <- bind_rows(missense_variants_list)

# Print results. Only use for small list of genes
# print(missense_variants_df)










                       



                       
# Load required libraries
library(dplyr)
library(readr)  # For read_csv

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read.csv("all_exons_canonical_and_noncanonical.csv", stringsAsFactors = FALSE)

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
  
  # Read the variant file (assuming CSV format)
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
  
  write.csv(output, "exon_gnomad_variant_counts.csv", row.names = FALSE)
  cat("Exon gnomAD variant counts and densities have been saved to 'exon_gnomad_variant_counts.csv'.\n")
} else {
  cat("No exons with variants were found in the gnomAD files.\n")
}














                       

# Load required libraries
library(dplyr)

# Read in the ClinVar and gnomAD exon variant count/density files
clinvar <- read.csv("exon_clinvar_variant_counts.csv", stringsAsFactors = FALSE)
gnomad  <- read.csv("exon_gnomad_variant_counts.csv", stringsAsFactors = FALSE)

# Add a column to indicate dataset
clinvar <- clinvar %>% mutate(dataset = "ClinVar")
gnomad  <- gnomad %>% mutate(dataset = "gnomAD")

# --- Attempt to Pair Data ---
# Merge by common identifiers: ensembl_exon_id, gene, and isoform_type.
paired_data <- merge(clinvar, gnomad, by = c("ensembl_exon_id", "gene", "isoform_type"),
                     suffixes = c(".clinvar", ".gnomad"))

if(nrow(paired_data) > 0) {
  cat("Paired exons found:", nrow(paired_data), "\n")
  
  # Paired test on variant density (if applicable)
  paired_test <- wilcox.test(paired_data$variant_density.clinvar, 
                             paired_data$variant_density.gnomad, 
                             paired = TRUE)
  cat("Paired Wilcoxon signed-rank test result:\n")
  print(paired_test)
  
  # --- Compute Differences ---
  # Calculate the density difference for each paired exon and keep only those with
  # higher ClinVar density than gnomAD (i.e., density_diff > 0).
  paired_data <- paired_data %>%
    mutate(density_diff = variant_density.clinvar - variant_density.gnomad,
           abs_diff = abs(density_diff)) %>%
    dplyr::filter(density_diff > 0)
  
  if(nrow(paired_data) > 0) {
    # Sort exons by the absolute difference in variant density (largest differences first)
    top_exons <- paired_data %>% arrange(desc(abs_diff)) %>% head(10)
    
    cat("Top exons with higher ClinVar variant density than gnomAD:\n")
    print(dplyr::select(top_exons, gene, ensembl_exon_id, isoform_type, 
                        variant_density.clinvar, variant_density.gnomad, density_diff))
    
    # Save the top exons to a separate file
    write.csv(top_exons, "top_exons_high_clinvar_density.csv", row.names = FALSE)
    cat("Top exons saved to 'top_exons_high_clinvar_density.csv'.\n")
  } else {
    cat("No paired exons found with higher ClinVar density than gnomAD.\n")
  }
  
} else {
  cat("No paired exons found. Falling back to unpaired tests.\n")
  
  # --- Unpaired Test Example: Canonical exons ---
  canon_clinvar <- clinvar %>% dplyr::filter(isoform_type == "canonical") %>% pull(variant_density)
  canon_gnomad  <- gnomad %>% dplyr::filter(isoform_type == "canonical") %>% pull(variant_density)
  
  if(length(canon_clinvar) > 0 && length(canon_gnomad) > 0) {
    unpaired_test <- wilcox.test(canon_clinvar, canon_gnomad, paired = FALSE)
    cat("Unpaired Wilcoxon rank-sum test result for canonical exons:\n")
    print(unpaired_test)
  } else {
    cat("Insufficient data for canonical exons for unpaired testing.\n")
  }
  
  # In this fallback, you won't have paired differences per exon.
}

# (Optional) Save the paired data with computed differences for further inspection
if(nrow(paired_data) > 0) {
  write.csv(paired_data, "paired_exon_variant_density_comparison.csv", row.names = FALSE)
  cat("Paired exon data with density differences have been saved to 'paired_exon_variant_density_comparison.csv'.\n")
}
