# Load required libraries
library(biomaRt)
library(dplyr)
library(tidyr)
library(readr)

# Connect to Ensembl BioMart using the "www" mirror
ensembl_mart <- useEnsembl(biomart = "ensembl", 
                           dataset = "hsapiens_gene_ensembl", 
                           mirror = "useast")

# Define the list of target genes
target_genes <- read_lines("gene_symbols.csv")

# Retrieve Ensembl Gene IDs
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = target_genes,
  mart = ensembl_mart
)

if (nrow(gene_info) == 0) {
  stop("No matching genes found in Ensembl. Check spelling or dataset availability.")
}

# Retrieve transcript data
transcript_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

# Identify canonical transcripts
canonical_transcripts <- transcript_data %>%
  group_by(ensembl_gene_id) %>%
  slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(ensembl_transcript_id)

# Retrieve exon data
exon_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
                 "exon_chrom_start", "exon_chrom_end", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

# Merge exon and gene info
exons_joined <- merge(exon_data, gene_info, by = "ensembl_gene_id")

# Rename columns
colnames(exons_joined)[colnames(exons_joined) == "external_gene_name"] <- "gene"
colnames(exons_joined)[colnames(exons_joined) == "chromosome_name"] <- "exon_chr"
colnames(exons_joined)[colnames(exons_joined) == "exon_chrom_start"] <- "exon_start"
colnames(exons_joined)[colnames(exons_joined) == "exon_chrom_end"] <- "exon_end"

# Convert exon start and end to numeric
exons_joined$exon_start <- as.numeric(exons_joined$exon_start)
exons_joined$exon_end   <- as.numeric(exons_joined$exon_end)

# Find exons unique to one isoform
unique_exons <- exons_joined %>%
  group_by(ensembl_exon_id, gene, exon_chr, exon_start, exon_end) %>%
  summarize(
    n_transcripts = n_distinct(ensembl_transcript_id),
    isoform_type = ifelse(any(ensembl_transcript_id %in% canonical_transcripts$ensembl_transcript_id),
                          "canonical", "non_canonical"),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_transcripts == 1)

# Now correctly check overlaps within each gene
non_overlapping_exons <- unique_exons %>%
  group_by(gene, exon_chr) %>%
  arrange(exon_start, .by_group = TRUE) %>%
  mutate(
    overlap = exon_start < lag(exon_end, default = -Inf)
  ) %>%
  ungroup() %>%
  filter(!overlap)

# Save outputs
write.csv(non_overlapping_exons, "non_overlapping_unique_exons.csv", row.names = FALSE)

canonical_non_overlapping_exons_only <- non_overlapping_exons %>%
  filter(isoform_type == "canonical")
write.csv(canonical_non_overlapping_exons_only, "canonical_non_overlapping_unique_exons.csv", row.names = FALSE)

cat("✔ All non-overlapping unique exons saved to 'non_overlapping_unique_exons.csv'\n")
cat("✔ Canonical non-overlapping unique exons saved to 'canonical_non_overlapping_unique_exons.csv'\n")














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

  # Keep only pathogenic and likely pathogenic variants
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

# Define input VCF file and read gene list
vcf_file <- "clinvar.vcf"
genes <- read_lines("gene_symbols.csv")

# Process each gene
gene_results <- lapply(genes, function(gene) extract_clinvar_missense(vcf_file, gene))




















# Count clinvar variants in each unique exon                       

# Load required libraries
library(dplyr)
library(tidyr)                       
library(readr)  # For read_csv/read_delim

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read_csv("non_overlapping_unique_exons.csv")

target_genes <- read_lines("gene_symbols.csv")

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
  
  write.csv(output, "whole_genome_non_overlapping_exon_clinvar_variant_counts.csv", row.names = FALSE)
  cat("Exon ClinVar variant counts and densities have been saved to 'whole_genome_non_overlapping_exon_clinvar_variant_counts.csv'.\n")
} else {
  cat("No exons with Pathogenic or Likely Pathogenic variants were found.\n")
}














# Get gnomAD missense variants
                       
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
    dplyr::select(variant_id, transcript_id, hgvsc, hgvsp, rsids)
    
  # Convert list-type columns to comma-separated strings
  variants <- variants %>%
    mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))
    
  # Split variant ID
  variants <- variants %>%
    separate(variant_id, into = c("Chromosome", "Position", "Reference", "Alternate"), sep = "-", remove = FALSE) %>%
    dplyr::select(Chromosome, Position, Reference, Alternate, transcript_id, hgvsc, hgvsp, rsids, variant_id)

  # Save to file
  write.csv(variants, file = file_name, row.names = FALSE)
  message("✅ Saved results to: ", file_name)
  return(variants)
}

# Read gene list
genes <- read_lines("gene_symbols.csv")

# Fetch missense variants for each gene, skipping existing files
missense_variants_list <- lapply(genes, get_missense_variants_gnomad)

# Combine results into a single dataframe (excluding NULLs from skipped genes)
missense_variants_df <- bind_rows(Filter(Negate(is.null), missense_variants_list))











                       




# GnomAd missense variant count in each unique exon

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)  # For read_csv

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read_csv("non_overlapping_unique_exons.csv")

target_genes <- read_lines("gene_symbols.csv")

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
  
  write.csv(output, "whole_genome_non_overlapping_exon_gnomad_variant_counts.csv", row.names = FALSE)
  cat("Exon gnomAD variant counts and densities have been saved to 'whole_genome_non_overlapping_exon_gnomad_variant_counts.csv'.\n")
} else {
  cat("No exons with variants were found in the gnomAD files.\n")
}













                       

library(dplyr)
library(readr)

# Load data
clinvar_df <- read_csv("whole_genome_non_overlapping_exon_clinvar_variant_counts.csv", show_col_types = FALSE) %>%
  mutate(gene = trimws(gene))
gnomad_df <- read_csv("whole_genome_non_overlapping_exon_gnomad_variant_counts.csv", show_col_types = FALSE) %>%
  mutate(gene = trimws(gene))
all_exons <- read_csv("non_overlapping_unique_exons.csv", show_col_types = FALSE)

# Merge exon-level variant counts
merged_variants <- full_join(clinvar_df, gnomad_df,
                             by = c("gene", "ensembl_exon_id", "exon_chr", "exon_start", "exon_end", "isoform_type"),
                             suffix = c("_clinvar", "_gnomad")) %>%
  mutate(
    variant_count_clinvar = replace_na(variant_count_clinvar, 0),
    variant_count_gnomad = replace_na(variant_count_gnomad, 0)
  )

# Calculate per-gene total counts
total_variant_counts <- merged_variants %>%
  group_by(gene) %>%
  summarize(
    total_clinvar_gene = sum(variant_count_clinvar, na.rm = TRUE),
    total_gnomad_gene = sum(variant_count_gnomad, na.rm = TRUE),
    .groups = "drop"
  )

# Merge totals into exon data
merged_variants <- merged_variants %>%
  left_join(total_variant_counts, by = "gene") %>%
  mutate(
    a = variant_count_clinvar,
    b = total_clinvar_gene - a,
    c = variant_count_gnomad,
    d = total_gnomad_gene - c
  )

# Perform Fisher's exact test
fisher_results <- merged_variants %>%
  rowwise() %>%
  mutate(
    fisher_output = list(
      tryCatch({
        mat <- matrix(c(a, b, c, d), nrow = 2)
        fisher.test(mat)
      }, error = function(e) NULL)
    ),
    fisher_p = if (!is.null(fisher_output)) fisher_output$p.value else NA_real_,
    odds_ratio = if (!is.null(fisher_output)) unname(fisher_output$estimate[[1]]) else NA_real_
  ) %>%
  ungroup()

# Apply FDR correction
fisher_results <- fisher_results %>%
  mutate(fdr = p.adjust(fisher_p, method = "fdr"))

# Select relevant columns
fisher_results_limited <- fisher_results %>%
  dplyr::select(
    gene, ensembl_exon_id, isoform_type,
    variant_count_clinvar = a,
    variant_count_gnomad = c,
    fisher_p, fdr, odds_ratio
  )

# Save full exon-level results
write_csv(fisher_results_limited, "whole_genome_non_overlapping_exon_fisher_enrichment_results.csv")
cat("✔ Fisher's test results saved: 'whole_genome_non_overlapping_exon_fisher_enrichment_results.csv'\n")

# Find significant genes (FDR < 0.05, best exon per gene)
significant_genes <- fisher_results_limited %>%
  filter(fdr < 0.05) %>%
  group_by(gene) %>%
  slice_min(order_by = fdr, with_ties = FALSE) %>%
  ungroup()

write_csv(significant_genes, "whole_genome_genes_with_significant_non_overlapping_exons.csv")
cat("✔ Significant genes saved: 'whole_genome_genes_with_significant_non_overlapping_exons.csv'\n")

# Prioritize exons (FDR < 0.05 and odds ratio > 2)
prioritized_exons <- fisher_results_limited %>%
  filter(fdr < 0.05, odds_ratio > 2) %>%
  mutate(avg_rank_score = rank(fdr) + rank(-odds_ratio)) %>%
  arrange(avg_rank_score)

write_csv(prioritized_exons, "whole_genome_non_overlapping_exons_prioritized_by_fdr_and_ratio.csv")
cat("✔ Prioritized exons saved: 'whole_genome_non_overlapping_exons_prioritized_by_fdr_and_ratio.csv'\n")
