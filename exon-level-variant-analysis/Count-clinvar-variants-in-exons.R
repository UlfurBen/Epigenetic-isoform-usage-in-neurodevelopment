# Load required libraries
library(dplyr)
library(readr)  # For read_csv/read_delim

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read.csv("all_exons_canonical_and_noncanonical.csv", stringsAsFactors = FALSE)

# Assumed columns: ensembl_exon_id, gene, exon_chr, exon_start, exon_end, isoform_type
# Convert exon coordinates to numeric
exons$exon_start <- as.numeric(exons$exon_start)
exons$exon_end   <- as.numeric(exons$exon_end)

# Define target genes (updated list)
target_genes <- c("BAHD1", "BAZ1B", "BRD2", "BRD4", "BRPF1", "CHD2", "EHMT1", 
                  "EP400", "KAT6A", "KDM4A", "MBD1", "MECP2", "PHF8", "SETD5", "TAF1")

# ---- 2. Process ClinVar Variant Files for Each Gene ----
# Variant files are assumed to be named "clinvar_result_[GENE].txt" and contain:
#   - Chromosome column: "GRCh38Chromosome"
#   - Location column (variant coordinate): "GRCh38Location"
#   - Classification column: "Germline.classification"
#
# For each exon, we count the number of variants (rows) that fall within its coordinates,
# then calculate variant density = variant_count / (exon_end - exon_start).

# Create an empty list to store exon variant counts
exon_variant_counts <- list()

for(g in target_genes) {
  variant_file <- paste0("clinvar_result_", g, ".txt")
  
  if(!file.exists(variant_file)) {
    message("File not found: ", variant_file, " ... skipping.")
    next
  }
  
  # Read the variant file (assuming tab-delimited)
  variants <- read.delim(variant_file, stringsAsFactors = FALSE, header = TRUE)
  
  # Clean up column names: trim any leading/trailing whitespace
  colnames(variants) <- trimws(colnames(variants))
  
  # (Optional) Debug: print column names
  # cat("Columns in", variant_file, ":\n")
  # print(colnames(variants))
  
  # Check for the required columns
  required_cols <- c("GRCh38Chromosome", "GRCh38Location", "Germline.classification")
  if(!all(required_cols %in% names(variants))) {
    warning("Variant file ", variant_file, " does not contain the required columns. Skipping.")
    next
  }
  
  # Filter for "Likely pathogenic" and "Pathogenic" variants
  variants <- variants %>%
    filter(`Germline.classification` %in% c("Likely pathogenic", "Pathogenic"))
  
  # Convert the location column ("GRCh38Location") to numeric.
  variants[["GRCh38Location"]] <- as.numeric(variants[["GRCh38Location"]])
  
  # Ensure the chromosome column is character
  variants[["GRCh38Chromosome"]] <- as.character(variants[["GRCh38Chromosome"]])
  
  # Get the exon data for the current gene
  gene_exons <- exons %>% filter(gene == g)
  
  if(nrow(gene_exons) == 0) {
    message("No exon data for gene ", g, " found in exon file.")
    next
  }
  
  # For each exon, count the number of variants (rows) that fall within its coordinates
  # and then calculate variant density.
  gene_exon_counts <- gene_exons %>%
    rowwise() %>%
    mutate(
      variant_count = sum(
        variants[["GRCh38Chromosome"]] == as.character(exon_chr) &
          variants[["GRCh38Location"]] >= exon_start &
          variants[["GRCh38Location"]] <= exon_end,
        na.rm = TRUE
      ),
      # Calculate exon length (make sure length > 0 to avoid division by zero)
      exon_length = exon_end - exon_start,
      variant_density = ifelse(exon_length > 0, variant_count / exon_length, NA)
    ) %>%
    ungroup()
  
  # Keep only exons with at least one variant.
  gene_exon_counts <- gene_exon_counts %>% filter(variant_count > 0)
  
  if(nrow(gene_exon_counts) > 0) {
    exon_variant_counts[[g]] <- gene_exon_counts
  }
}

# Combine results from all genes
final_results <- bind_rows(exon_variant_counts)

# Write results to an output CSV file if any exons have variants
if(nrow(final_results) > 0) {
  output <- final_results %>%
    select(gene, ensembl_exon_id, exon_chr, exon_start, exon_end, isoform_type, variant_count, variant_density)
  
  write.csv(output, "exon_clinvar_variant_counts.csv", row.names = FALSE)
  cat("Exon ClinVar variant counts and densities have been saved to 'exon_clinvar_variant_counts.csv'.\n")
} else {
  cat("No exons with Likely pathogenic or Pathogenic variants were found.\n")
}
