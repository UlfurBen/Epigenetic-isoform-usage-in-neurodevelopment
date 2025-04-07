# Load required libraries
library(dplyr)
library(readr)  # For read_csv

# ---- 1. Read Exon Data (already contains necessary columns) ----
exons <- read.csv("all_exons_canonical_and_noncanonical.csv", stringsAsFactors = FALSE)

# Assumed columns: ensembl_exon_id, gene, exon_chr, exon_start, exon_end, isoform_type
# Convert exon coordinates to numeric
exons$exon_start <- as.numeric(exons$exon_start)
exons$exon_end   <- as.numeric(exons$exon_end)

# Define target genes
target_genes <- c("A1BG","A1BG-AS1","A1CF","A1S9T","A2M","KMT2A")

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
  variant_file <- paste0("gnomAD_", g, ".tsv.bgz")
  
  if (!file.exists(variant_file)) {
    message("File not found: ", variant_file, " ... skipping.")
    next
  }
  
  # Read the variant file (assuming tsv format)
  variants <- read_tsv(variant_file, show_col_types = FALSE)
  
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
