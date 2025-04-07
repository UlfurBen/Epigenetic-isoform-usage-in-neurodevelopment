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
