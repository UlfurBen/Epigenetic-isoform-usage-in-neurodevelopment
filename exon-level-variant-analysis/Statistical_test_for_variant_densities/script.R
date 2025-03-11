# Load required libraries
library(dplyr)
library(ggplot2)

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
  
  # --- Optional Visualization (commented out) ---
  # Bland-Altman plot: differences vs. mean density for each paired exon
  # paired_data <- paired_data %>%
  #   mutate(mean_density = (variant_density.clinvar + variant_density.gnomad) / 2)
  # ggplot(paired_data, aes(x = mean_density, y = density_diff)) +
  #   geom_point() +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  #   theme_bw() +
  #   labs(title = "Bland-Altman Plot of Variant Density",
  #        x = "Mean Variant Density (ClinVar & gnomAD)",
  #        y = "Density Difference (ClinVar - gnomAD)")
  
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
