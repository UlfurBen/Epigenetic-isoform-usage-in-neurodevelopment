###############################################################
## 1) Parse MECP2 Domain Data from InterPro TSV File
###############################################################
library(dplyr)
library(readr)
library(tidyr)

get_mecp2_domains_from_tsv <- function(file_path = "entry-matching-P51608-mecp2.tsv") {
  # Read the InterPro TSV file
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  # Select relevant columns
  domain_df <- df %>%
    dplyr::select(Accession, Name, `Source Database`, Type, Matches)
  
  # Extract start and end positions from "Matches"
  domain_df <- domain_df %>%
    separate(Matches, into = c("Start", "End"), sep = "\\.\\.", convert = TRUE)
  
  # Convert to numeric
  domain_df <- domain_df %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
    arrange(Start)  # Sort by position
  
  # Save to CSV
  write.csv(domain_df, "mecp2_domains.csv", row.names = FALSE)
  message("✅ MECP2 domain data saved to: mecp2_domains.csv")
  
  return(domain_df)
}

###############################################################
## 2) Extract ClinVar Missense Variants for MECP2 (Updated)
###############################################################
library(dplyr)
library(readr)
library(stringr)

extract_clinvar_missense_manual <- function(vcf_file, gene = "MECP2") {
  # Read VCF file, skipping header lines (lines starting with '#')
  clinvar_data <- read_tsv(vcf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
  
  # Define VCF column names
  colnames(clinvar_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  
  # Filter variants related to the gene of interest
  gene_variants <- clinvar_data %>%
    dplyr::filter(str_detect(INFO, paste0("GENEINFO=", gene, ":"))) %>%
    dplyr::filter(str_detect(INFO, "MC=SO:0001583\\|missense_variant"))  # Keep only missense variants
  
  # Extract Clinical Significance (CLNSIG) from INFO column
  gene_variants <- gene_variants %>%
    mutate(CLNSIG = str_extract(INFO, "CLNSIG=[^;]+")) %>%
    mutate(CLNSIG = str_replace(CLNSIG, "CLNSIG=", ""))  # Remove "CLNSIG=" prefix
  
  # Extract Amino Acid Change Position from INFO column
  gene_variants <- gene_variants %>%
    mutate(PROTEIN_POS = str_extract(INFO, "g\\.[0-9]+[A-Z]>[A-Z]")) %>%
    mutate(PROTEIN_POS = str_extract(PROTEIN_POS, "[0-9]+")) %>%
    mutate(PROTEIN_POS = as.numeric(PROTEIN_POS))  # Convert to numeric
  
  # Select relevant columns
  gene_variants <- gene_variants %>%
    dplyr::select(CHROM, POS, REF, ALT, PROTEIN_POS, CLNSIG)
  
  # Save to CSV
  output_file <- "plotA_clinvar_mecp2.csv"
  write_csv(gene_variants, output_file)
  
  message("✅ Saved ", nrow(gene_variants), " variants to: ", output_file)
  return(gene_variants)
}

###############################################################
## 3) Retrieve gnomAD missense variants for MECP2 (Unmodified)
###############################################################
library(httr)
library(jsonlite)
library(tidyr)

get_missense_variants_gnomad <- function(gene, dataset = "gnomad_r3") {
  base_url <- "https://gnomad.broadinstitute.org/api"
  
  query_body <- list(query = paste0('
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
  
  response <- httr::POST(url = base_url, body = query_body, encode = "json")
  response_data <- rawToChar(response$content)
  parsed_data <- fromJSON(response_data, flatten = TRUE)
  
  if (!is.null(parsed_data$data$gene$variants)) {
    variants <- as.data.frame(parsed_data$data$gene$variants)
    
    if (!"consequence" %in% colnames(variants)) {
      message("Skipping gene: ", gene, " (no 'consequence' field found)")
      return(data.frame())
    }
    
    variants <- variants %>%
      dplyr::filter(consequence == "missense_variant") %>%
      dplyr::select(variant_id, transcript_id, hgvsc, hgvsp, rsids)
    
    variants <- variants %>% 
      mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))
    
    # Ensure variant_id exists before splitting
    if ("variant_id" %in% colnames(variants)) {
      variants <- variants %>%
        tidyr::separate(variant_id, into = c("CHROM", "POS", "REF", "ALT"), sep = "-", remove = FALSE, extra = "drop")
    } else {
      stop("Error: 'variant_id' column not found in gnomAD response.")
    }
    
    return(variants)
  } else {
    message("No variants found for gene: ", gene)
    return(data.frame())
  }
}


###############################################################
## 4) Combine & Save All Variant Data to a Single CSV
###############################################################
retrieve_and_save_mecp2_data <- function(
    gene          = "MECP2",
    clinvar_vcf   = "clinvar.vcf",
    gnomad_dataset = "gnomad_r3",
    domain_tsv    = "entry-matching-P51608-mecp2.tsv"
) {
  # ✅ Parse InterPro TSV instead of retrieving Pfam data
  domain_data <- get_mecp2_domains_from_tsv(domain_tsv)

  # Retrieve ClinVar missense variants
  clinvar_data <- extract_clinvar_missense_manual(clinvar_vcf) %>%
    mutate(Source = "ClinVar")

  # Retrieve gnomAD missense variants
  gnomad_data  <- get_missense_variants_gnomad(gene, dataset = gnomad_dataset) %>%
    mutate(Source = "gnomAD")

  # Ensure `POS` column is numeric in both datasets before merging
  clinvar_data <- clinvar_data %>%
    mutate(POS = as.numeric(POS))

  gnomad_data <- gnomad_data %>%
    mutate(POS = as.numeric(POS))

  # ✅ Save gnomAD data separately
  gnomad_file <- "plotB_gnomad_mecp2.csv"
  write.csv(gnomad_data, gnomad_file, row.names = FALSE)
  message("✅ gnomAD missense variants saved to: ", gnomad_file)

  # Combine ClinVar & gnomAD variants
  all_variants <- dplyr::bind_rows(
    clinvar_data %>% dplyr::select(CHROM, POS, REF, ALT, Source, PROTEIN_POS, CLNSIG),
    gnomad_data %>%
      mutate(PROTEIN_POS = NA, CLNSIG = NA) %>%
      dplyr::select(CHROM, POS, REF, ALT, Source, PROTEIN_POS, CLNSIG)
  ) %>%
    arrange(POS)  # Sort by position

  # Save merged variant data
  variants_file <- "mecp2_variants.csv"
  write.csv(all_variants, variants_file, row.names = FALSE)
  message("✅ MECP2 variant data saved to: ", variants_file)

  return(list(domains = domain_data, clinvar_variants = clinvar_data, gnomad_variants = gnomad_data, combined_variants = all_variants))
}




###############################################################
## Execute the retrieval
###############################################################
results_list <- retrieve_and_save_mecp2_data(
  gene         = "MECP2", 
  clinvar_vcf  = "clinvar.vcf",
  gnomad_dataset = "gnomad_r3",
  domain_tsv    = "entry-matching-P51608-mecp2.tsv"
)

# Expected outputs:
#   1) ✅ "mecp2_domains.csv" (MECP2 domain data from InterPro TSV)
#   2) ✅ "plotA_clinvar_mecp2.csv" (ClinVar variants with protein positions)
#   3) ✅ "plotB_gnomad_mecp2.csv" (gnomAD missense variants)
#   4) ✅ "mecp2_variants.csv" (combined ClinVar + gnomAD dataset)
