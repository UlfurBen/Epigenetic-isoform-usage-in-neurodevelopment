# This code uses api to access gnomad missense variants with protein position and domains for that gene.
# Domain data can be found on https://www.ebi.ac.uk/interpro/
# The code counts number of gnomad variants in each domain.
# The code then outputs gnomad variant count for each domain

###############################################################
## Required Libraries
###############################################################
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(httr)
library(jsonlite)

###############################################################
## 1) Parse CHD3 Domain Data from InterPro TSV File
###############################################################
get_chd3_domains_from_tsv <- function(file_path = "entry-matching-Q03164.tsv") {
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  domain_df <- df %>%
    dplyr::select(Accession, Name, `Source Database`, Type, Matches) %>%
    mutate(Matches = str_extract(Matches, "[0-9]+\\.\\.[0-9]+")) %>%
    separate(Matches, into = c("Start", "End"), sep = "\\.\\.", convert = TRUE) %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
    arrange(Start)
  
  write.csv(domain_df, "kmt2a_domains.csv", row.names = FALSE)
  message("✅ kmt2a domain data saved to: chd3_domains.csv")
  
  return(domain_df)
}

###############################################################
## 2) Retrieve gnomAD Missense Variants for kmt2a
###############################################################
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
    
    variants <- variants %>%
      filter(consequence == "missense_variant") %>%
      select(variant_id, transcript_id, hgvsc, hgvsp, rsids) %>%
      mutate(across(where(is.list), ~ sapply(., paste, collapse = ","))) %>%
      separate(variant_id, into = c("CHROM", "POS", "REF", "ALT"), sep = "-", remove = FALSE, extra = "drop")
    
    return(variants)
  } else {
    message("⚠️ No variants found for gene: ", gene)
    return(data.frame())
  }
}

###############################################################
## 3) Extract Protein Position from hgvsp
###############################################################
extract_protein_position <- function(hgvsp_col) {
  # Match number between letters (e.g., p.Arg123Cys or p.R123C)
  position <- str_extract(hgvsp_col, "(?<=p\\.[A-Za-z]+)[0-9]+(?=[A-Za-z]+$)")
  as.numeric(position)
}

###############################################################
## 4) Map Variants to Domains and Count
###############################################################
map_variants_to_domains <- function(variants_df, domains_df) {
  variants_df <- variants_df %>%
    mutate(protein_position = extract_protein_position(hgvsp)) %>%
    filter(!is.na(protein_position))
  
  # Initialize count column
  domains_df$variant_count <- 0
  
  # Count variants in each domain
  for (i in 1:nrow(domains_df)) {
    domain_start <- domains_df$Start[i]
    domain_end <- domains_df$End[i]
    
    count <- sum(variants_df$protein_position >= domain_start & 
                 variants_df$protein_position <= domain_end, na.rm = TRUE)
    
    domains_df$variant_count[i] <- count
  }
  
  write.csv(domains_df, "kmt2a_domain_variant_counts.csv", row.names = FALSE)
  message("✅ Variant counts per domain saved to: kmt2a_domain_variant_counts.csv")
  
  return(domains_df)
}

###############################################################
## 5) Run Full Pipeline
###############################################################
retrieve_and_analyze_chd3_data <- function(
    gene = "kmt2a",
    gnomad_dataset = "gnomad_r3",
    domain_tsv = "entry-matching-Q12873-chd3.tsv"
) {
  domains <- get_chd3_domains_from_tsv(domain_tsv)
  variants <- get_missense_variants_gnomad(gene, gnomad_dataset)
  
  write.csv(variants, "chd3_variants.csv", row.names = FALSE)
  message("✅ Raw variant data saved to: chd3_variants.csv")
  
  final_result <- map_variants_to_domains(variants, domains)
  return(final_result)
}

###############################################################
## 6) Execute
###############################################################
results <- retrieve_and_analyze_chd3_data()
