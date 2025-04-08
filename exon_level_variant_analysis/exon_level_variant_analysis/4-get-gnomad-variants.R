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

genes <- c("KMT2A")

# Fetch missense variants for each gene and save to individual files
missense_variants_list <- lapply(genes, get_missense_variants_gnomad)

# Combine results into a single dataframe
missense_variants_df <- bind_rows(missense_variants_list)

# Print results. Only use for small list of genes
# print(missense_variants_df)
