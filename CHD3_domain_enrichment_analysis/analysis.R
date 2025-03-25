###############################################################
## 1) Parse chd3 Domain Data from InterPro TSV File
###############################################################
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

get_chd3_domains_from_tsv <- function(file_path = "entry-matching-Q12873-chd3.tsv") {
  # Read the InterPro TSV file
  df <- read_tsv(file_path, show_col_types = FALSE)

  # Select relevant columns
  domain_df <- df %>%
    dplyr::select(Accession, Name, `Source Database`, Type, Matches)

  # Extract the first start..end match from "Matches"
  domain_df <- domain_df %>%
    mutate(Matches = str_extract(Matches, "[0-9]+\\.\\.[0-9]+")) %>%
    separate(Matches, into = c("Start", "End"), sep = "\\.\\.", convert = TRUE)

  # Convert to numeric and sort
  domain_df <- domain_df %>%
    mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
    arrange(Start)

  # Save to CSV
  write.csv(domain_df, "chd3_domains.csv", row.names = FALSE)
  message("✅ chd3 domain data saved to: chd3_domains.csv")

  return(domain_df)
}

###############################################################
## 2) Extract ClinVar Missense Variants for chd3 (Updated)
###############################################################
extract_clinvar_missense_manual <- function(vcf_file, gene = "chd3") {
  # Read VCF file, skipping header lines
  clinvar_data <- read_tsv(vcf_file, comment = "#", col_names = FALSE, show_col_types = FALSE)

  # Define column names
  colnames(clinvar_data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

  # Debug print
  message("🔎 Total variants read: ", nrow(clinvar_data))

  # Case-insensitive gene filtering
  gene_variants <- clinvar_data %>%
    dplyr::filter(str_detect(INFO, regex(paste0("GENEINFO=", gene, ":"), ignore_case = TRUE))) %>%
    dplyr::filter(str_detect(INFO, "MC=SO:0001583\\|missense_variant"))

  message("✅ Missense variants for ", gene, ": ", nrow(gene_variants))

  # Extract Clinical Significance (CLNSIG)
  gene_variants <- gene_variants %>%
    mutate(CLNSIG = str_extract(INFO, "CLNSIG=[^;]+")) %>%
    mutate(CLNSIG = str_replace(CLNSIG, "CLNSIG=", ""))

  # Extract Protein Position
  gene_variants <- gene_variants %>%
    mutate(PROTEIN_POS = str_extract(INFO, "g\\.[0-9]+[A-Z]>[A-Z]")) %>%
    mutate(PROTEIN_POS = str_extract(PROTEIN_POS, "[0-9]+")) %>%
    mutate(PROTEIN_POS = as.numeric(PROTEIN_POS))

  # Select relevant columns
  gene_variants <- gene_variants %>%
    dplyr::select(CHROM, POS, REF, ALT, PROTEIN_POS, CLNSIG)

  # Save to CSV
  output_file <- "plotA_clinvar_chd3.csv"
  write_csv(gene_variants, output_file)

  message("✅ Saved ", nrow(gene_variants), " variants to: ", output_file)
  return(gene_variants)
}

###############################################################
## 3) Retrieve gnomAD missense variants for chd3 (Unmodified)
###############################################################
library(httr)
library(jsonlite)

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
      message("⚠️ Skipping gene: ", gene, " (no 'consequence' field found)")
      return(data.frame())
    }

    variants <- variants %>%
      dplyr::filter(consequence == "missense_variant") %>%
      dplyr::select(variant_id, transcript_id, hgvsc, hgvsp, rsids)

    variants <- variants %>% 
      mutate(across(where(is.list), ~ sapply(., paste, collapse = ",")))

    if ("variant_id" %in% colnames(variants)) {
      variants <- variants %>%
        tidyr::separate(variant_id, into = c("CHROM", "POS", "REF", "ALT"), sep = "-", remove = FALSE, extra = "drop")
    } else {
      stop("❌ Error: 'variant_id' column not found in gnomAD response.")
    }

    return(variants)
  } else {
    message("⚠️ No variants found for gene: ", gene)
    return(data.frame())
  }
}

###############################################################
## 4) Combine & Save All Variant Data to a Single CSV
###############################################################
retrieve_and_save_chd3_data <- function(
  gene = "chd3",
  clinvar_vcf = "clinvar.vcf",
  gnomad_dataset = "gnomad_r3",
  domain_tsv = "entry-matching-Q12873-chd3.tsv"
) {
  # 1. Get InterPro domains
  domain_data <- get_chd3_domains_from_tsv(domain_tsv)

  # 2. Get ClinVar missense variants
  clinvar_data <- extract_clinvar_missense_manual(clinvar_vcf, gene = gene) %>%
    mutate(Source = "ClinVar")

  # 3. Get gnomAD missense variants
  gnomad_data <- get_missense_variants_gnomad(gene, dataset = gnomad_dataset) %>%
    mutate(Source = "gnomAD")

  # Convert POS to numeric in both
  clinvar_data <- clinvar_data %>%
    mutate(POS = as.numeric(POS))

  gnomad_data <- gnomad_data %>%
    mutate(POS = as.numeric(POS))

  # 4. Combine datasets
  all_variants <- dplyr::bind_rows(
    clinvar_data %>% dplyr::select(CHROM, POS, REF, ALT, Source, PROTEIN_POS, CLNSIG),
    gnomad_data %>%
      mutate(PROTEIN_POS = NA, CLNSIG = NA) %>%
      dplyr::select(CHROM, POS, REF, ALT, Source, PROTEIN_POS, CLNSIG)
  ) %>%
    arrange(POS)

  # Save combined CSV
  variants_file <- "chd3_variants.csv"
  write.csv(all_variants, variants_file, row.names = FALSE)
  message("✅ chd3 variant data saved to: ", variants_file)

  return(list(domains = domain_data, variants = all_variants))
}

###############################################################
## Execute the retrieval
###############################################################
results_list <- retrieve_and_save_chd3_data(
  gene = "chd3",
  clinvar_vcf = "clinvar.vcf",
  gnomad_dataset = "gnomad_r3",
  domain_tsv = "entry-matching-Q12873-chd3.tsv"
)
