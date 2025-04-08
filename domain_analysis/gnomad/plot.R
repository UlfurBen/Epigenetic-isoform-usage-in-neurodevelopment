###############################################################################
# COMBINED SCRIPT: gnomAD Missense Variants Lollipop Plot with InterPro Domains
###############################################################################

# Load required libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(httr)
library(jsonlite)
library(ggplot2)

###############################################################################
# 1) Helper: Load and Parse InterPro Domain TSV (if exists)
###############################################################################

get_kmt2a_domains_from_tsv <- function(file_path = "entry-matching-Q03164.tsv") {
  if (!file.exists(file_path)) {
    stop("❌ Domain file not found: ", file_path)
  }
  
  message("✅ Domain file found, processing...")
  
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  domain_df <- df %>%
    filter(!is.na(Matches)) %>%
    separate_rows(Matches, sep = ",") %>%
    mutate(
      Match_Pos = str_extract(Matches, "[0-9]+\\.\\.[0-9]+"),
      Start = as.numeric(str_extract(Match_Pos, "^[0-9]+")),
      End = as.numeric(str_extract(Match_Pos, "[0-9]+$"))
    ) %>%
    select(Name, Type, Start, End) %>%
    arrange(Start)
  
  write_csv(domain_df, "kmt2a_domains.csv")
  message("✅ Parsed domain data saved to: kmt2a_domains.csv")
  return(domain_df)
}

###############################################################################
# 2) Fetch gnomAD missense variants for KMT2A
###############################################################################

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
    stop("❌ No missense variants found for gene: ", gene)
  }
}

###############################################################################
# 3) Extract protein position from hgvsp
###############################################################################

extract_protein_position <- function(hgvsp_col) {
  pos_str <- str_extract(hgvsp_col, "p\\.[A-Za-z]+([0-9]+)")
  as.numeric(str_extract(pos_str, "[0-9]+"))
}


###############################################################################
# 4) Main Execution Pipeline
###############################################################################

run_gnomad_lollipop_pipeline <- function(
    gene = "KMT2A",
    gnomad_dataset = "gnomad_r3",
    domain_tsv = "entry-matching-Q03164.tsv"
) {
  # Load domain data
  domains_df <- get_kmt2a_domains_from_tsv(domain_tsv)
  
  # Fetch gnomAD missense variants
  variants_df <- get_missense_variants_gnomad(gene, gnomad_dataset)
  
  write_csv(variants_df, "kmt2a_variants.csv")
  message("✅ Raw gnomAD variants saved to: kmt2a_variants.csv")
  
  # Extract protein positions
  variants_df <- variants_df %>%
    mutate(PROTEIN_POS = extract_protein_position(hgvsp)) %>%
    filter(!is.na(PROTEIN_POS))
  
  # Count variants at each position
  variant_counts <- variants_df %>%
    group_by(PROTEIN_POS) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Determine height for domain bars
  domain_height <- max(variant_counts$Count, na.rm = TRUE) * 0.2
  
  # Plot lollipop
  plot <- ggplot() +
    
    # Add domain boxes (gray, no legend)
    geom_rect(
      data = domains_df,
      aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0),
      fill = "gray70",
      alpha = 0.4,
      color = "black"
    ) +
    
    # Lollipop stems
    geom_segment(
      data = variant_counts,
      aes(x = PROTEIN_POS, xend = PROTEIN_POS, y = 0, yend = Count),
      color = "black"
    ) +
    
    # Lollipop heads
    geom_point(
      data = variant_counts,
      aes(x = PROTEIN_POS, y = Count),
      color = "black",
      size = 3
    ) +
    
    scale_x_continuous("Amino Acid Position", expand = c(0, 0)) +
    scale_y_continuous("Number of Variants", limits = c(-domain_height, max(variant_counts$Count) + 1)) +
    ggtitle(paste("gnomAD Missense Variants in", gene)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    )
  
  ggsave("lollipop_gnomAD_kmt2a.png", plot, width = 10, height = 4, dpi = 300)
  message("✅ Lollipop plot saved to: lollipop_gnomAD_kmt2a.png")
}

###############################################################################
# 5) Run everything
###############################################################################

run_gnomad_lollipop_pipeline()
