###############################################################################
# COMBINED SCRIPT: gnomAD Missense Variants Lollipop Plot with Pfam Domains (labeled)
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
# 1) Helper: Load and Parse Pfam Domain TSV
###############################################################################

get_pfam_domains_from_tsv <- function(file_path = "entry-matching-P18507.tsv") {
  if (!file.exists(file_path)) {
    stop("❌ Domain file not found: ", file_path)
  }
  
  message("✅ Domain file found, extracting Pfam domains...")
  
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  domain_df <- df %>%
    filter(`Source Database` == "pfam", !is.na(Matches)) %>%
    separate_rows(Matches, sep = ",") %>%
    mutate(
      Match_Pos = str_extract(Matches, "[0-9]+\\.\\.[0-9]+"),
      Start = as.numeric(str_extract(Match_Pos, "^[0-9]+")),
      End = as.numeric(str_extract(Match_Pos, "[0-9]+$"))
    ) %>%
    dplyr::select(Name, Type, Start, End) %>%
    arrange(Start)
  
  write_csv(domain_df, "GABRG2_pfam_domains.csv")
  message("✅ Parsed Pfam domain data saved to: GABRG2_pfam_domains.csv")
  return(domain_df)
}

###############################################################################
# 2) Fetch gnomAD missense variants for GABRG2
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
      dplyr::select(variant_id, transcript_id, hgvsc, hgvsp, rsids) %>%
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
    gene = "GABRG2",
    gnomad_dataset = "gnomad_r3",
    domain_tsv = "entry-matching-P18507.tsv"
) {
  # Load domain data (Pfam only)
  domain_data <- get_pfam_domains_from_tsv(domain_tsv)
  
  # Fetch gnomAD missense variants
  variants_df <- get_missense_variants_gnomad(gene, gnomad_dataset)
  
  write_csv(variants_df, "GABRG2_gnomad_variants.csv")
  message("✅ Raw gnomAD variants saved to: GABRG2_gnomad_variants.csv")
  
  # Extract protein positions
  variants_df <- variants_df %>%
    mutate(PROTEIN_POS = extract_protein_position(hgvsp)) %>%
    filter(!is.na(PROTEIN_POS))
  
  # Count variants per protein position
  variant_counts <- variants_df %>%
    group_by(PROTEIN_POS) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Set heights for plotting
  domain_height <- max(variant_counts$Count, na.rm = TRUE) * 0.2
  label_offset <- domain_height * 0.6
  
  # Build plot
  plot <- ggplot() +
    
    # Pfam domain rectangles
    geom_rect(
      data = domain_data,
      aes(xmin = Start, xmax = End, ymin = -domain_height, ymax = 0),
      fill = "gray70",
      color = "black",
      alpha = 0.4,
      inherit.aes = FALSE
    ) +
    
    # Domain name labels (angled below boxes)
    geom_text(
      data = domain_data,
      aes(x = (Start + End) / 2, y = -domain_height - label_offset, label = Name),
      size = 3,
      angle = 45,
      hjust = 0.5,
      vjust = 0,
      inherit.aes = FALSE
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
    scale_y_continuous(
      "Number of Variants",
      limits = c(-domain_height - label_offset * 2.5, max(variant_counts$Count) + 1),
      expand = c(0, 0)
    ) +
    ggtitle(paste("gnomAD Missense Variants in", gene, "(Pfam Domains Only)")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
    )
  
  # Save plot
  ggsave("lollipop_gnomAD_GABRG2_Pfam_labeled.png", plot, width = 12, height = 4.5, dpi = 300)
  message("✅ Lollipop plot with Pfam domain names saved as: lollipop_gnomAD_GABRG2_Pfam_labeled.png")
}

###############################################################################
# 5) Run everything
###############################################################################

run_gnomad_lollipop_pipeline()
