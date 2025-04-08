# Load required libraries
library(biomaRt)
library(dplyr)

# Connect to Ensembl BioMart using the "useast" mirror
ensembl_mart <- useEnsembl(biomart = "ensembl", 
                           dataset = "hsapiens_gene_ensembl", 
                           mirror = "useast")

# List attributes containing "external_gene_name" for verification
# all_attrs <- listAttributes(ensembl_mart)
# gene_attrs <- all_attrs[grep("external_gene_name", all_attrs$name, ignore.case = TRUE), ]
# print(gene_attrs)

# Define the list of target genes (using their external gene names)
target_genes <- c("KMT2A")

# Retrieve Ensembl Gene IDs for target genes using external_gene_name
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = target_genes,
  mart = ensembl_mart
)

# Check the retrieved gene_info column names
print(names(gene_info))

# Stop if no genes were found
if (nrow(gene_info) == 0) {
  stop("No matching genes found in Ensembl. Check spelling or dataset availability.")
}

# Retrieve transcript data (including transcript length)
transcript_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

# Identify canonical transcripts (longest transcript per gene)
canonical_transcripts <- transcript_data %>%
  group_by(ensembl_gene_id) %>%
  slice_max(transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(ensembl_transcript_id)

# Retrieve exon data for all target genes
exon_data <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
                 "exon_chrom_start", "exon_chrom_end", "chromosome_name"),
  filters = "ensembl_gene_id",
  values = gene_info$ensembl_gene_id,
  mart = ensembl_mart
)

# Merge gene info with exon data to include external gene names
exons_joined <- merge(exon_data, gene_info, by = "ensembl_gene_id")

# Print column names to check what we have
print(names(exons_joined))

# Rename columns using base R
colnames(exons_joined)[colnames(exons_joined) == "external_gene_name"] <- "gene"
colnames(exons_joined)[colnames(exons_joined) == "chromosome_name"] <- "exon_chr"
colnames(exons_joined)[colnames(exons_joined) == "exon_chrom_start"] <- "exon_start"
colnames(exons_joined)[colnames(exons_joined) == "exon_chrom_end"] <- "exon_end"

# (Optional) Verify the new names
print(names(exons_joined))

# Convert exon start and end positions to numeric
exons_joined$exon_start <- as.numeric(exons_joined$exon_start)
exons_joined$exon_end   <- as.numeric(exons_joined$exon_end)

# Group by exon and count how many distinct transcripts each exon appears in.
# Also, annotate isoform type based on whether any of the transcripts is canonical.
unique_exons <- exons_joined %>%
  group_by(ensembl_exon_id, gene, exon_chr, exon_start, exon_end) %>%
  summarize(
    n_transcripts = n_distinct(ensembl_transcript_id),
    isoform_type = ifelse(any(ensembl_transcript_id %in% canonical_transcripts$ensembl_transcript_id),
                          "canonical", "non_canonical"),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_transcripts == 1)  # Keep only exons unique to a single isoform

# Save all unique exons (both canonical and non-canonical)
write.csv(unique_exons, "all_exons_canonical_and_noncanonical.csv", row.names = FALSE)

# Save only canonical unique exons
canonical_exons_only <- unique_exons %>% filter(isoform_type == "canonical")
write.csv(canonical_exons_only, "canonical_unique_exons.csv", row.names = FALSE)

cat("✔ All unique exons saved to 'all_exons_canonical_and_noncanonical.csv'\n")
cat("✔ Canonical unique exons saved to 'canonical_unique_exons.csv'\n")
