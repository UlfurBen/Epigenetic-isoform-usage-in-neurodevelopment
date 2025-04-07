import hail as hl
import pandas as pd
import os

# Initialize Hail with local threads
hl.init(master="local[60]")

# Load your list of gene symbols
genes = pd.read_csv("gene_symbols_test.csv", header=None).values.flatten().tolist()

# Path to gnomAD Hail Table
gnomad_path = "/proj/hpcdata/Mimir/kimberl/gnomad_bigData/gnomad.exomes.v4.1.sites.ht"

# Read gnomAD HT
ht = hl.read_table(gnomad_path)

# Filter to missense variants
ht = ht.filter(ht.vep.most_severe_consequence == "missense_variant")

# Loop through each gene and export
for gene in genes:
  try:
  print(f"Processing {gene}...")

gene_filtered = ht.filter(
  hl.any(lambda tx: tx.gene_symbol == gene, ht.vep.transcript_consequences)
)

if gene_filtered.count() == 0:
  print(f"No variants found for {gene}. Skipping.")
continue

# Explode transcript consequences to flatten nested info
exploded = gene_filtered.explode(gene_filtered.vep.transcript_consequences)

# Select relevant fields
flattened = exploded.select(
  Chromosome=exploded.locus.contig,
  Position=hl.int64(exploded.locus.position),
  Reference=exploded.alleles[0],
  Alternate=exploded.alleles[1],
  transcript_id=exploded.vep.transcript_consequences.transcript_id,
  hgvsc=exploded.vep.transcript_consequences.hgvsc,
  hgvsp=exploded.vep.transcript_consequences.hgvsp,
  rsids=exploded.rsid,
  consequence=exploded.vep.transcript_consequences.consequence_terms
)

# Export to compressed TSV file
out_path = f"gnomAD_{gene}.tsv.bgz"
flattened.export(out_path)
print(f"Saved: {out_path}")

except Exception as e:
  print(f"Error with {gene}: {e}")
