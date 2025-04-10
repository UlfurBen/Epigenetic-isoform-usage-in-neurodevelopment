import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Load the CSV
df = pd.read_csv("EM_genes_exon_clinvar_gnomad_ratio_output.csv")

# Recompute ratio
df["clinvar_gnomad_ratio"] = df["variant_count_clinvar"] / (df["variant_count_gnomad"] + 1e-6)

# Map chromosome values: 'X' → 23, 'Y' → 24
df["exon_chr_numeric"] = df["exon_chr"].replace({'X': 23, 'Y': 24}).astype(int)

# Sort by chromosome and exon start
df_sorted = df.sort_values(by=["exon_chr_numeric", "exon_start"]).reset_index(drop=True)

# Plot
plt.figure(figsize=(18, 6))
plt.scatter(df_sorted.index, df_sorted["clinvar_gnomad_ratio"], alpha=0.7, label="Exons")

# Add reference line at y = 1
plt.axhline(1.0, color='red', linestyle='--', label="ClinVar = gnomAD")

# Use log scale on y-axis
plt.yscale("log")
plt.gca().yaxis.set_major_formatter(ticker.LogFormatter())  # Make ticks show nicely

# Annotate exons with ratio > 10^4
for i, row in df_sorted.iterrows():
    if row["clinvar_gnomad_ratio"] > 1e4:
        label = f"{row['gene']}"
        plt.annotate(
            label,
            (i, row["clinvar_gnomad_ratio"]),
            textcoords="offset points",
            xytext=(0, 5),
            ha='center',
            fontsize=8,
            rotation=45
        )

# Styling
plt.xlabel("Exons (sorted by genomic position)")
plt.ylabel("ClinVar / gnomAD variant ratio (log scale)")
plt.title("ClinVar/gnoMAD Variant Ratios Across Exons (Log Scale, Genomic Order)")
plt.grid(True, which="both", axis='y', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.tight_layout(pad=3.0)


# Save and show
plt.savefig("EM_genes_clinvar_gnomad_ratio_by_genomic_position_log_annotated.png", dpi=300)
plt.show()

