import pandas as pd
import matplotlib.pyplot as plt

# Load the saved CSV
df = pd.read_csv("whole_genome_exon_clinvar_gnomad_ratio_output.csv")

# Recompute ratio to be sure
df["clinvar_gnomad_ratio"] = df["variant_count_clinvar"] / (df["variant_count_gnomad"] + 1e-6)

# Sort by ratio descending
df_sorted = df.sort_values("clinvar_gnomad_ratio", ascending=False).reset_index(drop=True)

# Plot setup
plt.figure(figsize=(14, 6))
plt.scatter(df_sorted.index, df_sorted["clinvar_gnomad_ratio"], alpha=0.7)

# Add horizontal reference line at y=1
plt.axhline(1.0, color='red', linestyle='--', label="ClinVar = gnomAD")

# Find top 10 exons with ratio > 1
top_exons = df_sorted[df_sorted["clinvar_gnomad_ratio"] > 1].head(10)

# Annotate those top 10 exons
for i, row in top_exons.iterrows():
    plt.annotate(
        row["ensembl_exon_id"],
        (i, row["clinvar_gnomad_ratio"]),
        textcoords="offset points",
        xytext=(0, 5),
        ha='center',
        fontsize=8,
        rotation=45
    )

# Labels and styling
plt.xlabel("Exons (sorted by ratio)")
plt.ylabel("ClinVar / gnomAD variant ratio")
plt.title("Top Exons with Elevated ClinVar/gnomAD Variant Ratios")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save and show
plt.savefig("whole_genome_clinvar_gnomad_ratio_top10.png", dpi=300)
plt.show()
