import pandas as pd
import matplotlib.pyplot as plt

# Load the saved CSV
df = pd.read_csv("exon_clinvar_gnomad_ratio_output.csv")

# Compute ratio again in case not present
df["clinvar_gnomad_ratio"] = df["variant_count_clinvar"] / (df["variant_count_gnomad"] + 1e-6)

# Sort for better x-axis
df_sorted = df.sort_values("clinvar_gnomad_ratio", ascending=False).reset_index(drop=True)

# Plot
plt.figure(figsize=(14, 6))
plt.scatter(df_sorted.index, df_sorted["clinvar_gnomad_ratio"], alpha=0.7)
plt.axhline(1.0, color='red', linestyle='--', label="ClinVar = gnomAD")
plt.xlabel("Exons (sorted by ratio)")
plt.ylabel("ClinVar / gnomAD variant ratio")
plt.title("Ratio of ClinVar to gnomAD Variants per Exon")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
