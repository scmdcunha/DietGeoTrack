import pandas as pd
from pathlib import Path

# Paths
blast_file = Path("results/blast/arthropoda.blast.top3.tsv")
taxonomy_file = Path("results/blast/ncbi_taxonomy_lookup.csv")
output_file = Path("results/blast/arthropoda_top3_taxonomy.csv")

# Read BLAST file
blast_df = pd.read_csv(blast_file, sep="\t")

# Read taxonomy file
taxonomy_df = pd.read_csv(taxonomy_file)

# Check required columns
required_blast_cols = {"qseqid", "sseqid", "pident"}
required_tax_cols = {"Accession ID", "Order", "Family", "Genus", "Species"}

if not required_blast_cols.issubset(blast_df.columns):
    raise ValueError(f"Missing columns in BLAST file: {required_blast_cols - set(blast_df.columns)}")

if not required_tax_cols.issubset(taxonomy_df.columns):
    raise ValueError(f"Missing columns in taxonomy file: {required_tax_cols - set(taxonomy_df.columns)}")

# Merge based on accession number
merged_df = pd.merge(
    blast_df,
    taxonomy_df,
    left_on="sseqid",
    right_on="Accession ID",
    how="left"
)

# Select and rename/reorder columns
output_df = merged_df[[
    "qseqid", "sseqid", "pident", "Order", "Family", "Genus", "Species"
]]

# Sort by query and identity for easier inspection
output_df = output_df.sort_values(by=["qseqid", "pident"], ascending=[True, False])

# Save output
output_df.to_csv(output_file, index=False)
print(f"Table saved as {output_file}")
