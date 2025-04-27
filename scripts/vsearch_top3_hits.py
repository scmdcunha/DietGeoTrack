#!/usr/bin/env python3

# Select top 3 VSEARCH hits per query based on percent identity

import os
import pandas as pd
from pathlib import Path

input_file = Path("results/vsearch/arthropoda.vsearch.tsv")
output_file = Path("results/vsearch/arthropoda.vsearch.top3.tsv")

# Read VSEARCH results (no header)
columns = ["query", "target", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

df = pd.read_csv(input_file, sep='\t', header=None, names=columns)

# Sort by query and percent identity (descending)
df_sorted = df.sort_values(['query', 'pident'], ascending=[True, False])

# Select top 3 hits per query
df_top3 = df_sorted.groupby('query').head(3)

# Ensure the output directory exists
output_file.parent.mkdir(parents=True, exist_ok=True)

# Save to output file
df_top3.to_csv(output_file, sep='\t', index=False, header=False)
