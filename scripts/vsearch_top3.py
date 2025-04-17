#!/usr/bin/env python3

# Select top 3 VSEARCH hits per query based on percent identity

import os
import pandas as pd

INPUT = "results/vsearch/results.tsv"
OUTPUT_DIR = "results/vsearch/filtered"
OUTPUT = os.path.join(OUTPUT_DIR, "top3_hits.tsv")

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Read VSEARCH results (no header)
colnames = ["query", "target", "pident", "length", "mismatch", "gapopen", "qstart", "qend" , "qend", "sstart", "send", "evalue", "bitscore"
]

df = pd.read_csv(INPUT, sep='\t', header=None, names=colnames)

# Sort by query and percent identity (descending)
df_sorted = df.sort_values(['query', 'pident'], ascending=[True, False])

# Select top 3 hits per query
df_top3 = df_sorted.groupby('query').head(3)

# Save to output file
df_top3.to_csv(OUTPUT, sep= '\t', index=False, header=False)
