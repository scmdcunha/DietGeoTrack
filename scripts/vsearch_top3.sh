#!/bin/bash

# Paths
QUERY="data/samples/guano_samples.fasta"
DB="data/reference/arthropoda.fasta"
OUT="results/vsearch/results/vsearch_top3_hits.tsv"

# Execute vsearch
vsearch --usearch_global "$QUERY" \
        --db "$DB" \
        --id 0.95 \
        --maxhits 3 \
        --blast6out "$OUT" \
        --threads 4

echo "vsearch completed.: result saved in: \$OUT"
