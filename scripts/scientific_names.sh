#!/bin/bash

# Paths
HIT_IDS="/results/vsearch/results/species"
FASTA="/data/reference/arthropoda.fasta"

# Output file
OUTPUT="species_names.txt"

# Clean previous output file
> "$OUTPUT"

# For each ID, search the FASTA file and extract the scientific name
while read -r id; do
    # Search the line with matching header
    line=$(grep -m 1 "^>$id " "$FASTA")

    if [[ $line ]]; then
        # Extract the 2nd and 3rd words (scientific name)
        species=$(echo "$line" | awk '{print $2, $3}')
        echo "$species" >> "$OUTPUT"
    fi
done < "$HIT_IDS"

# Remove duplicates and sort
sort "$OUTPUT" | uniq > tmp && mv tmp "$OUTPUT"

echo "Scientific names extracted to $OUTPUT"
