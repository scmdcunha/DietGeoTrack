#!/bin/bash

# Get the absolute path of the script directory
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Define absolute paths for HIT_IDS and FASTA based on the new directory structure
HIT_IDS="$SCRIPT_DIR/../results/vsearch/results/species/species_ids.txt"
FASTA="$SCRIPT_DIR/../data/reference/arthropoda.fasta"

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
