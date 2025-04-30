from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
import pandas as pd
from pathlib import Path
import csv
import time

# Set email for NCBI API
Entrez.email = "saracmc21@gmail.com"

# Initialize NCBI Taxa
ncbi = NCBITaxa()

def read_ids(file_path):
    """Read accession IDs from a file."""
    with open(file_path, 'r') as file:
        ids = [line.strip() for line in file if line.strip()]
    return ids

def fetch_taxonomy(accession_id):
    """
    Fetches the taxonomy information for a given accession ID from the NCBI database.

    Args:
        accession_id (str): The accession ID to fetch taxonomy information for.

    Returns:
        dict: A dictionary with taxonomic information, including order, family, genus, species, etc.
        or None if no data is found.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        organism = record.annotations.get("organism", "Unknown")
        if organism == "Unknown":
            return None

        # Get TaxID with ete3
        name_to_taxid = ncbi.get_name_translator([organism])
        if organism not in name_to_taxid:
            return None

        taxid = name_to_taxid[organism][0]
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)

        taxonomy = {"Accession ID": accession_id}
        for taxid in lineage:
            rank = ranks[taxid]
            if rank in ["order", "family", "genus"]:
                taxonomy[rank.capitalize()] = names[taxid]

        taxonomy["Species"] = organism

        print(f"Taxonomy for {accession_id}: {taxonomy}")
        return taxonomy

    except Exception as e:
        print(f"Error fetching ID {accession_id}: {e}")
        return None

def save_to_csv(results, output_file):
    """
    Saves the list of results (taxonomic information) to a CSV file.

    Args:
        results (list): A list of dictionaries containing taxonomic information.
        output_file (str): The path to the output CSV file.
    """
    with open(output_file, 'a', newline='') as csvfile:
        fieldnames = ["Accession ID", "Order", "Family", "Genus", "Species"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # If file is empty, write the header first
        if csvfile.tell() == 0:
            writer.writeheader()

        writer.writerows(results)

if __name__ == "__main__":
    ids = read_ids("results/blast/arthropoda.blast.top3.unique.ids.txt")
    print(f"Fetched {len(ids)} accession IDs.")

    # Define batch size (100 IDs per batch)
    batch_size = 100
    results = []

    for i, accession_id in enumerate(ids, start=1):
        taxonomy_info = fetch_taxonomy(accession_id)
        if taxonomy_info:
            print(f"{accession_id}: {taxonomy_info['Species']}")

            results.append(taxonomy_info)

        # Every batch_size IDs, save to CSV
        if i % batch_size == 0 or i == len(ids):  # Last batch
            save_to_csv(results, "results/blast/arthropoda_taxonomy.csv")
            results = []
        # Pause to avoid rate limit
        time.sleep(0.34)

    print("Process completed.")
