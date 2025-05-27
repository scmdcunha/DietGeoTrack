from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
from pathlib import Path
import pandas as pd
import csv
import time
import argparse

# Initialize NCBI Taxa
ncbi = NCBITaxa()

def read_sseqids_from(file_path):
    """Read accession IDs from a resulting BLAST/VSEARCH output file."""
    df = pd.read_csv(file_path, sep='\t')
    return df['sseqid'].drop_duplicates().tolist()

def fetch_taxonomy(accession_id):
    """
    Fetches taxonomy for a given accession ID from NCBI and returns a dictionary.
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
        # Add ranks if they exist
        for taxid in lineage:
            rank = ranks.get(taxid)
            if rank in ["order", "family", "genus"]:
                taxonomy[rank.capitalize()] = names.get(taxid, "NA")

        taxonomy["Species"] = organism
        return taxonomy

    except Exception as e:
        print(f"[{accession_id}] Error: {type(e).__name__} - {e}")
        return None

def save_to_csv(results, output_file):
    """Saves list of taxonomy dicts to CSV."""
    with open(output_file, 'a', newline='') as csvfile:
        fieldnames = ["Accession ID", "Order", "Family", "Genus", "Species"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if csvfile.tell() == 0:
            writer.writeheader()

        writer.writerows(results)

def main(input_file, output_file):
    ids = read_sseqids_from(input_file)
    print(f"Fetched {len(ids)} accession IDs.")

    batch_size = 100
    results = []
    failed_ids_path = Path(output_file).with_name("failed_ids.txt")

    for i, accession_id in enumerate(ids, start=1):
        taxonomy_info = fetch_taxonomy(accession_id)
        if taxonomy_info:
            print(f"{accession_id}: {taxonomy_info.get('Species', 'Unknown')}")
            results.append(taxonomy_info)
        else:
            with open(failed_ids_path, 'a') as fail_log:
                fail_log.write(accession_id + "\n")

        if i % batch_size == 0 or i == len(ids):
            save_to_csv(results, output_file)
            results = []

        time.sleep(0.34)  # prevent rate-limiting

    print("Process completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch taxonomy from NCBI for accession IDs.")
    parser.add_argument("--input", required=True, help="Path to input file with accession IDs")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    args = parser.parse_args()
    Entrez.email = args.email
    main(args.input, args.output)
