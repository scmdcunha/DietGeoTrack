from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
from pathlib import Path
import pandas as pd
import csv
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from tqdm import tqdm
import sys

# Initialize NCBI Taxa
ncbi = NCBITaxa()
csv_lock = Lock()
fail_lock = Lock()

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
        return {"error": True, "accession_id": accession_id, "message": str(e)}

def save_to_csv(results, output_file):
    """Saves list of taxonomy dicts to CSV."""
    with csv_lock, open(output_file, 'a', newline='') as csvfile:
        fieldnames = ["Accession ID", "Order", "Family", "Genus", "Species"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if csvfile.tell() == 0:
            writer.writeheader()
        writer.writerows(results)

def main(input_file, output_file, threads):
    ids = read_sseqids_from(input_file)
    print(f"Fetched {len(ids)} accession IDs.")

    failed_ids_path = Path(output_file).with_name("failed_ids.txt")
    results = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(fetch_taxonomy, acc_id): acc_id for acc_id in ids}

        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching taxonomy", file=sys.stdout):

            result = future.result()
            if result is None:
                with fail_lock, open(failed_ids_path, 'a') as fail_log:
                    fail_log.write(futures[future] + "\n")
            elif "error" in result:
                with fail_lock, open(failed_ids_path, 'a') as fail_log:
                    fail_log.write(result.get("accession_id", "UNKNOWN") + "\n")
            else:
                results.append(result)

            # Save in batches of 100
            if len(results) >= 100:
                save_to_csv(results, output_file)
                results = []

    # Save remaining
    if results:
        save_to_csv(results, output_file)

    print("Process completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch taxonomy from NCBI for accession IDs.")
    parser.add_argument("--input", required=True, help="Path to input file with accession IDs")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default=4)")
    args = parser.parse_args()

    Entrez.email = args.email
    main(args.input, args.output, args.threads)
