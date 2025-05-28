from requests import api
from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
from pathlib import Path
import pandas as pd
import csv
import argparse
from multiprocessing import Pool, Manager, cpu_count
from threading import Lock
from tqdm import tqdm
import sys
from urllib.error import HTTPError
import time

# Initialize NCBI Taxa
ncbi = NCBITaxa()
csv_lock = Lock()
fail_lock = Lock()

def read_sseqids_from(file_path):
    """Read accession IDs from a resulting BLAST/VSEARCH output file."""
    df = pd.read_csv(file_path, sep='\t')
    return df['sseqid'].drop_duplicates().tolist()

def fetch_taxonomy_wrapper(args):
    """
    Fetches taxonomy for a given accession ID from NCBI and returns a dictionary.
    """
    accession_id, email, api_key = args
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        time.sleep(0.12)  # ~8.3 requests/second

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
        print(f"[ERROR] Accession: {accession_id} | {e}")
        time.sleep(0.5)
        return {"error": True, "accession_id": accession_id, "message": str(e)}

def save_to_csv(results, output_file):
    """Saves list of taxonomy dicts to CSV."""
    fieldnames = ["Accession ID", "Order", "Family", "Genus", "Species"]
    with open(output_file, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if csvfile.tell() == 0:
            writer.writeheader()
        writer.writerows(results)

def main(input_file, output_file, email, api_key=None, threads=cpu_count()):
    ids = read_sseqids_from(input_file)
    print(f"Fetched {len(ids)} accession IDs.")

    if api_key and threads > 10:
        print("Reducing threads to 10 (NCBI limit with API key).")
        threads = 10
    failed_ids_path = Path(output_file).with_name("failed_ids.txt")
    results = []
    failed = []

    with Pool(processes=threads) as pool:
        args = [(acc_id, email, api_key) for acc_id in ids]
        for result in tqdm(pool.imap_unordered(fetch_taxonomy_wrapper, args), total=len(args), desc="Fetching taxonomy", file=sys.stdout):
            if result is None or "error" in result:
                failed.append(result.get("accession_id", "UNKNOWN") if result else "UNKNOWN")
            else:
                results.append(result)

            if len(results) >= 100:
                save_to_csv(results, output_file)
                results = []

    # Save remaining
    if results:
        save_to_csv(results, output_file)

    if failed:
        with open(failed_ids_path, 'w') as fail_log:
            for fid in failed:
                fail_log.write(fid + '\n')
    print("Process completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch taxonomy from NCBI for accession IDs.")
    parser.add_argument("--input", required=True, help="Path to input file with accession IDs")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--api_key", required=False, help="API key for NCBI Entrez")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of processes (default=max cores)")
    args = parser.parse_args()

    main(args.input, args.output, args.email, args.api_key, args.threads)
