from Bio import Entrez, SeqIO
from ete3 import NCBITaxa
from pathlib import Path
import pandas as pd
import csv
import argparse
from multiprocessing import Pool, cpu_count
from threading import Lock
from tqdm import tqdm
import sys
import time

# Initialize NCBI Taxa
ncbi = NCBITaxa()
csv_lock = Lock()
fail_lock = Lock()

def read_accession_query_pairs(file_path):
    """Read qseqid and sseqid pairs from BLAST output file."""
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                  "qend", "sstart", "send", "evalue", "bitscore"]
    pairs = df.drop_duplicates(subset=["qseqid", "sseqid"])[["qseqid", "sseqid"]]
    return pairs.values.tolist()  # returns list of [qseqid, sseqid]


def fetch_taxonomy_wrapper(args):
    """
    Fetches taxonomy for a given accession ID from NCBI and returns a dictionary.
    """
    query_id, accession_id, email, api_key = args
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

        taxonomy = {"Query ID": query_id, "Accession ID": accession_id}
        # Add ranks if they exist
        for taxid in lineage:
            rank = ranks.get(taxid)
            if rank in ["order", "family", "genus"]:
                taxonomy[rank.capitalize()] = names.get(taxid, "NA")

        taxonomy["Species"] = organism
        return taxonomy

    except Exception as e:
        print(f"[ERROR] Query: {query_id} | Accession: {accession_id} | {e}")
        time.sleep(0.5)
        return {"error": True, "query_id": query_id, "accession_id": accession_id, "message": str(e)}

def save_to_csv(results, output_file):
    """Saves list of taxonomy dicts to CSV."""
    fieldnames = ["Query ID", "Accession ID", "Order", "Family", "Genus", "Species"]
    with open(output_file, 'a', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        if csvfile.tell() == 0:
            writer.writeheader()
        writer.writerows(results)


def main(input_file, output_file, email, api_key=None, threads=cpu_count()):
    pairs = read_accession_query_pairs(input_file)
    print(f"Fetched {len(pairs)} (query, accession) pairs.")

    if api_key and threads > 10:
        print("Reducing threads to 10 (NCBI limit with API key).")
        threads = 10
    failed_ids_path = Path(output_file).with_name("failed_ids.txt")
    results = []
    failed = []

    with Pool(processes=threads) as pool:
        args = [(query_id, accession_id, email, api_key) for query_id, accession_id in pairs]
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
