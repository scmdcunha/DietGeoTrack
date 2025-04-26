from Bio import Entrez, SeqIO
import csv
import time

# Set email for NCBI API
Entrez.email = ""

def read_ids(file_path):
    """
    Reads accession IDs from a specified file and returns them as a list.

    Args:
        file_path (str): The path to the file containing the accession IDs.

    Returns:
        list: A list of accession IDs as strings.

    Example:
        ids = read_ids("results/blast/arthropoda.blast.top3.unique.ids.txt")
    """
    with open(file_path, 'r') as file:
        ids = [line.strip() for line in file if line.strip()]
    return ids

def fetch_taxonomy(accession_id):
    """
    Fetches the taxonomy information (organism name) for a given accession ID from the NCBI database.

    Args:
        accession_id (str): The accession ID to fetch taxonomy information for.

    Returns:
        str: The organism name, or "Unknown" if no organism is found.

    Example:
        organism = fetch_taxonomy("AB007981.1")
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        organism = record.annotations.get("organism", "Unknown")
        return organism

    except Exception as e:
        print(f"Error fetching ID {accession_id}: {e}")
        return "Error"

def save_to_csv(results, output_file):
    """
    Saves the list of results (accession ID and organism) to a CSV file.

    Args:
        results (list): A list of dictionaries with "Accession ID" and "Organism".
        output_file (str): The path to the output CSV file.

    Example:
        save_to_csv(results, "results/blast/arthropoda_top3_taxonomy.csv")
    """
    with open(output_file, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["Accession ID", "Organism"])
        writer.writerows(results)

if __name__ == "__main__":
    ids = read_ids("results/blast/arthropoda.blast.top3.unique.ids.txt")
    print(f"Fetched {len(ids)} accession IDs.")

    # Define batch size (100 IDs per batch)
    batch_size = 100
    results = []

    for i, accession_id in enumerate(ids, start=1):
        organism = fetch_taxonomy(accession_id)
        print(f"{accession_id}: {organism}")

        results.append({
            "Accession ID": accession_id,
            "Organism": organism
        })

        # Every batch_size IDs, save to CSV
        if i % batch_size == 0 or i == len(ids):  # Last batch
            save_to_csv(results, "results/blast/arthropoda_top3_taxonomy.csv")
            results = []  # Clear results for the next batch

        # Pause to avoid rate limit
        time.sleep(0.34)

    print("Process completed.")
