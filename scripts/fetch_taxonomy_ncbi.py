from Bio import Entrez, SeqIO
import csv
import time

# Set email for NCBI API
Entrez.email = "saracmc21@gmail.com"

# Function to read accession IDs from a file
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

        # Fetching the taxonomy data from the annotations
        taxonomy = record.annotations.get("taxonomy", [])

        # Ensure taxonomy has all levels (fill with "Unknown" if not available)
        taxonomy_info = {
            "Accession ID": accession_id,
            "Domain": taxonomy[0] if len(taxonomy) > 0 else "Unknown",
            "Phylum": taxonomy[1] if len(taxonomy) > 1 else "Unknown",
            "Class": taxonomy[2] if len(taxonomy) > 2 else "Unknown",
            "Order": taxonomy[3] if len(taxonomy) > 3 else "Unknown",
            "Family": taxonomy[4] if len(taxonomy) > 4 else "Unknown",
            "Genus": taxonomy[5] if len(taxonomy) > 5 else "Unknown",
            "Species": taxonomy[6] if len(taxonomy) > 6 else "Unknown"
            }

        print(f"Taxonomy for {accession_id}: {taxonomy_info}")
        return taxonomy_info

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
        fieldnames = ["Accession ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
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
            results = []  # Clear results for the next batch

        # Pause to avoid rate limit
        time.sleep(0.34)

    print("Process completed.")
