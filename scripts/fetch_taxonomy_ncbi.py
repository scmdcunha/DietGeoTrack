from Bio import Entrez
import csv

# Set email for NCBI API
Entrez.email = "saracmc21@gmail.com"

# Function to read accession IDs from a file
def read_ids(file_path):
    with open(file_path, 'r') as file:
        ids = [line.strip() for line in file if line.strip()]
    return ids

if __name__ == "__main__":
    ids = read_ids("results/blast/arthropoda.blast.top3.unique.ids.txt")
    print(ids)
