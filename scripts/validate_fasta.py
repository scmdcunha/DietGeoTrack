# scripts/validate_fasta.py

import sys
from Bio import SeqIO

def validate_fasta(fasta_path):
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
        if not records:
            print(f"Error: No sequences found in FASTA file {fasta_path}")
            sys.exit(1)
    except Exception as e:
        print(f"Error: Invalid FASTA format in {fasta_path}\n{e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python validate_fasta.py <fasta_file>")
        sys.exit(1)

    validate_fasta(sys.argv[1])
