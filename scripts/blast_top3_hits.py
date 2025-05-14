import pandas as pd
from pathlib import Path

# Paths to input and output files
input_file = Path("results/blast/arthropoda.blast.tsv")
output_file = Path("results/blast/arthropoda.blast.top3.tsv")

# Column names according to BLAST outfmt 6
columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def main():
    # Read the BLAST results
    df = pd.read_csv(input_file, sep="\t", names=columns)

    # Sort by query ID and percent identity (descending),
    # then keep top 3 hits per query
    top3 = (
        df.sort_values(by=["qseqid", "pident"], ascending=[True, False])
        .groupby("qseqid")
        .head(3)
    )

    # Save the result
    output_file.parent.mkdir(parents=True, exist_ok=True)
    top3.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
