from Bio import SeqIO
import sys

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

with open(output_fasta, 'w') as output_file:
    for record in SeqIO.parse(input_fasta, 'fasta'):
        record.seq = record.seq.upper()
        SeqIO.write(record, output_file, 'fasta')
