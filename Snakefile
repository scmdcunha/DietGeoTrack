from pathlib import Path

SAMPLES = "data/guano_samples.fasta"
REFERENCE = "data/arthropoda.fasta"
BLAST_DB = "data/arthropoda.blastdb"

rule all:
    input:
        "results/blast/arthropoda.blast.tsv"

rule make_blast_db:
    input:
        fasta=REFERENCE
    output:
        expand("data/arthropoda.blastdb.{ext}", ext=[
            "nin", "nsq", "nhr", "ndb", "njs", "ntf", "not", "nto"
        ])
    params:
        db_out=BLAST_DB
    shell:
        """
        makeblastdb -in {input.fasta} -dbtype nucl -out {params.db_out}
        """

rule run_blast:
    input:
        query=SAMPLES,
        db=BLAST_DB
    output:
        "results/blast/arthropoda.blast.tsv"
    params:
        perc_identity=95 # percent identity threshold
    threads: 4
    shell:
        """
        blastn \
            -query {input.query} \
            -db {input.db} \
            -out {output} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -perc_identity {params.perc_identity} \
            -num_threads {threads}
        """
