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
        db_nin="data/arthropoda.blastdb.nin",
        db_nsq="data/arthropoda.blastdb.nsq",
        db_nhr="data/arthropoda.blastdb.nhr",
        db_ndb="data/arthropoda.blastdb.ndb",
        db_njs="data/arthropoda.blastdb.njs",
        db_ntf="data/arthropoda.blastdb.ntf",
        db_not="data/arthropoda.blastdb.not",
        db_nto="data/arthropoda.blastdb.nto"
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
