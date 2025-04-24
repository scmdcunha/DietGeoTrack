from pathlib import Path

SAMPLES = "data/guano_samples.fasta"
REFERENCE = "data/arthropoda.fasta"
BLAST_DB = "data/arthropoda.blastdb"

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
