from pathlib import Path

SAMPLES = "data/guano_samples.fasta"
REFERENCE = "data/arthropoda.fasta"
BLAST_DB_PREFIX = "data/arthropoda.blastdb"

rule all:
    input:
        "results/blast/arthropoda.blast.tsv",
        "results/blast/arthropoda.blast.top3.tsv",
        "results/blast/arthropoda.blast.top3.unique.ids.txt",
        "results/vsearch/arthropoda.vsearch.tsv"


rule make_blast_db:
    input:
        fasta=REFERENCE
    output:
        touch("data/.blastdb_created")
    shell:
        """
        makeblastdb -in {input.fasta} -dbtype nucl -out {BLAST_DB_PREFIX}
        touch {output}
        """

rule run_blast:
    input:
        query=SAMPLES,
        db_flag="data/.blastdb_created"
    output:
        "results/blast/arthropoda.blast.tsv"
    params:
        db_prefix=BLAST_DB_PREFIX,
        perc_identity=95
    threads: 4
    shell:
        """
        blastn \
            -query {input.query} \
            -db {params.db_prefix} \
            -out {output} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -perc_identity {params.perc_identity} \
            -num_threads {threads}
        """
rule filter_top3_blast_hits:
    input:
        "results/blast/arthropoda.blast.tsv"
    output:
        "results/blast/arthropoda.blast.top3.tsv"
    shell:
        """
        python3 scripts/blast_top3_hits.py {input} {output}
        """
rule extract_unique_ids:
    input:
        "results/blast/arthropoda.blast.top3.tsv"
    output:
        "results/blast/arthropoda.blast.top3.unique.ids.txt"
    shell:
        """
            cut -f2 {input} | sort | uniq > {output}
        """
rule run_vsearch:
    input:
        query=SAMPLES,
        db=REFERENCE
    output:
        "results/vsearch/arthropoda.vsearch.tsv"
    params:
        identity=0.95
    threads: 4
    shell:
        """
        vsearch --usearch_global {input.query} \
             --db {input.db} \
             --id {params.identity} \
             --threads {threads} \
             --userout {output} \
             --userfields query+target+id+length+mismatch+gapopen+sstart+send+bitscore
        """
