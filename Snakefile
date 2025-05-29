import yaml
from pathlib import Path

configfile: "config.yaml"

SAMPLES = config["samples"]
REFERENCE = config["reference"]
BLAST_DB_PREFIX = config["blast_db_prefix"]
RESULTS = config["results_dir"]
CACHE_DIR = config["cache_dir"]
PERC_IDENTITY = config["perc_identity"]
ref_lat = config["ref_lat"]
ref_lon = config["ref_lon"]

rule all:
    input:
        f"{RESULTS}/arthropoda.blast.tsv",
        f"{RESULTS}/arthropoda.blast.top3.tsv",
        f"{RESULTS}/species_without_occurrence.csv"

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
        f"{RESULTS}/arthropoda.blast.tsv"
    params:
        db_prefix=BLAST_DB_PREFIX,
        perc_identity=PERC_IDENTITY
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
        f"{RESULTS}/arthropoda.blast.tsv"
    output:
        f"{RESULTS}/arthropoda.blast.top3.tsv"
    shell:
        """
        python3 scripts/blast_top3_hits.py {input} {output}
        """

rule fetch_taxonomy:
    input:
        f"{RESULTS}/arthropoda.blast.top3.tsv"
    output:
        f"{RESULTS}/arthropoda.taxonomy.csv"
    shell:
        """
        python3 scripts/fetch_taxonomy_ncbi.py \
        --input {input} \
        --output {output} \
        --email {config[email]} \
        --api_key {config[api_key]} \
        --threads {config[threads]}
        """

rule fetch_occurrences:
    input:
        f"{RESULTS}/arthropoda.taxonomy.csv"
    output:
        f"{RESULTS}/arthropoda.occurrences.csv",
        f"{RESULTS}/species_without_occurrence.csv"
    shell:
        """
        python3 scripts/gbif_closest_occurrence.py \
        --input_csv {input} \
        --column Species \
        --output_csv {output[0]} \
        --no_occurrences_csv {output[1]} \
        --ref_lat {config[ref_lat]} \
        --ref_lon {config[ref_lon]} \
        --radius {config[radius]} \
        --cache_dir {CACHE_DIR} \
        --threads {config[threads]}
        """
