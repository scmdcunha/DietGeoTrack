import os

# This Snakefile defines the workflow for metabarcoding taxonomic validation.
# It performs the following main steps:
# 1. Creates necessary output directories.
# 2. Builds a local BLAST nucleotide database from a reference FASTA file.
# 3. Runs BLASTn of query sequences against the local database, filtering hits by minimum identity.
# 4. Selects the top N hits per query based on bit score.
# 5. Fetches taxonomic information for the top BLAST hits using an external script.
# 6. Retrieves species occurrence data from GBIF for the identified taxa.
# 7. Calculates a combined validation score integrating BLAST identity, geographic distance, and occurrence dates.
#
# The configuration parameters and paths are read from config.yaml.

configfile: "config.yaml"

rule all:
    input:
        # Final output file containing the validation scores
        "results/scores/final_scores.csv"

rule create_dirs:
    run:
        # Create directories for results if they don't exist
        os.makedirs("results/blast", exist_ok=True)
        os.makedirs("results/db", exist_ok=True)
        os.makedirs("results/taxonomy", exist_ok=True)
        os.makedirs("results/occurrences", exist_ok=True)
        os.makedirs("results/scores", exist_ok=True)

rule makeblastdb:
    input:
        fasta=config["reference_fasta"]  # Reference FASTA for local BLAST database
    output:
        db_nsq="results/db/arthropoda_local_db.nsq"  # Indicator file that BLAST DB is created
    params:
        db_prefix="results/db/arthropoda_local_db"  # Prefix for BLAST DB files
    shell:
        # Create a nucleotide BLAST database from the reference FASTA
        "makeblastdb -in {input.fasta} -dbtype nucl -out {params.db_prefix}"

rule run_blast:
    input:
        query=config["blast_query"],                 # Query sequences to BLAST
        db_nsq="results/db/arthropoda_local_db.nsq" # BLAST DB index file
    output:
        "results/blast/blast_output.tsv"  # BLAST results in tabular format
    params:
        db="results/db/arthropoda_local_db",
        identity=config["min_identity"],          # Minimum percent identity threshold for filtering
        threads=config["blast_threads"],          # Number of threads for BLAST
        max_targets=config["blast_max_targets"]   # Max number of target sequences per query
    shell:
        # Run BLASTn with filtering on minimum identity using awk on column 3
        "blastn -query {input.query} -db {params.db} "
        "-outfmt 6 -max_target_seqs {params.max_targets} "
        "-num_threads {params.threads} -evalue 1e-5 "
        "| awk '$3 >= {params.identity}' > {output}"

rule top_hits:
    input:
        "results/blast/blast_output.tsv"
    output:
        "results/blast/top_hits.tsv"
    params:
        top_n=config["top_hits"]  # Number of top hits to keep per query
    shell:
        # Select top N hits per query with highest percent identity (column 3)
        "awk '{{print $1}}' {input} | sort -u | while read id; do "
        "grep -w \"$id\" {input} | sort -k3,3nr | head -n {params.top_n}; done > {output}"

rule fetch_taxonomy:
    input:
        "results/blast/top_hits.tsv"
    output:
        "results/taxonomy/taxonomy.csv"
    params:
        script="scripts/fetch_taxonomy_ncbi.py",
        email=config["ncbi_email"],
        api_key_arg=f"--api_key {config['ncbi_api_key']}" if config.get("ncbi_api_key", "") else "",
        threads=config["taxonomy_threads"]
    shell:
        # Run the taxonomy fetching script to retrieve taxonomic info for BLAST hits
        "micromamba run -n metabarcoding python {params.script} --input {input} --output {output} "
        "--email {params.email} {params.api_key_arg} --threads {params.threads}"

rule fetch_occurrences:
    input:
        "results/taxonomy/taxonomy.csv"
    output:
        gbif="results/occurrences/occurrences.csv",
        no_occ="results/occurrences/no_occurrences.csv"
    params:
        script="scripts/gbif_occurrences.py",
        column=config["occurrence_column"],
        lat=config["ref_lat"],
        lon=config["ref_lon"],
        radius=config["radius"],
        top_n=config["occ_top_n"],
        min_year_arg=f"--min_year {config['min_year']}" if config.get("min_year", "") else "",
        cache=config["cache_dir"],
        threads=config["occ_threads"]
    shell:
        # Run occurrences fetching script querying GBIF API with parameters
        "micromamba run -n metabarcoding python {params.script} --input_csv {input} --column {params.column} "
        "--output_csv {output.gbif} --no_occurrences_csv {output.no_occ} "
        "--ref_lat {params.lat} --ref_lon {params.lon} --radius {params.radius} "
        "--top_n {params.top_n} {params.min_year_arg} "
        "--cache_dir {params.cache} --threads {params.threads}"

rule calculate_score:
    input:
        blast="results/blast/top_hits.tsv",
        taxonomy="results/taxonomy/taxonomy.csv",
        gbif="results/occurrences/occurrences.csv"
    output:
        scores="results/scores/final_scores.csv",
        missing="results/scores/missing_data.csv"
    params:
        script="scripts/calculate_score.py",
        w_identity=config["w_identity"],
        w_distance=config["w_distance"],
        w_date=config["w_date"],
    shell:
        # Calculate combined score integrating identity, distance, and date using the scoring script
        """
        micromamba run -n metabarcoding python {params.script} \
        --blast {input.blast} \
        --taxonomy {input.taxonomy} \
        --gbif {input.gbif} \
        --output {output.scores} \
        --missing {output.missing} \
        --w_identity {params.w_identity} \
        --w_distance {params.w_distance} \
        --w_date {params.w_date} \
        """
