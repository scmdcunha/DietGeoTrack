
# Diet Geo Track

Diet Geo Track is a bioinformatics pipeline designed to identify arthropods (or other prey) consumed by animals through the analysis of metabarcoding sequence data, such as those obtained from guano samples or similar environmental DNA sources. The main goal is to process and analyze Next-Generation Sequencing (NGS) data to assign taxonomic identities to query sequences and then validate these identifications geographically using occurrence data from GBIF (Global Biodiversity Information Facility). This ensures the taxonomic assignments make ecological and biogeographical sense.

## Project Overview

- **Input:** FASTA sequences obtained by metabarcoding or similar NGS approaches from dietary samples.
- **Reference database:** Arthropod COI gene sequences (or any relevant reference sequences) downloaded from NCBI or other sources.
- **Taxonomic assignment:** Performed by sequence alignment using BLAST, selecting the top hits per query sequence filtered by percent identity thresholds.
- **Geographic validation:** Top hits are compared against GBIF occurrence records to filter improbable identifications based on geographic proximity and occurrence date.
- **Automation:** The pipeline is managed using Snakemake and containerized with Docker and Micromamba to ensure reproducibility and ease of use.

## Tools and Technologies

- Snakemake
- BLAST+
- VSEARCH
- GBIF API
- Docker & Micromamba
- Python (with Biopython, pandas, ete3)
- pyGDAL for geographic calculations

## Requirements

- Python 3.8+
- BLAST+ and VSEARCH installed in your environment
- Docker installed (if you want to use containerization)

## Installation and Setup

### 1. Clone the repository

```bash
git clone https://github.com/scmdcunha/diet-geo-track.git
cd diet-geo-track
```

### 2. Install Docker

If you do not have Docker installed or are unfamiliar with it, please follow the official Docker installation guide:

- [Docker Installation Guide](https://docs.docker.com/get-docker/)

### 3. Build the Docker image

Inside the project directory, build the Docker container with:

```bash
docker build -t diet-geo-track .
```

### 4. Run the Docker container

To start the container and mount the current directory inside it (so you can access files):

```bash
docker run -it --rm -v "$(pwd)":/app diet-geo-track bash
```


### 5. Configure your parameters in `config.yaml`

Edit the `config.yaml` file to set your specific parameters (see below for detailed explanations and an example).

### 6. Run the pipeline with Snakemake

Run the workflow specifying the number of CPU cores and your config file:

```bash
snakemake --cores 4 --configfile config.yaml
```

---

## Configuration Parameters

Before running the pipeline, configure the `config.yaml` file with the following parameters:

- **reference_fasta**: Path to the FASTA file with reference sequences (e.g., arthropod COI sequences).
- **blast_query**: Path to the FASTA file containing your query sequences.
- **blast_threads**: Number of CPU threads to use for BLAST.
- **blast_max_targets**: Maximum number of BLAST hits to retrieve per query.
- **min_identity**: Minimum percent identity (%) threshold to consider a hit valid.
- **top_hits**: Number of top hits (by percent identity) to retain per query.
- **ncbi_email**: Your email address required by NCBI Entrez API.
- **ncbi_api_key**: Optional NCBI API key to increase request limits. Obtain it at: https://www.ncbi.nlm.nih.gov/account/settings/
- **taxonomy_threads**: Number of parallel threads for taxonomy fetching.
- **occurrence_column**: Column name in input files containing species scientific names for GBIF queries.
- **ref_lat** and **ref_lon**: Latitude and longitude of the geographic reference point to filter GBIF occurrences.
- **radius**: Search radius (km) around the reference point for GBIF occurrence search.
- **occ_top_n**: Number of closest GBIF occurrences to consider per species.
- **min_year**: Minimum year of occurrence data to accept.
- **cache_dir**: Directory to cache API responses.
- **occ_threads**: Number of threads for GBIF occurrence queries.
- **w_identity**, **w_distance**, **w_date**: Weights for identity, geographic distance, and occurrence date used in the final scoring. Must sum to 1.

---

## Example `config.yaml`

```yaml
reference_fasta: data/arthropoda.fasta
blast_query: data/queries.fasta

blast_threads: 4
blast_max_targets: 10
min_identity: 97
top_hits: 3

ncbi_email: "your_email@example.com"        # Required by NCBI Entrez API
ncbi_api_key: ""                            # Optional: get your API key at https://www.ncbi.nlm.nih.gov/account/settings/
taxonomy_threads: 4

occurrence_column: Species
ref_lat: 40.3397
ref_lon: -7.6120
radius: 20
occ_top_n: 1
min_year: 2010
cache_dir: cache
occ_threads: 4

w_identity: 0.5
w_distance: 0.3
w_date: 0.2
```

---

## Additional Scripts

- **fetch_taxonomy_ncbi.py**
  Fetches taxonomy information from NCBI for accession numbers.
  Usage example:
  ```bash
  python3 scripts/fetch_taxonomy_ncbi.py --input results/blast/accessions.txt --output results/taxonomy.csv --email your_email@example.com
  ```

- **gbif_occurrences.py**
  Retrieves the closest GBIF occurrences for species names from a CSV file.
  Usage example:
  ```bash
  python3 scripts/gbif_occurrences.py --input_csv data/species.csv --column Species --output_csv results/gbif_occurrences.csv --ref_lat 40.3397 --ref_lon -7.6120 --radius 20
  ```

---

## References and Documentation

- [NCBI Entrez API Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [GBIF API Documentation](https://www.gbif.org/developer/summary)
- [BLAST+ Documentation](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
- [VSEARCH Documentation](https://github.com/torognes/vsearch/wiki)
- [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/)
- [Docker Documentation](https://docs.docker.com/)

---

## License

This project is licensed under the GNU General Public License v3.0 (GPLv3).

---

## Author

Developed by Sara Cristina Marques da Cunha (scmdcunha).
