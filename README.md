# Diet Geo Track

**Diet Geo Track** is a bioinformatics pipeline designed to identify arthropods consumed by bats through the analysis of guano samples collected in Serra da Estrela, Portugal. The project focuses on processing and analyzing sequence data obtained via Next-Generation Sequencing (NGS), with the goal of validating taxonomic identifications geographically using GBIF occurrence data. This ensures that the identified species make ecological and biogeographical sense.

> ⚠️ **Note:** This project is in its early development phase. The current focus is on testing and comparing results from BLAST and VSEARCH. The Dockerfile and Snakefile are still under construction.

---

## Project Overview

- **Input:** FASTA sequences (~5700) obtained through NGS of bat guano samples. These sequences were provided externally and likely generated via metabarcoding (the lab work was not performed in this project).
- **Reference database:** Arthropod COI gene sequences downloaded from the NCBI database using:
    ```
    arthropoda[organism] AND COI[gene]
    ```
- **Taxonomic assignment:** Performed with `vsearch` and `BLAST`, selecting the **top 3 hits per query** with **identity thresholds ≥ 95%**. These are then filtered to assign taxonomic ranks:
    - ≥ 99% identity → species
    - ≥ 97% → genus
    - ≥ 95% → family
    - ≥ 90% → order
- **Geographic validation:** Top hits are compared against GBIF occurrence records to filter out improbable identifications based on geography.
- **Automation goal:** The pipeline is being developed with Snakemake and containerized via Docker and Micromamba, to ensure reproducibility and ease of use for biologists and researchers.

---

## Tools and Technologies

- [Snakemake](https://snakemake.readthedocs.io/)
- [vsearch](https://github.com/torognes/vsearch)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- GBIF API
- [Docker](https://www.docker.com/)
- [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
- [Python](https://www.python.org/) (with pandas, Biopython, ete3)

---

## Usage

### 1. Requirements

- Python 3.8+
- [Biopython](https://biopython.org/)
- [ete3](http://etetoolkit.org/)
- [pandas](https://pandas.pydata.org/)
- BLAST+ and VSEARCH installed (for the main pipeline)

### 2. Fetch Taxonomy for Accession Numbers

You can use the provided script to fetch taxonomy information from NCBI for a list of accession numbers.

**Script:** `fetch_taxonomy_ncbi.py`

**Arguments:**
- `--input` : Path to a text file with one accession number per line.
- `--output` : Path to the output CSV file.
- `--email` : Your email address (required by NCBI Entrez).

**Example usage:**
    ```
    python3 scripts/fetch_taxonomy.py --input results/blast/arthropoda.blast.top3.unique.ids.txt
    --output results/blast/ncbi_taxonomy_lookup.csv --email your_email@domain.com
    ```

**Output:**
- A CSV file with columns: Accession ID, Order, Family, Genus, Species.
- A `failed_ids.txt` file (in the same folder as the output) with accession numbers that could not be resolved.

### 3. Fetch Closest GBIF Occurrences for Species

You can use the provided script to retrieve the closest GBIF occurrence for each species in a CSV file, based on a defined geographic reference point (e.g., location of sample collection).

**Script:** `gbif_occurrences.py`

**Arguments:**
- `--input_csv` : Path to the input CSV file containing species names.
- `--column` : Name of the column containing the scientific names.
- `--output_csv` : Path to the output CSV file where the closest occurrences will be saved.
- `--no_occurrences_csv` : Path to the CSV file for species without nearby occurrences (default: `no_occurrences.csv`).
- `--ref_lat` : Latitude of the reference point (default: `40.3397` – Serra da Estrela).
- `--ref_lon` : Longitude of the reference point (default: `-7.6120`).
- `--radius` : Search radius (in kilometers) around the reference point (default: `20`).
- `--cache_dir` : Directory to cache GBIF API responses (default: `cache`).
- `--threads` : Number of parallel threads to use (default: `5`).

**Example usage:**

    ```
    python3 scripts/gbif_closest_occurrence.py
    --input_csv data/species_list.csv
    --column species
    --output_csv results/gbif/closest_occurrences.csv
    --no_occurrences_csv results/gbif/species_without_occurrence.csv
    --ref_lat 40.3397
    --ref_lon -7.6120
    --radius 50
    --cache_dir cache
    --threads 6
    ```

**Output:**
- A CSV file with the closest GBIF occurrence for each species, including coordinates, distance (in km), date, country, locality, and dataset key.
- A CSV file (e.g. `species_without_occurrence.csv`) listing species with no occurrences found within the specified radius.
- Cached JSON responses stored in the directory specified by `--cache_dir`.


---

## Getting Started (coming soon)

Once the Dockerfile and workflow are finalized, instructions will be added here for:

- Cloning the repository
- Building the Docker container
- Running the pipeline with Snakemake
- Interpreting the results

---

## Author

This project was developed by **scmdcunha (Sara Cristina Marques da Cunha)**.

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**.
