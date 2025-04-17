# Diet Geo Track

**Diet Geo Track** is a bioinformatics pipeline designed to identify arthropods consumed by bats through the analysis of guano samples collected in Serra da Estrela, Portugal. The project focuses on processing and analyzing sequence data obtained via Next-Generation Sequencing (NGS), with the goal of validating taxonomic identifications geographically using GBIF occurrence data. This ensures that the identified species make ecological and biogeographical sense.

> ⚠️ **Note:** This project is in its early development phase. The current focus is on testing and comparing results from BLAST and VSEARCH. The Dockerfile and Snakefile are still under construction.

## Project Overview

- **Input**: FASTA sequences (~5000) obtained through NGS of bat guano samples. These sequences were provided externally and likely generated via metabarcoding (the lab work was not performed in this project).
- **Reference database**: Arthropod COI gene sequences downloaded from the NCBI database using:
    ```
    arthropoda[organism] AND COI[gene]
    ```
- **Taxonomic assignment**: Performed with `vsearch` and `BLAST`, selecting the **top 3 hits per query** with **identity thresholds ≥ 95%**. These are then filtered to assign taxonomic ranks:
    
    - ≥ 99% identity → species
        
    - ≥ 97% → genus
        
    - ≥ 95% → family
        
    - ≥ 90% → order

- **Geographic validation**: Top hits are compared against GBIF occurrence records to filter out improbable identifications based on geography.
    
- **Automation goal**: The pipeline is being developed with Snakemake and containerized via Docker and Micromamba, to ensure reproducibility and ease of use for biologists and researchers.
    
## Tools and Technologies

- [Snakemake](https://snakemake.readthedocs.io/)
    
- [vsearch](https://github.com/torognes/vsearch)
    
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    
- GBIF API
    
- [Docker](https://www.docker.com/)
    
- [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
    
- [Python](https://www.python.org/) and pandas
    
## Getting Started (coming soon)

Once the Dockerfile and workflow are finalized, instructions will be added here for:

- Cloning the repository
    
- Building the Docker container
    
- Running the pipeline with Snakemake
    
- Interpreting the results
    
## Author

This project was developed by **scmdcunha (Sara Cristina Marques da Cunha)** as part of a of a curricular internship.

## License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**.