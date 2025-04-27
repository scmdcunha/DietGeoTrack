FROM snakemake/snakemake:v9.1.7

SHELL ["bash", "-l", "-c"]
WORKDIR /data

# Creates micromamba environment with necessary packages
RUN micromamba create -n metabarcoding -y \
    -c conda-forge -c bioconda \
    vsearch=2.30.0 \
    blast=2.16.0 \
    pandas=2.2.3 \
    gdal=3.7.2

# Activates the environment
ENV PATH=/opt/conda/envs/metabarcoding/bin:$PATH
