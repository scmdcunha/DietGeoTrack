#
# DietGeoTrack Dockerfile
#

# Pull base image
FROM snakemake/snakemake:v9.1.7
SHELL ["/bin/bash", "-c"]
WORKDIR /data

# Install all necessary software
RUN MAMBA_NO_BANNER=1 micromamba create -n metabarcoding -y \
    -c bioconda -c conda-forge \
    blast=2.13.0 \
    vsearch=2.22.1 \
    pandas=2.2.2 \
    plotly=5.21.0 \
    folium=0.16.0 \
    requests=2.31.0 \
    tqdm=4.66.2 \
    biopython=1.83 \
    ete3=3.1.3 \
    gdal=3.6.2 \
    python=3.10

ENV PATH="/opt/conda/envs/metabarcoding/bin:$PATH"

CMD ["/bin/bash"]
