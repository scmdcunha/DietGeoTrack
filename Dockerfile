FROM snakemake/snakemake:v9.1.7

WORKDIR /data

# Copy environment.yml into the container
COPY environment.yml /tmp/environment.yml

# Create metabarcoding environment with micromamba
RUN micromamba create -y -n metabarcoding -f /tmp/environment.yml -c conda-forge -c bioconda \
    && micromamba clean --all --yes

# Ensures the environment is in the PATH
ENV PATH="/opt/conda/envs/metabarcoding/bin:$PATH"

# Use bash by default
CMD ["/bin/bash"]
