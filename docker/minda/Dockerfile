# Use the Miniconda base image
FROM quay.io/condaforge/mambaforge:24.9.2-0
LABEL authors="asherpreskasteinberg"

# Clone the Minda repository
RUN git clone https://github.com/shahcompbio/minda.git && \
    cd minda && git checkout 3e4881b && \
    git rev-parse --short HEAD > /opt/minda_version.txt && \
    git rev-parse HEAD > /opt/minda_version_full.txt && \
    mamba env create -f environment.yml -n minda && \
    mamba clean -a -y

# Configure the container to activate the environment and switch to the minda directory on start
# Activate conda environment
RUN echo "conda activate minda" >> ~/.bashrc
ENV PATH /opt/conda/envs/minda/bin:$PATH

# Install additional packages in the "minda" environment
RUN mamba install -n minda -c conda-forge intervaltree -y
