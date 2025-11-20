FROM quay.io/condaforge/mambaforge:24.9.2-0

# Clone Wakhan and create conda environment (git is pre-installed in mambaforge)
WORKDIR /opt/wakhan
RUN git clone https://github.com/KolmogorovLab/Wakhan.git && \
    cd Wakhan && git checkout 43c3197 && \
    git rev-parse --short HEAD > /opt/wakhan_version.txt && \
    git rev-parse HEAD > /opt/wakhan_version_full.txt && \
    mamba env create -f environment.yml -n wakhan && \
    mamba clean -a -y

# Create wrapper script that activates conda env and runs wakhan
RUN echo '#!/bin/bash\nsource /opt/conda/etc/profile.d/conda.sh\nconda activate wakhan\npython3 /opt/wakhan/Wakhan/wakhan.py "$@"' > /usr/bin/wakhan && \
    chmod +x /usr/bin/wakhan

# Set default command to use conda environment
SHELL ["conda", "run", "-n", "wakhan", "/bin/bash", "-c"]

# Add metadata labels
LABEL org.opencontainers.image.source="https://github.com/KolmogorovLab/Wakhan" \
      org.opencontainers.image.description="Wakhan: Copy number aberration analysis tool"
