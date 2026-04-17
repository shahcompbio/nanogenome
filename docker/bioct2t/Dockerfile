# Use tidyverse base
FROM rocker/tidyverse:4.5.1
LABEL authors="asherpreskasteinberg"

# install BiocT2T and dependencies for T2T gene annotation
RUN R -e 'BiocManager::install("GenomicFeatures")'
RUN R -e 'devtools::install_github("vjcitn/BiocT2T@4a244e4")'
RUN R -e 'BiocT2T::install_early_t2t_txdb()'

# Set default to bash for script execution flexibility
CMD ["bash"]
