# Use tidyverse base
FROM rocker/tidyverse:4.5.1
LABEL authors="asherpreskasteinberg"

# install biocmanager
RUN R -e 'install.packages("BiocManager")'
# install karyoploteR
RUN R -e 'BiocManager::install("karyoploteR")'
RUN R -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'
RUN R -e 'BiocManager::install("BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0")'

# Set default to bash for script execution flexibility
CMD ["bash"]
