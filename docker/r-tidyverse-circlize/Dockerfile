# Use tidyverse base
FROM rocker/tidyverse:4.5.1
LABEL authors="asherpreskasteinberg"

# install r-circos
RUN R -e 'install.packages("circlize")'
# Set default to bash for script execution flexibility
CMD ["bash"]
