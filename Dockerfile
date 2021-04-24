FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.4
WORKDIR build
# build steps go here
# remember to clear the build directory!
RUN apt-get update
RUN apt-get install -y r-base r-base-dev
# RUN apt-get install -y curl  ## not necessary for this, but needed for rtracklayer
# RUN apt-get install --fix-missing -y libxml2-dev  ## fails due to package issues
RUN Rscript -e "install.packages('optparse', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('XML')"
RUN Rscript -e "BiocManager::install('RCurl')"
RUN Rscript -e "BiocManager::install('GenomicRanges')"

WORKDIR /app
RUN mkdir -p /xchip/tcga/Tools/absolute/releases/v1.5/
COPY v1.5/ /xchip/tcga/Tools/absolute/releases/v1.5/
COPY src/*.py /usr/local/bin/