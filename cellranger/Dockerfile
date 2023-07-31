# image adapted from https://github.com/hisplan/docker-cellranger/

# setup OS
## no reason to use this other than it's recent
FROM ubuntu:22.04
RUN apt-get update

ARG CELLRANGER_VERSION
ARG DOWNLOAD_URL

# get wget
RUN apt-get update && apt-get install --yes build-essential wget

# cell ranger binaries
RUN wget -O cellranger-${CELLRANGER_VERSION}.tar.gz ${DOWNLOAD_URL} \
    && tar -xzf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && rm -rf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && mv cellranger-${CELLRANGER_VERSION} /opt/

## add to path
ENV PATH=$PATH:/opt/cellranger-${CELLRANGER_VERSION}

WORKDIR /
VOLUME /rootvol
WORKDIR /rootvol
#CMD ["/bin/bash"]
ENTRYPOINT ["cellranger"]