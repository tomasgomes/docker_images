FROM ubuntu:22.04
RUN apt-get update

#nextflow
WORKDIR /opt
RUN apt-get update && apt-get install -y procps default-jre default-jdk wget
## only up to this version to work with kallisto_nextflow.nf pipeline
##ENV NXF_VER=22.10.6
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
## add nextflow to path
ENV PATH=$PATH:/opt

WORKDIR /
VOLUME /rootvol
WORKDIR /rootvol
ENTRYPOINT ["nextflow"]