FROM ubuntu:22.04
RUN apt-get update

# sra-tools
RUN apt-get update && apt-get install -y wget curl
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.7/setup-apt.sh \
  && bash setup-apt.sh
RUN /bin/bash -c 'source /etc/profile.d/sra-tools.sh'
## make sure it really is in the PATH
ENV PATH=$PATH:/usr/local/ncbi/sra-tools/bin/

#nextflow
WORKDIR /opt
RUN apt-get update && apt-get install -y procps default-jre default-jdk
## only up to this version to work with kallisto_nextflow.nf pipeline
##ENV NXF_VER=22.10.6
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
## add nextflow to path
ENV PATH=$PATH:/opt

# python basics
## more stuff to install bc it may be important
RUN update-ca-certificates
RUN apt-get install -y python3-dev python3-pip python-is-python3
RUN pip install requests

WORKDIR /
VOLUME /rootvol
WORKDIR /rootvol
ENTRYPOINT ["nextflow"]