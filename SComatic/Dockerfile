FROM ubuntu:22.04
RUN apt-get update

# avoid interactive prompts (I think)
ENV DEBIAN_FRONTEND=noninteractive

# basic stuff
## setting timezone required to avoid interactive prompt later
ENV TZ=Europe/Lisbon
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update \
 && apt-get install --yes build-essential apt-utils libcurl4-openssl-dev \
 git-all autoconf automake wget datamash gzip perl zlib1g-dev libbz2-dev \
 liblzma-dev libssl-dev libncurses5-dev make libncursesw5-dev \
 gcc-multilib cmake libhdf5-dev bzip2 keyutils sqlite tk zstd
 
# samtools and bedtools
RUN apt-get update && apt-get install -y samtools tabix  bedtools

# python basics
RUN update-ca-certificates
RUN apt-get update && \
 apt-get install -y libssl-dev libpng-dev libboost-all-dev libxml2-dev \
 openjdk-8-jdk python3-dev python3-pip libfftw3-dev libgsl-dev

# R
RUN apt-get update \
 && apt-get install -y --no-install-recommends ed less locales vim-tiny \
 ca-certificates fonts-texgyre \
 && rm -rf /var/lib/apt/lists/* # buildkit
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 \
 && /usr/sbin/update-locale LANG=en_US.UTF-8 # buildkit
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
#ENV R_BASE_VERSION=4.3.1 # NOT USED
RUN apt-get update && apt-get install -y --no-install-recommends \
 libopenblas0-pthread littler r-cran-docopt r-cran-littler r-base r-base-dev \
 r-base-core r-recommended \
 && chown root:staff "/usr/local/lib/R/site-library" \
 && chmod g+ws "/usr/local/lib/R/site-library" \
 && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
 && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
 && ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
 && ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r \
 && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
 && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
 && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
 && rm -rf /var/lib/apt/lists/* # buildkit
 
## SComatic
WORKDIR /opt
RUN git clone https://github.com/cortes-ciriano-lab/SComatic.git
WORKDIR SComatic

## Install requirements from package lists
### R
RUN Rscript r_requirements_install.R
### python
RUN pip3 cache purge
RUN pip3 install --upgrade pip && \
 pip3 install --no-cache-dir about-time numpy numpy-groupies pandas \
 scipy pybedtools pysam rpy2

## unpack other files
RUN gunzip PoNs/PoN.scRNAseq.hg38.tsv.gz
RUN gunzip PoNs/PoN.scATACseq.hg38.tsv.gz 
RUN gunzip RNAediting/AllEditingSites.hg38.txt.gz

#nextflow
WORKDIR /opt
RUN apt-get update && apt-get install -y procps default-jre default-jdk
## only up to this version to work with kallisto_nextflow.nf pipeline
ENV NXF_VER=22.10.6
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
## add nextflow to path
ENV PATH=$PATH:/opt

WORKDIR /
VOLUME /rootvol
WORKDIR /rootvol
RUN mv /opt/SComatic ./
#ENTRYPOINT ["nextflow"]




