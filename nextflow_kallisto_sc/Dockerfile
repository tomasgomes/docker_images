FROM ubuntu:22.04
RUN apt-get update

# kallisto - see https://github.com/stevetsa/Docker-kallisto/blob/master/Dockerfile
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 zlib1g-dev \
 cmake \
 libhdf5-dev \
 git-all \
 autoconf \
 automake \
 wget \
 libcurl4-openssl-dev
## Get latest source from releases
WORKDIR /opt
RUN git clone https://github.com/pachterlab/kallisto.git
WORKDIR kallisto
RUN git clone https://github.com/samtools/htslib
RUN rm -rf -f build
RUN rm -rf /ext/htslib
RUN cp -r htslib /ext/
WORKDIR /opt/kallisto/ext/htslib
RUN autoconf && autoheader
WORKDIR /opt/kallisto
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN make
RUN make install


# bustools
WORKDIR /opt
RUN git clone --single-branch --branch devel https://github.com/BUStools/bustools.git
RUN mkdir bustools/build
RUN cd bustools/build && cmake .. && make && make install


# python basics
## more stuff to install bc it may be important
RUN update-ca-certificates
RUN apt-get install -y libssl-dev libpng-dev libboost-all-dev libxml2-dev openjdk-8-jdk python3-dev python3-pip libfftw3-dev libgsl-dev
RUN apt-get install -y llvm-11
RUN apt-get install -y libgeos-dev
## python stuff
RUN LLVM_CONFIG=/usr/lib/llvm-11/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy umap-learn


# R
RUN apt-get update && apt-get install -y --no-install-recommends ed less locales vim-tiny ca-certificates fonts-texgyre && rm -rf /var/lib/apt/lists/* # buildkit
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && /usr/sbin/update-locale LANG=en_US.UTF-8 # buildkit
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
#ENV R_BASE_VERSION=4.3.1 # NOT USED
RUN apt-get update && apt-get install -y --no-install-recommends libopenblas0-pthread littler r-cran-docopt r-cran-littler r-base r-base-dev r-base-core r-recommended && chown root:staff "/usr/local/lib/R/site-library" && chmod g+ws "/usr/local/lib/R/site-library" && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r && ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r && ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r && rm -rf /tmp/downloaded_packages/ /tmp/*.rds && rm -rf /var/lib/apt/lists/* # buildkit


# R packages
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'DropletUtils', 'Biobase', 'DelayedArray'))"
RUN R --no-echo --no-restore --no-save -e "install.packages(c('Seurat', 'Matrix'))"
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"


#nextflow
WORKDIR /opt
RUN apt-get update && apt-get install -y procps default-jre default-jdk
ENV NXF_VER=22.10.6
RUN wget -qO- https://get.nextflow.io | bash
RUN chmod +x nextflow
## add nextflow to path
ENV PATH=$PATH:/opt

