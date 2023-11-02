## Description

This image is used to run a kallisto-bustools nextflow pipeline. The entrypoint of the image is `nextflow`. The nextflow pipeline can be found [here](https://github.com/tomasgomes/quantification_nextflow).

TODO:

-   test in Lobo

-   add spaceranger to allow for Visium

-   add support for velocity pipeline

-   update to latest nextflow

## Build image

```{bash}
docker build -t nextflow_kallisto_sc ./
# before building, go to
## https://www.10xgenomics.com/support/software/space-ranger/downloads
## (or your preferred version) to get the download URL
## and set the version accordingly
docker build -t nextflow_kallisto_sc --build-arg SPACERANGER_VERSION="2.1.1" --build-arg DOWNLOAD_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz?Expires=1698957977&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=mlXpuDhn5JpA9afdYzp7Z79NstFSrZC3FfaEvQKCLYbJWdHEkx~rJTDZ8fRUaUtHQ32qttvN-idGoo2HIAGOA4aEqfElGBPLPiQ1s4rGe-8JqCIYeHivPvdxQYvqlSkarlpTmOz7N9-0iLvDp~EDVRrFNtDcCtJX4MUDU9pgsXMpWQYnU6SnS-1ovFKaanIVNKNX2zwk8gRDrvOBDqQibjUXDsX2afw54hgY7uF5fY-c5f8zmTcccPEflYQMZmPcSaweniQ00O~j3dy8X5H8eqlfb20z7RCa2z2A5uUanqhrEAPTQp42KIq6mPyQAd-L0fTW1eOyZoIWw2~SqPIdCA__" .
```

## Run

```{bash}
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run projects/quantification_nextflow/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/737K-august-2016.txt --samplename "test" --outdir ./ --protocol sc5pe --reads "data/published/HC1/*.fastq.gz" --cores 2
```
