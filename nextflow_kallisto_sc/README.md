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
docker build -t nextflow_kallisto_sc --build-arg SPACERANGER_VERSION="2.1.1" --build-arg DOWNLOAD_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz?Expires=1699385788&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=OSQP0ezQFhch7urKgPG2KKGD4m6QSgcw9kL9y~zxHXqYk8mM9NWh-i0bd~I5ok1NpY9Kwkm4ffG2uKywweLkWpNfWrdm3iTivww11wHU9Abs69lQvDAo3anNfRUwDOKewN2~wEGbs2S~moVGVtq49PwQ1WisBwEu7ScURcr3IVAs-xZeSdAJXc8GsJQOUL-e7JlK4BYiGTwIgjV4w2ZMwsMA3E6ln93dULdY3cU7WN5VzjBoz8O~2f2JNIBe4AAc6vhvia5irh375igiyP3SFiDgTTW0nO8FK1I4YcF~oNgm18FvCZ6HE7cSIHxvRFjAPG~Xxy6XGM7D1i22aDa5qg__" .
```

## Run

```{bash}
# 10xv3
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/10xv3_whitelist.txt --samplename "test_10xv3" --outdir ./ --protocol 10xv3 --reads "test_datasets/10xv3/Brain_Tumor_3p_fastqs/*.fastq.gz" --cores 2

# visium fresh-frozen


#visium FFPE


# sc5pe - OLD, TEST AGAIN
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/737K-august-2016.txt --samplename "test" --outdir ./ --protocol sc5pe --reads "data/published/HC1/*.fastq.gz" --cores 2
```
