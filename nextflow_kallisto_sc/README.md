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
docker build -t nextflow_kallisto_sc --build-arg SPACERANGER_VERSION="2.1.1" --build-arg DOWNLOAD_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz?Expires=1699477452&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=VYvXIKsxFaG7waviCSQlE7ySLybjx7kA-lpG38l~umwcFYmPcmejxnIETDL3Hckgm6TIfQK9dfH5A1LvgBqTIHXkNCiIurVyzDAnkx0NHcJ3c63Gkb0-YWt0KZ8XtDXUnk8tgBbL-sAPYu1bhrJ3ET7D5UFT3VNd3I55MSZ00t0QBAA6iFK5jSpolwCVOBrOhoZOuPBwPB94jPXHoZPIH8HS2bbMCg85hNpYKEvtgmcEh269UUU~cwLZtHDYJMOKhye4i6JawYs4pFM28mOqSD2CPT~Uni~UH9EgssbhejG2qmQJE4rK~Fz3mjxnb1mODiJ~TwyrHMgqeGkqvlbniA__" .
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
