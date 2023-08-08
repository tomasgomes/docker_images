## Description

This image is used to run a kallisto-bustools nextflow pipeline. The entrypoint of the image is `nextflow`.

TODO:

-   test in Lobo

-   add spaceranger to allow for Visium

-   add support for velocity pipeline

## Build image

```{bash}
docker build -t nextflow_kallisto_sc ./
```

## Run

```{bash}
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run projects/quantification_nextflow/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/737K-august-2016.txt --samplename "test" --outdir ./ --protocol sc5pe --reads "data/published/HC1/*.fastq.gz" --cores 2
```
