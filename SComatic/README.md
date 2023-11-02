## Description

An image to run the SComatic pipeline.

TODO:

-   test in Lobo

## Build image

```{bash}
# this may occasionally fail for reasons unknown, complaining about package
## hashes. if so, just try again
docker build -t scomatic_pipeline ./
```

## Run

```{bash}
# with volume mounted in /rootvol/localdir; dry run
docker run --rm -i -t -v $(pwd):/rootvol/localdir scomatic_pipeline
# proper run
docker run --rm -i -t -v $(pwd):/rootvol/localdir scomatic_pipeline run localdir/docker_images/SComatic/SComatic_pipeline.nf --samplename testchr1 --outdir localdir/testchr1/ --genome /rootvol/localdir/data/references/human/refdata-gex-GRCh38-2020-A/fasta/genome.fa --bam /rootvol/localdir/bam/possorted_genome_bam_chr1.bam --ct /rootvol/localdir/projects/Hong2021_SS_blood/results/ct_df.tsv
```
