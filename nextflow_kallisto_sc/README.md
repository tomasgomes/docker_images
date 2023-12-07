## Description

This image is used to run a kallisto-bustools nextflow pipeline. The entrypoint of the image is `nextflow`. The nextflow pipeline can be found [here](https://github.com/tomasgomes/quantification_nextflow).

TODO:

-   add support for velocity pipeline

    -   should run both the normal and velocity pipeline
  
-   fix MT genes

-   add scanpy option

    -   need to do emptyDrops first and save results independently
  
    -   option for creating scanpy (Seurat should always run as default)
  
    -   scvelo setup for RNA velocity
  
    -   account for spatial as well

-   add bustools umicorrect

- streamline code for RNA velocity

## Build image

```{bash}
docker build -t nextflow_kallisto_sc ./
# before building, go to
## https://www.10xgenomics.com/support/software/space-ranger/downloads
## (or your preferred version) to get the download URL
## and set the version accordingly
docker build -t nextflow_kallisto_sc --build-arg SPACERANGER_VERSION="2.1.1" --build-arg DOWNLOAD_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz?Expires=1701729264&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=FZx5Wa0rXXktHv-C~q6mKlr07E2d5THecp8Be-iF0uefEdnUU-0iM4vw-2WhES2Or36Z3GO7gqcC7CGavvFF5cpAXGUkvjobYZGj62tRY8POoSekm-ii7tz4ez4XAijINmm60Km~KSBHBlUb~yMEmRfTxG17uXI4TjfhomATd7mWAbld5aYjxBQmeY-Fk2pzUgBD9Ne4uRinxlSjkm8ICMJKLKw7yslKK0LxJA4Ti9j8oP8s4D3osTxWMiPeYTb~2bAwTaMy8PBDnC2i9hWExYKR3K7pZI1-miXHT3V~56IM27H4J-xgLQ-DUpeVQ99TatzIw4ohRZbvsQYxC7ZIWw__" .
```

## Pulling the image

```{bash}
# a script to pull it on a computing is required when the image is too large, as it may fill the /tmp folder
sbatch singularity_pull_img.sh
```

## Run

```{bash}
# 10xv3
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/10xv3_whitelist.txt --samplename "test_10xv3" --outdir ./ --protocol 10xv3 --reads "test_datasets/10xv3/Brain_Tumor_3p_fastqs/*.fastq.gz" --cores 4

# visium fresh-frozen
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/mouse/gencode.vM10.transcripts_simp.fa.gz --transindex gencode.vM10.transcripts_simp.kalid --t2g data/references/mouse/gencode.vM10.transcripts_t2g.txt --white data/references/technical/barcode_whitelists/visium-v1_whitelist_kallisto.txt --samplename "test_visium_FF" --outdir ./ --protocol visiumv1 --reads "test_datasets/visium_FF/Visium_Mouse_Olfactory_Bulb_fastqs/*.fastq.gz" --images "V10N30-322" --imagear "A1" --imagef test_datasets/visium_FF/Visium_Mouse_Olfactory_Bulb_image.tiff --imageal test_datasets/visium_FF/Visium_Mouse_Olfactory_Bulb_alignment_file.json --cores 2

#visium FFPE
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/barcode_whitelists/visium-v1_whitelist_kallisto.txt --samplename "test_visium_FFPE" --outdir ./ --protocol visiumv1 --reads "test_datasets/visium_FFPE/Visium_FFPE_Human_Ovarian_Cancer_fastqs/*.fastq.gz" --images "V10L13-020" --imagear "D1" --imagef test_datasets/visium_FFPE/Visium_FFPE_Human_Ovarian_Cancer_image.jpeg --cores 2

# sc5pe
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --white data/references/technical/737K-august-2016.txt --samplename "test_sc5pe" --outdir ./ --protocol sc5pe --reads "data/published/Hong2021/HC1/*.fastq.gz" --cores 2

# Bulk, one sample using quant
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --samplename "test_bulk" --outdir ./ --protocol bulk_quant --cores 4 --reads "./test_datasets/SS2/fastq/18689_2390*.fastq.gz"

# ParseBio (uses SPLIT-SEQ protocol)


# SS2 - HAS NOT PASSED TESTING
## REQUIRES ABSOLUTE PATH IN BATCH FILE (/rootvol/...)
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --transindex human_Ens109_GRCh38p13.kalid --t2g data/references/human/human_Ens109_GRCh38p13_t2g.txt --samplename "test_plate" --outdir ./ --protocol batch --cores 4 --reads "./test_datasets/SS2/kallisto_batch.txt"

# RNA velocity
## if the genome ends in .fa.gz, the index must be .fa.gz.fai (but not gzipped!)
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --velomode true --genome data/references/human/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf data/references/human/genome/Homo_sapiens.GRCh38.110.gtf.gz --overhang 100 --cores 4 --white data/references/technical/10xv3_whitelist.txt --protocol 10xv3 --reads "test_datasets/10xv3/Brain_Tumor_3p_fastqs/*.fastq.gz"
```




