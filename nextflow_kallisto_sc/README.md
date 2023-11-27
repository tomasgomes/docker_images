## Description

This image is used to run a kallisto-bustools nextflow pipeline. The entrypoint of the image is `nextflow`. The nextflow pipeline can be found [here](https://github.com/tomasgomes/quantification_nextflow).

TODO:

-   test in Lobo

-   add support for velocity pipeline

-   add bustools umicorrect

## Build image

```{bash}
docker build -t nextflow_kallisto_sc ./
# before building, go to
## https://www.10xgenomics.com/support/software/space-ranger/downloads
## (or your preferred version) to get the download URL
## and set the version accordingly
docker build -t nextflow_kallisto_sc --build-arg SPACERANGER_VERSION="2.1.1" --build-arg DOWNLOAD_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.1.tar.gz?Expires=1700166299&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=bnKLW16EhQ3IOGxf5jIFp9kuc9lgsF4fBFWX8R0-neFe34FPLEErTMyrBTNFHClBwDkwH96hL1lNTc67Ig4M6Bdkmkp5isAKfIHpOGYZ8RzQU-b4YSSEn4xaBG1r-u7uaQUBjrKjuUGwcueGyRJ22XSaHlAeTA-85RvH5yajJYCp9socvGnQDUyTEEtG~s8AjnwNyeHhd2jM2YsRFAIKG9lYqD0OjhFC7Hp98efP0AsfBdKwmGWY9xlCOLwEyW4um7~1b9jKAOHprqIXcmUPl61fXoh-NplytIg-O5K9goGUitRy~JSDX-ZR6RX~UMD7LvBxJBJZn11JNy-7Q-34xw__" .
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
docker run --rm -t -v $(pwd):/rootvol nextflow_kallisto_sc run docker_images/nextflow_kallisto_sc/kallisto_pipeline.nf --transcriptome data/references/human/human_Ens109_GRCh38p13.fa.gz --velomode true --genome data/references/human/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --gtf data/references/human/genome/Homo_sapiens.GRCh38.110.gtf.gz --overhang 100 --cores 4
```




