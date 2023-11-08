## Description

An image to run cellranger.

TODO:

-   test in Lobo

## Build image

```{bash}
# before building, go to
## https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest/
## (or your preferred version) to get the download URL
## and set the version accordingly
docker build -t cellranger --build-arg CELLRANGER_VERSION="7.1.0" --build-arg DOWNLOAD_URL="https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.1.0.tar.gz?Expires=1690832643&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2OTA4MzI2NDN9fX1dfQ__&Signature=CZz4kX1Vwr3FgxHEMUDfScbzwTdZysMkyENMAqBIot8P5lrVrcvfbxTHBEWmgpVAp0kti0Z6JSIZyurj3ZcquJyZT~2odH5PodCdLnV3USBrCEu50gS83e-EyEAwtQQlu6L-rLJYb59b2F5Je76F9nQlUT7~MPOR3eWqUdl-ha6Pz183uIW6Ubm1h5g2qjzh4GXPKIURCtwrTx3qS356pGd1Muqza5BlPFOmD5kcSpEaVvibohPPgay0NIhTc~9MJv6PKMctQWu~WR9qbW1mr25~J83X9LwFYfzXOtd8IojnU6Kw-lhD4V34iIoRc1PXtGUNIPg-QNpwyV9xMsLGcQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" .
```

## Run

```{bash}
# no args
docker run --rm -i -t -v $(pwd):/rootvol cellranger
# with args
docker run --rm -t -v $(pwd):/rootvol cellranger count --id=HC1_CR \
--fastqs=data/published/HC1 --sample=Hong2021_HC1 --localcores=2 \
--transcriptome=data/references/human/refdata-gex-GRCh38-2020-A
# make a reference
docker run --rm -t -v $(pwd):/rootvol cellranger mkref --genome=eGFP --fasta=eGFP.fa --genes=eGFP.gtf
```
