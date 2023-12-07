## Description

This folder contains a Docker image to use the nf-core pipeline to download data from public repositories.

## Build image

```{bash}
docker build -t data_download_nextflow ./
```

## Run

```{bash}
docker run --rm -t -v $(pwd):/rootvol data_download_nextflow run run nf-core/fetchngs --input ids.csv --outdir /rootvol/
```



