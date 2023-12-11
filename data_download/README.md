## Description

This folder contains a Docker image to use the nf-core pipeline to download data from public repositories. This will always have to pull nf-core/fetchngs at the start (though this is quick).

## Build image

```{bash}
docker build -t data_download_nextflow ./
```

## Run

```{bash}
docker run --rm -t -v $(pwd):/rootvol data_download_nextflow run nf-core/fetchngs --max_memory 7GB --max_cpus 2 --input ids.csv --outdir /rootvol/
```

## Singularity pull

```{bash}
singularity pull data_download_nextflow.sif docker://tomasgomes/data_download_nextflow:0.1
```

## Running on the cluster

```{bash}
# before running this, the image has to be pulled to the local nextflow cache (~/.nextflow)
# singularity run --mount type=bind,src=$(pwd),dst=/rootvol /mnt/beegfs/singularity/images/data_download_nextflow.sif pull nf-core/fetchngs
sbatch sbatch_download.sh 
```