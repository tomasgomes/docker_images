## Description

An example of how to run souporcell. The package is already organised into a Docker image and Singularity, with a full python pipeline, so this folder only provides examples of how to run it. See more in the [GitHub repo](https://github.com/wheaton5/souporcell/tree/master), including how to download the files required for the example run.

TODO:

-   test with the --known_genotypes and --known_genotypes_sample_names options

## Get image

```{bash}
singularity pull shub://wheaton5/souporcell
```

## Run

```{bash}
# -B binds a local path to a (dummy) path in the container
## useful for accessing local files
## multiple bindings are possible in a comma-separated list
singularity exec -B $PWD:/dummy /mnt/beegfs/apptainer/images/souporcell_latest.sif souporcell_pipeline.py -i /dummy/possorted_genome_bam.bam -b /dummy/barcodes.tsv -f /dummy/reference.fasta -t 16 -o output_dir_name -k 4
```
