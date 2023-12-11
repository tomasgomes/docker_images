#!/bin/bash
#SBATCH --job-name=download_data
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=6
#SBATCH --nodelist=compute-21
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

# before running this, the image has to be pulled to the local nextflow cache (~/.nextflow)
# singularity run --mount type=bind,src=$(pwd),dst=/rootvol /mnt/beegfs/singularity/images/data_download_nextflow.sif pull nf-core/fetchngs

# This is the code to run
singularity run --mount type=bind,src=$(pwd),dst=/rootvol /mnt/beegfs/singularity/images/data_download_nextflow.sif run nf-core/fetchngs --max_memory 31GB --max_cpus 6 --input /rootvol/ids.csv --outdir /rootvol/