#!/bin/bash
#SBATCH --job-name=test_wrapper
#SBATCH --time=0:01:30
#SBATCH --ntasks=1
#SBATCH --mem=100M
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=compute-1
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

## make sure --mem is enough to run the programme, otherwise it will freeze until time runs out

# This is the code to run
singularity run --mount type=bind,src=$(pwd),dst=/rootvol \
        /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run \
        /rootvol/kallisto_pipeline.nf


