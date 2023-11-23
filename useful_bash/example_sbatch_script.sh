#!/bin/bash
#SBATCH --job-name=test_wrapper
#SBATCH --time=0:02:00
#SBATCH --ntasks=1
#SBATCH --mem=20M
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=compute-1
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

# This is the code
#srun echo "Hello" && sleep 20 && echo "world"
singularity run --mount type=bind,src=$(pwd),dst=/rootvol /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run /rootvol/kallisto_pipeline.nf