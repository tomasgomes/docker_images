#!/bin/bash
#SBATCH --job-name=test_wrapper
#SBATCH --time=0:01:30
#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=compute-2
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

## make sure --mem is enough to run the programme, otherwise it will freeze until time runs out

# This is the code to run
#singularity run --mount type=bind,src=$(pwd),dst=/rootvol \
#/mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run /rootvol/kallisto_pipeline.nf
singularity run --mount type=bind,src=$(pwd),dst=/rootvol \
        /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run \
        /rootvol/kallisto_pipeline.nf --transcriptome \
        /rootvol/human_Ens109_GRCh38p13.fa.gz --velomode true \
        --genome /rootvol/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        --gtf /rootvol/genome/Homo_sapiens.GRCh38.110.gtf.gz \
        --overhang 100 --cores 4

