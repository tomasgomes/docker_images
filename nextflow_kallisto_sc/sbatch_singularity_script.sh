#!/bin/bash
#SBATCH --job-name=test_wrapper
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=12
#SBATCH --nodelist=compute-2
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

## make sure --mem is enough to run the programme, otherwise it will freeze until time runs out
## velocity takes longer to run, likely because of larger indices and extra steps

# This is the code to run
#singularity run --mount type=bind,src=$(pwd),dst=/rootvol \
#/mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run /rootvol/kallisto_pipeline.nf
singularity run --mount type=bind,src=$(pwd),dst=/rootvol \
        /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run \
        /rootvol/kallisto_pipeline.nf --transcriptome \
        /rootvol/human_Ens109_GRCh38p13.fa.gz --velomode true \
        --genome /rootvol/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        --gtf /rootvol/genome/Homo_sapiens.GRCh38.110.gtf.gz \
        --overhang 100 --cores 8 --white /rootvol/10xv3_whitelist.txt \
        --protocol 10xv3 --samplename "test_velo_10xv3" \
        --t2g /rootvol/human_Ens109_GRCh38p13_t2g.txt \
        --reads "/rootvol/Brain_Tumor_3p_fastqs/*.fastq.gz"

