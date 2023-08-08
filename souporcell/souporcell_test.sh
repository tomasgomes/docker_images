#!/bin/bash
#SBATCH --job-name=souporcell_testing
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=30000
#SBATCH --cpus-per-task=16
#SBATCH --output=souporcell_testing.out
#SBATCH --error=souporcell_testing.out

srun singularity exec -B $PWD:/dummy /mnt/beegfs/apptainer/images/souporcell_latest.sif souporcell_pipeline.py \
-i /dummy/A.merged.bam -b /dummy/GSM2560245_barcodes.tsv \
-f /dummy/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
-t 16 -o /dummy/demux_data_test -k 4