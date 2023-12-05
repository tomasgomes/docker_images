#!/bin/bash
#SBATCH --job-name=getnex
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=compute-2
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt


singularity pull nextflow_kallisto_sc.sif docker://tomasgomes/nextflow_kallisto_sc:0.3