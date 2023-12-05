# docker reference
https://dockerlabs.collabnix.com/docker/cheatsheet/

# slurm reference
https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/slurm-cheatsheat.html#gsc.tab=0

# pull image from Docker Hub
## images will be stored in /mnt/beegfs/singularity/images
singularity pull nextflow_kallisto_sc.sif docker://tomasgomes/nextflow_kallisto_sc:0.2

# running docker container
singularity run /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif

# running docker container with volume mounted (to access a local file)
singularity run --mount type=bind,src=$(pwd),dst=/rootvol /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run /rootvol/kallisto_pipeline.nf

# running the above command with sbatch
sbatch example_sbatch_script.sh

# which is equivalent to the following (except this will leave the terminal unusable)
srun --job-name=test_wrapper --time=0:01:20 --ntasks=1 --mem=100M --cpus-per-task=1 --nodelist=compute-1 -o outfile.txt -e errfile.txt singularity run --mount type=bind,src=$(pwd),dst=/rootvol /mnt/beegfs/singularity/images/nextflow_kallisto_sc.sif run /rootvol/kallisto_pipeline.nf ...
