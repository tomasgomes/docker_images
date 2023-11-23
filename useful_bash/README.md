## Description

This folder contains several files that serve as general guides for using Docker, Singularity, Bash, etc.

Many of the files will be tailored to iMM's Lobo server.

File descriptions:

 - `.bash_profile`: a file with a collection of useful things to add to your server's `.bash_profile`. **WARNING**: if you already have a `.bash_profile` file in your home (`~`) folder, do not simply replace it with the one presented here, since you may lose other things set up in your previous file. Instead, append the suggested lines to your original `.bash_profile`.
 
 - `step_by_step_singularity_use.md`: simple annotated commands of how to use singularity on the server, including with sbatch.
 
 - `example_sbatch_script`: an example script including multiple options for running sbatch.
