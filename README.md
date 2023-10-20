# Docker images

A central repository for all main docker images developed.

Images can be found in the following DockerHub repositories:

-   <https://hub.docker.com/u/tomasgomes>

## Useful docker commands

### Build as an image

```{bash}
docker build -t nextflow_kallisto_sc ./
```

### Open interactive image (console)

```{bash}
# this will be interactive, and remove the container after exit
# will also map the pwd into rootvol
docker run --rm -i -t -v $(pwd):/rootvol nextflow_kallisto_sc
```

### Push an image to dockerhub

```{bash}
# login to Docker ($DOCKERT should be your token/password)
docker login -u tomasgomes -p $DOCKERT
# tag the image you want to upload (i.e. give it a version)
docker tag nextflow_kallisto_sc tomasgomes/nextflow_kallisto_sc:0.1
# check if image was tagged (on the list)
docker images
# push it to dockerhub
docker push tomasgomes/nextflow_kallisto_sc:0.1
```

## TODO

-   base interactive computational image
-   SComatic image
-   SRA download image
-   base jupyter/python image (based off of existing RStudio image)
-   single-cell/spatial jupyter/python image (based off of existing RStudio image)
-   CellPhoneDB image
-   
