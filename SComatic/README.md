## Description

An image to run the SComatic pipeline.

TODO:

-   test in Lobo

## Build image

```{bash}
# this may occasionally fail for reasons unknown, complaining about package
## hashes. if so, just try again
docker build -t scomatic_pipeline ./
```

## Run

```{bash}
# no args
docker run --rm -i -t scomatic_pipeline
```
