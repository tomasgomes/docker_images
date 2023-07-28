## Description

An image to run cellranger.

## Build image
```{bash}
docker build \
    --tag ${image_name}:${version} \
    --build-arg DOWNLOAD_URL=${download_url} \
    --build-arg VDJ_REFERENCE_VERSION=${vdj_ref_version} \
    --build-arg CELLRANGER_VERSION=${version} ./
```

## Run
```{bash}

```