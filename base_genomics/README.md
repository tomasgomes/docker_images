## Description

An image with basic genomic tools, that can be run interactively.

Softawre included:

 - samtools
 - bedtools
 - fastqc
 - sratools
 - picardtools
 - bcftools
 - bedops

TODO:

-   test in Lobo

## Build image

```{bash}
# this may occasionally fail for reasons unknown, complaining about package
## hashes. if so, just try again
docker build -t base_genomic ./
```
