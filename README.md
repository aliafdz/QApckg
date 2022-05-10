# QApckg
This package provides a set of functions for NGS data processing, 
  quality analysis, filtering and demultiplexing. These functions are designed
  to be applied in consecutive order on Miseq raw data to obtain a set of 
  intersected haplotypes for each evaluated sample. With this consensus haplotypes 
  different kinds of computations can be made, i.e genotyping, variant calling and
  quasispecies diversity.

## How to install this package

`QApckg` is not available in public repositories. To install this package, download
  the `QApckg_0.1.0.tar.gz` file and execute this code: 
  
```{r}
install.packages("./QApckg_0.1.0.tar.gz", repos = NULL, type = "source")
```
