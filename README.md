# QApckg
## Presentation
This repository has been created to include the original code developed by Alicia Aranda Fernandez in the completion of Bioinformatics and Biostatistics master's thesis, resulting in `QApckg` R package. 

As with any other package, the [R](./R) folder includes the original code for each of the defined functions, while the [man](./man) folder (generated automatically) incorporates all their documentation. The DESCRIPTION (package metadata) and NAMESPACE (list of created functions and dependencies) required files are also presented. 

In addition, the repository includes the [QApckg_0.1.0.pdf](./QApckg_0.1.0.pdf) document corresponding to the package reference manual, and the compressed [QApckg_0.1.0.tar.gz](./QApckg_0.1.0.tar.gz) file which allows the local installation of the package.

It is important to note that all the presented code is derived from a subset of the scripts included in the following repository: <https://github.com/aliafdz/QA_genotipat_pipeline>

## Package description
This package provides a set of functions for NGS data processing, quality analysis, filtering and demultiplexing. These functions are designed to be applied in consecutive order on Miseq raw data to obtain a set of intersected haplotypes for each evaluated sample. With this consensus haplotypes different kinds of computations can be made, i.e genotyping, variant calling and quasispecies diversity.

## How to install this package

`QApckg` is not available in public repositories. To install this package, download
  the `QApckg_0.1.0.tar.gz` file and execute this code: 
  
```{r}
install.packages("./QApckg_0.1.0.tar.gz", repos = NULL, type = "source")
```
