# eSCAN

eSCAN: Scan Regulatory Regions for Aggregate Association Testing using Whole Genome Sequencing Data

## Description

eSCAN (or “Scan the Enhancers”, with “enhancers” as a shorthand for any potential regulatory regions in the genome) is an R package for performing genome-wide assessment of rare variant enhancer regions in sequencing studies, combining the advantages of dynamic window selection with the advantages of increasing genomic annotation information, including chromatin accessibility, histone markers, and 3D chromatin conformation.

## Dependencies

eSCAN links to R packages [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html), and also imports R packages [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [GENESIS](https://bioconductor.org/packages/release/bioc/html/GENESIS.html). These dependencies should be installed before installing eSCAN.

## Installation

```
devtools::install_github("yingxi-kaylee/eSCAN")
```

## Version

The current version is 0.1.0 (June 20, 2020).

## License

This software is licensed under GPLv3.
