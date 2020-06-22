# eSCAN

eSCAN: Scan Regulatory Regions for Aggregate Association Testing using Whole Genome Sequencing Data

## Description

eSCAN (or “Scan the Enhancers”, with “enhancers” as a shorthand for any potential regulatory regions in the genome) is an R package for performing genome-wide assessment of rare variant enhancer regions that are significantly associated with a phenotype, combining the advantages of dynamic window selection with the advantages of increasing genomic annotation information, including chromatin accessibility, histone markers, and 3D chromatin conformation. eSCAN defines test windows across the genome by genomic annotation either specified by users or provided by eSCAN package such that each window marks putative regulatory region(s). eSCAN then searches the defined windows using fastSKAT and detect significant regions by an empirical/analytical threshold which adjusts for multiple testing of all the searching windows across the genome, of different sizes, including some overlapping windows.


## Installation

```
install.packages("devtools")
devtools::install_github("yingxi-kaylee/eSCAN")
```

## Usage

Please see the [eSCAN user manual](doc/eSCAN-manual.pdf) for detailed usage of eSCAN package. See the following section for a working example.

## eSCAN Working Example

Under construction. Coming soon!

## Version

The current version is 0.1.1 (June 21, 2020).

## License

This software is licensed under GPLv3.

## Citation

eSCAN