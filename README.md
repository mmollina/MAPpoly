[![Build Status](https://travis-ci.org/mmollina/MAPpoly.svg?branch=master)](https://travis-ci.org/mmollina/MAPpoly) [![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Introduction

`mappoly` (v. 0.1.0) is an under development R package to construct genetic maps in autopolyploids with even ploidy levels. In its current version, `mappoly` can handle ploidy levels up to 8 when using hidden Markov models (HMM), and up to 12 when using the two-point simplification. All the two-point based functions are fast enough to run on standard computers. However, we strongly recommend to use high-performance computation for HMM-based analysis, especially for ploidy levels higher than 4. 

We assume the genotypic data is available and in the format required by `mappoly`. In a future version, this document will include instructions about genotype calling and `vcf` files. The derivation of the HMM used in `mappoly` can be found in [Mollinari and Garcia 2018](https://doi.org/10.1101/415232 ). 

`mappoly` is not available in CRAN, but you can install it from Git Hub. Within R, you need to install and load the package `devtools`:

```R
install.packages("devtools")
```
To install `mappoly` from Git Hub use

```R
devtools::install_github("mmollina/mappoly")
```

# Tutorial

1. [Building a genetic map in an hexaploid full-sib population using MAPpoly](https://mmollina.github.io/MAPpoly/)

# Acknowledgment 

This package is part of the of the [Genomic Tools for Sweetpotato Improvement project](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP), founded by [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/).
