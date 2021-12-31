[![R-CMD-check](https://github.com/mmollina/MAPpoly/workflows/R-CMD-check/badge.svg)](https://github.com/mmollina/mappoly/actions)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mmollina/mappoly?branch=master&svg=true)](https://ci.appveyor.com/project/mmollina/mappoly)
![Development](https://img.shields.io/badge/development-active-blue.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/github/mmollina/MAPpoly/branch/master/graphs/badge.svg)](https://codecov.io/github/mmollina/MAPpoly)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mappoly)](https://cran.r-project.org/package=mappoly)
[![R-universe PolyVerse Status Badge](https://polyploids.r-universe.dev/badges/mappoly)](https://polyploids.r-universe.dev/badges/mappoly)
[![CRAN_monthly_downloads](https://cranlogs.r-pkg.org/badges/mappoly)](https://cranlogs.r-pkg.org/badges/mappoly)

<!-- ![](https://raw.githubusercontent.com/mmollina/MAPpoly/master/mappoly_hexsticker.png) -->

# MAPpoly <img src="https://raw.githubusercontent.com/mmollina/MAPpoly/main/hex.png" align="right" width="150" />

MAPpoly (v. 0.3.0) is an R package to construct genetic maps in autopolyploids with even ploidy levels. In its current version, MAPpoly can handle ploidy levels up to 8 when using hidden Markov models (HMM), and up to 12 when using the two-point simplification. When dealing with large numbers of markers (> 10,000), we strongly recommend using high-performance computation. 

![](https://raw.githubusercontent.com/mmollina/MAPpoly/master/mappoly.gif)


In its current version, MAPpoly can handle the following types of datasets:

1. CSV files 
2. MAPpoly files
    - Dosage based
    - Probability based
3. [fitPoly](https://CRAN.R-project.org/package=fitPoly) files
4. VCF files

MAPpoly also is capable of importing objects generated by the following R packages 

1. [updog](https://CRAN.R-project.org/package=updog)
2. [polyRAD](https://CRAN.R-project.org/package=polyRAD)
3. [polymapR](https://CRAN.R-project.org/package=polymapR)
    - Datasets
    - Maps

The mapping strategy is based on using pairwise recombination fraction estimation as the first source of information to position allelic variants in specific homologues sequentially. For situations where pairwise analysis has limited power, the algorithm relies on the multilocus likelihood obtained through a hidden Markov model (HMM). The derivation of the HMM used in MAPpoly can be found in [Mollinari and Garcia, 2019](https://doi.org/10.1534/g3.119.400378). 

# Installation

## From CRAN (stable version)

To install MAPpoly from the The Comprehensive R Archive Network (CRAN) use

```R
install.packages("mappoly")
```

## From GitHub (development version)

You can install the development version from Git Hub. Within R, you need to install `devtools`:

```R
install.packages("devtools")
```

If you are using Windows, you must install the the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install MAPpoly from Git Hub use

```R
devtools::install_github("mmollina/mappoly", dependencies=TRUE)
```

For further QTL analysis, we recommend our [QTLpoly](https://github.com/guilherme-pereira/QTLpoly) package. QTLpoly is an under development software to map quantitative trait loci (QTL) in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model [Pereira et al. 2020](https://doi.org/10.1534/genetics.120.303080). 


# Workflow
![](https://raw.githubusercontent.com/mmollina/MAPpoly/main/MAPpoly_workflow.png)

# Vignettes
* [Building a genetic map in a tetraploid potato full-sib population using MAPpoly](https://rpubs.com/mmollin/tetra_mappoly_vignette)
* [Building a genetic map in an hexaploid full-sib population using MAPpoly](https://mmollina.github.io/tutorials/hexa_fake/haxaploid_map_construction.html)
* Real datasets
  * [Hexaploid sweetpotato VCF dataset (Beauregard x Tanzania) obtained using VCF2SM](https://github.com/mmollina/MAPpoly_vignettes/tree/master/data/BT)
  * [Tetraploid potato with dosage call in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/SolCAP_dosage)
  * [Tetraploid potato with dosage call in CSV format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/tetra_solcap.csv)
  * [Tetraploid potato with dosage probabilities in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/SolCAP)
  * [Tetraploid potato in CSV format obtained using ClusterCall](https://raw.githubusercontent.com/mmollina/B2721_map/master/cluster_call/B2721_CC.csv)
  * [Compressed tetraploid potato with dosage probabilities obtained using fitPoly](https://github.com/mmollina/SCRI/raw/main/data/fitpoly_tetra_call/B2721_scores.zip)
* Simulated datasets
   * [Hexaploid simulation with dosage call in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/hexafake)
   * [Hexaploid simulation with dosage probabilities in MAPpoly format](https://github.com/mmollina/MAPpoly_vignettes/blob/master/data/hexafake_geno_dist)
   
  
# Related software

* [Polyverse](https://polyploids.r-universe.dev/ui#builds) - the polyploid R universe (a Lindsay Clark's initiative)
```R
# Enable this universe
options(repos = c(
    polyploids = 'https://polyploids.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

# Install some packages
install.packages('MAPpoly')
```

* Variant Calling
  *  [GBSapp: An automated pipeline for variant calling and filtering.](https://github.com/bodeolukolu/GBSapp)

* Simulations
  * [PedigreeSim: Simulation of genetic marker data in diploid and polyploid pedigreed populations.](https://www.wur.nl/en/show/Software-PedigreeSim.htm)

* Genotype calling
  * [ClusterCall: Automated tetraploid genotype calling by hierarchical clustering](https://potatobreeding.cals.wisc.edu/software/)
  * [fitPoly: Genotype Calling for Bi-Allelic Marker Assays](https://CRAN.R-project.org/package=fitPoly)
  * [polyRAD: Genotype Calling with Uncertainty from Sequencing Data in Polyploids and Diploids](https://CRAN.R-project.org/package=polyRAD)
  * [SuperMASSA: Graphical Bayesian inference tool for genotyping polyploids](https://bitbucket.org/orserang/supermassa)
  * [updog: Flexible Genotyping for Polyploids](https://CRAN.R-project.org/package=updog)
  * [VCF2SM: Python script that integrates VCF files and SuperMASSA](https://github.com/guilherme-pereira/vcf2sm)
 
* Genetic mapping in polyploids
  * [MDSMap: High Density Genetic Linkage Mapping using Multidimensional Scaling](https://CRAN.R-project.org/package=MDSMap)
  * [polymapR: Linkage Analysis in Outcrossing Polyploids](https://CRAN.R-project.org/package=polymapR)
  * [TetraploidSNPMap: Linkage maps and mapping QTLs for autotetraploid species, using SNP dosage data.](https://www.bioss.ac.uk/knowledge/tetraploidmap/)
  
  
* Haplotype reconstruction
  * [TetraOrigin:haplotype reconstruction in a full-sib tetraploid family](https://github.com/chaozhi/TetraOrigin)
  * [PolyOriginR:haplotype reconstruction in polyploid multiparental populations](https://github.com/chaozhi/PolyOriginR)

* QTL mapping
  * [QTLpoly: QTL mapping in full-sib families of outcrossing autopolyploid species based on a random-effect multiple QTL model](https://cran.r-project.org/package=qtlpoly)
  * [diaQTL: QTL analysis of diploid and autotetraploid diallel populations](https://github.com/jendelman/diaQTL)
  * [polyqtlR: QTL analysis and exploration of meiotic patterns in autopolyploid bi-parental F1 populations.](https://cran.r-project.org/web/packages/polyqtlR/index.html)

* Visualization
  * [VIEWpoly: integrate, visualize and explore results from genetic analysis, together with genomic information for autopolyploids](https://cran.r-project.org/package=viewpoly)

# Miscellaneous
* [Supplementary scripts for Mollinari and Garcia (2019)](https://github.com/mmollina/Autopolyploid_Linkage)
* [Miscellaneous scripts](https://github.com/mmollina/MAPpoly_vignettes/blob/master/README.md)

# Articles referencing MAPpoly

1. Using probabilistic genotypes in linkage analysis of polyploids. ([Liao et al., 2021](https://doi.org/10.1007/s00122-021-03834-x))
2. Discovery of a major QTL for root-knot nematode *Meloidogyne incognita* resistance in cultivated sweetpotato *Ipomoea batatas*. ([Oloka, et al., 2021](https://doi.org/10.1007/s00122-021-03797-z))
3. Quantitative trait locus mapping for common scab resistance in a tetraploid potato full-sib population. ([Pereira et al., 2021](https://doi.org/10.1094/PDIS-10-20-2270-RE))
4. The recombination landscape and multiple QTL mapping in a Solanum tuberosum cv.'Atlantic'-derived F1 population. ([Pereira et al., 2021](https://doi.org/10.1101/2020.08.24.265397))
5. High-Resolution Linkage Map and QTL Analyses of Fruit Firmness in Autotetraploid Blueberry ([Cappai et al., 2020](https://doi.org/10.3389/fpls.2020.562171))
6. When a phenotype is not the genotype: Implications of phenotype misclassification and pedigree errors in genomics-assisted breeding of sweetpotato *Ipomoea batatas* (L.) Lam.([Gemenet et al., 2020](https://doi.org/10.1101/747469 ))
7. Quantitative trait loci and differential gene expression analyses reveal the genetic basis for negatively associated beta-carotene and starch content in hexaploid sweetpotato [*Ipomoea batatas* (L.) Lam.] ([Gemenet et al., 2020](https://doi.org/10.1007/s00122-019-03437-7))
8. Multiple QTL Mapping in Autopolyploids: A Random-Effect Model Approach with Application in a Hexaploid Sweetpotato Full-Sib Population. ([Pereira et al., 2020](https://doi.org/10.1534/genetics.120.303080))
9. Unraveling the Hexaploid Sweetpotato Inheritance Using Ultra-Dense Multilocus Mapping. ([Mollinari et al., 2020](https://doi.org/10.1534/g3.119.400620)).

# Acknowledgment

This package has been developed as part of the [Genomic Tools for Sweetpotato Improvement project](https://sweetpotatogenomics.cals.ncsu.edu/) (GT4SP) and [SweetGAINS](https://cgspace.cgiar.org/handle/10568/106838), both funded by [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/). Its continuous improvement is made possible by [Tools for polyploids](https://www.polyploids.org/), funded by USDA NIFA Specialty Crop Research Initiative Award 
