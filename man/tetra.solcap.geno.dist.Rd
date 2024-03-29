% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tetra.solcap.geno.dist}
\alias{tetra.solcap.geno.dist}
\title{Autotetraploid potato dataset with genotype probabilities.}
\format{
An object of class \code{mappoly.data} which contains a
    list with the following components:
\describe{
    \item{ploidy}{ploidy level = 4}
    \item{n.ind}{number individuals = 160}
    \item{n.mrk}{total number of markers = 4017}
    \item{ind.names}{the names of the individuals}
    \item{mrk.names}{the names of the markers}
    \item{dosage.p1}{a vector containing the dosage in
      parent P for all \code{n.mrk} markers}
    \item{dosage.p2}{a vector containing the dosage in
      parent Q for all \code{n.mrk} markers}
    \item{chrom}{a vector indicating which sequence each marker
      belongs. Zero indicates that the marker was not assigned to any
      sequence}
    \item{genome.pos}{Physical position of the markers into the
      sequence}
    \item{prob.thres = 0.95}{probability threshold to associate a marker 
                             call to a dosage. Markers with maximum genotype 
                             probability smaller than 'prob.thres' are considered 
                             as missing data for the dosage calling purposes}
    \item{geno}{a data.frame 
      containing the probability distribution for each combination of
      marker and offspring. The first two columns represent the marker
      and the offspring, respectively. The remaining elements represent
      the probability associated to each one of the possible
      dosages}
      \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
      for each individual (columns). Missing data are represented by 
      \code{ploidy_level + 1 = 5}}
    \item{n.phen}{There are no phenotypes in this simulation}
    \item{phen}{There are no phenotypes in this simulation}
}
}
\usage{
tetra.solcap.geno.dist
}
\description{
A dataset of the B2721 population which derived from a cross between 
two tetraploid potato varieties: Atlantic × B1829-5. The population comprises 160 
offsprings genotyped with the SolCAP Infinium 8303 potato array. The original data 
set can be found in [The Solanaceae Coordinated Agricultural Project (SolCAP) webpage](http://solcap.msu.edu/potato_infinium.shtml) 
The dataset also contains the genomic order of the SNPs from the Solanum 
tuberosum genome version 4.03. The genotype calling was performed using the
fitPoly R package. Although this dataset contains the
probability distribution of the genotypes, 
it is essentially the same dataset found in \code{\link[mappoly]{tetra.solcap}}
}
\keyword{datasets}
