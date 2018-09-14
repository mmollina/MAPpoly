#' Simulated autohexaploid dataset.
#'
#' A dataset of an hipotetical autohexaploid full-sib population 
#' containing three homology groups
#'
#' @format An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#' \describe{
#'     \item{m}{ploidy level = 6}
#'     \item{n.ind}{number individuals = 300}
#'     \item{n.mrk}{total number of markers = 1500}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.q}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{sequence}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence.}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}.}
#'     \item{n.phen}{There is no phenotypes in this simulation}
#'     \item{phen}{There is no phenotypes in this simulation}
#' }
"hexafake"

#' Simulated autohexaploid dataset with genotype probabilities.
#'
#' A dataset of an hipotetical autohexaploid full-sib population 
#' containing three homology groups. This data set contains the
#' probability distribution of teh genotypes and 2% of missing data, 
#' but is essentialy the same data set found in \code{\link[mappoly]{hexafake}}
#'
#' @format An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#' \describe{
#'     \item{m}{ploidy level = 6}
#'     \item{n.ind}{number individuals = 300}
#'     \item{n.mrk}{total number of markers = 1500}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.q}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{sequence}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence.}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}.}
#'     \item{geno}{a data.frame 
#'       containing the probability distribution for each combination of
#'       marker and offspring. The first two columns represent the marker
#'       and the offspring, respectively. The remaining elements represent
#'       the probability associated to each one of the possible
#'       dosages.}
#'     \item{n.phen}{There is no phenotypes in this simulation}
#'     \item{phen}{There is no phenotypes in this simulation}
#' }
"hexafake.geno.dist"
