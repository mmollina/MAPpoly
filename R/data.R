#' Simulated autohexaploid dataset.
#'
#' A dataset of an hipotetical autohexaploid full-sib population 
#' containing three homology groups of 
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
#'     \item{geno}{a data.frame containing probability distribution
#'                 for each combination of marker and offspring}
#'     \item{n.phen}{There is no phenotypes in this simulation}
#'     \item{phen}{There is no phenotypes in this simulation}
#' }
"hexafake"
