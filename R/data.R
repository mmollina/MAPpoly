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
#'     \item{sequence}{a vector indicating the sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1 = 7}}
#'     \item{n.phen}{There are no phenotypes in this simulation}
#'     \item{phen}{There are no phenotypes in this simulation}
#'     \item{chisq.pval}{vector containing p-values for all markers associated to 
#'                       the chi-square test for the expected segregation patterns 
#'                       under Mendelian segregation}
#' }
"hexafake"

#' Simulated autohexaploid dataset with genotype probabilities.
#'
#' A dataset of an hipotetical autohexaploid full-sib population 
#' containing three homology groups. This data set contains the
#' probability distribution of the genotypes and 2\% of missing data, 
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
#'       sequence}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{prob.thres = 0.95}{probability threshold to associate a marker 
#'                              call to a dosage. Markers with maximum genotype 
#'                              probability smaller than 'prob.thres' are considered 
#'                              as missing data for the dosage calling purposes}
#'     \item{geno}{a data.frame 
#'       containing the probability distribution for each combination of
#'       marker and offspring. The first two columns represent the marker
#'       and the offspring, respectively. The remaining elements represent
#'       the probability associated to each one of the possible
#'       dosages}
#'       \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1 = 7}}
#'     \item{n.phen}{There are no phenotypes in this simulation}
#'     \item{phen}{There are no phenotypes in this simulation}
#' }
"hexafake.geno.dist"

#' Resulting maps from \code{\link[mappoly]{hexafake}}
#'
#' A list containing three linkage groups estimated using the procedure available in  
#' [MAPpoly's tutorial](https://mmollina.github.io/MAPpoly/#estimating_the_map_for_a_given_order)
#' 
#' @format A list containing three objects of class \code{mappoly.map}, each one 
#' representing one linkage group in the simulated data. 
#' 
"maps.hexafake"


#' Autotetraploid potato dataset.
#'
#' A data set of the B2721 population which derived from a cross between 
#' two tetraploid potato varieties: Atlantic × B1829-5. The population comprises 160 
#' offsprings genotyped with the SolCAP Infinium 8303 potato array. The original data 
#' set can be found in [The Solanaceae Coordinated Agricultural Project (SolCAP) webpage](http://solcap.msu.edu/potato_infinium.shtml) 
#' The data set also contains the genomic order of the SNPs from the Solanum 
#' tuberosum genome version 4.03. The genotype calling was performed using the
#' fitPoly R package.
#'
#' @format An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#' \describe{
#'     \item{m}{ploidy level = 4}
#'     \item{n.ind}{number individuals = 160}
#'     \item{n.mrk}{total number of markers = 4017}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.q}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{sequence}{a vector indicating the sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1 = 5}}
#'     \item{n.phen}{There are no phenotypes in this simulation}
#'     \item{phen}{There are no phenotypes in this simulation}
#'     \item{chisq.pval}{vector containing p-values for all markers associated to 
#'                       the chi-square test for the expected segregation patterns 
#'                       under Mendelian segregation}
#' }
"tetra.solcap"

#' Autotetraploid potato dataset with genotype probabilities.
#'
#' A data set of the B2721 population which derived from a cross between 
#' two tetraploid potato varieties: Atlantic × B1829-5. The population comprises 160 
#' offsprings genotyped with the SolCAP Infinium 8303 potato array. The original data 
#' set can be found in [The Solanaceae Coordinated Agricultural Project (SolCAP) webpage](http://solcap.msu.edu/potato_infinium.shtml) 
#' The data set also contains the genomic order of the SNPs from the Solanum 
#' tuberosum genome version 4.03. The genotype calling was performed using the
#' fitPoly R package. Although this data set contains the
#' probability distribution of the genotypes, 
#' it is essentialy the same data set found in \code{\link[mappoly]{tetra.solcap}}
#'
#' @format An object of class \code{mappoly.data} which contains a
#'     list with the following components:
#' \describe{
#'     \item{m}{ploidy level = 4}
#'     \item{n.ind}{number individuals = 160}
#'     \item{n.mrk}{total number of markers = 4017}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.q}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{sequence}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{prob.thres = 0.95}{probability threshold to associate a marker 
#'                              call to a dosage. Markers with maximum genotype 
#'                              probability smaller than 'prob.thres' are considered 
#'                              as missing data for the dosage calling purposes}
#'     \item{geno}{a data.frame 
#'       containing the probability distribution for each combination of
#'       marker and offspring. The first two columns represent the marker
#'       and the offspring, respectively. The remaining elements represent
#'       the probability associated to each one of the possible
#'       dosages}
#'       \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1 = 5}}
#'     \item{n.phen}{There are no phenotypes in this simulation}
#'     \item{phen}{There are no phenotypes in this simulation}
#' }
"tetra.solcap.geno.dist"

#' @name solcap.maps
#' @aliases solcap.map
#' @aliases solcap.mds.map
#' @aliases solcap.prior.map
#' @aliases solcap.err.map
#' 
#' @title Solcap Maps
#' 
#' Resulting maps from \code{\link[mappoly]{tetra.solcap}} and \code{\link[mappoly]{tetra.solcap.geno.dist}} datasets.
#'
#' @format All maps have the same format: a list containing 12 objects of class \code{mappoly.map}, each one 
#' representing one linkage group. 
#' \describe{
#'     \item{solcap.map}{genetic map using genomic order and dosage call}
#'     \item{solcap.mds.map}{genetic map using \code{\link[mappoly]{mds_mappoly} order and dosage call}}
#'     \item{solcap.prior.map}{genetic map using genomic order and prior probability distribution}
#'     \item{solcap.err.map}{genetic map using genomic order, dosage call and global calling error}     
#' }
#' 
#' @rdname solcap.maps
"solcap.dose.map"
#' @rdname solcap.maps
"solcap.mds.map"
#' @rdname solcap.maps
"solcap.prior.map"
#' @rdname solcap.maps
"solcap.err.map"

