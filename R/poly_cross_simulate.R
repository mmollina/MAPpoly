#' Simulate an autopolyploid full-sib population
#'
#' Simulate an autopolyploid full-sib population with one or two 
#' informative parents under random chromosome segregation.
#'
#' \code{hom.allele.p} and \code{hom.allele.q} are lists of vectors
#'  containing linkage phase configurations. Each vector contains the
#'  numbers of the homologous chromosomes in which the alleles are
#'  located. For instance, a vector containing \eqn{(1,3,4)} means that
#'  the marker has three doses located in the chromosomes 1, 3 and 4. For
#'  zero doses, use 0.
#'  For more sophisticated simulations, we strongly recommend using PedigreeSim V2.0
#'  \url{https://www.wur.nl/en/show/Software-PedigreeSim.htm}
#'
#' @param ploidy ploidy level. Must be an even number
#' 
#' @param rf.vec vector containing the recombination fractions between
#'     adjacent markers. If a single recombination fraction is
#'     provided, it is repeated \eqn{n.mrk-1} times
#'     
#' @param n.mrk number of markers
#' 
#' @param n.ind number of individuals in the offspring
#' 
#' @param hom.allele a list containing the linkage phase information for both parents
#' 
#' @param draw if \code{TRUE}, draws a graphical representation of the
#'     parental map, including the linkage phase configuration, in a
#'     pdf output (default = FALSE)
#'     
#' @param file name of the output file. It is ignored if
#'     \code{draw = TRUE}
#'     
#' @param seed random number generator seed (default = NULL)
#'
#' @param width the width of the graphics region in inches (default = 12)
#'
#' @param height the height of the graphics region in inches (default = 6)
#'
#' @param prob.P a vector indicating the proportion of preferential
#'     pairing in parent P (currently ignored)
#'
#' @param prob.Q a vector indicating the proportion of preferential
#'     pairing in parent Q (currently ignored)
#'
#' @return an object of class \code{mappoly.data}. See
#'     \code{\link[mappoly]{read_geno}} for more information
#'
#' @examples
#'     h.temp <- sim_homologous(ploidy = 6, n.mrk = 20, max.d = 3, max.ph = 3, seed = 123)
#'     fake.poly.dat <- poly_cross_simulate(ploidy = 6, rf.vec = .05, n.mrk = 20,
#'                                   n.ind = 200, h.temp, seed = 123)
#'     plot(fake.poly.dat)
#'                                    
#'                                   
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'
#' @export
poly_cross_simulate <- function(ploidy, rf.vec, n.mrk,
                              n.ind, hom.allele,
                              draw = FALSE,
                              file = "output.pdf",
                              seed = NULL,
                              width = 12,
                              height = 6,
                              prob.P = NULL,
                              prob.Q = NULL)
{
  if(length(rf.vec) == 1) rf.vec <- rep(rf.vec, n.mrk-1)
  x <- sim_cross_two_informative_parents(ploidy,
                                       n.mrk,
                                       rf.vec,
                                       n.ind,
                                       hom.allele$hom.allele.p,
                                       hom.allele$hom.allele.q,
                                       seed = seed,
                                       prob.P = NULL,
                                       prob.Q = NULL)
  if(draw == TRUE)
    draw_cross(ploidy,round(rf.vec,4) ,hom.allele$hom.allele.p,hom.allele$hom.allele.q,
               file = file, width = width, height = height)
  geno <- t(x[[1]])
  ind.names <- colnames(geno)
  mrk.names <- rownames(geno)
  structure(list(ploidy = ploidy,
                 n.ind = n.ind,
                 n.mrk = n.mrk,
                 ind.names = ind.names,
                 mrk.names = mrk.names,
                 dosage.p1 = hom.allele$p,
                 dosage.p2 = hom.allele$q,
                 chrom = NA,
                 genome.pos = NA,
                 geno.dose = geno,
                 nphen = 0,
                 phen = NULL),
            class = "mappoly.data")
}
