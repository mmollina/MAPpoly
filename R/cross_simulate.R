#' Simulate an autopolyploid full-sib population
#'
#' Simulate an autopolyploid full-sib population with one or two 
#' informative parents under random chromosome segregation.
#'
#' \code{parental.phases.p} and \code{parental.phases.q} are lists of vectors
#'  containing linkage phase configurations. Each vector contains the
#'  numbers of the homologous chromosomes in which the alleles are
#'  located. For instance, a vector containing \eqn{(1,3,4)} means that
#'  the marker has three doses located in the chromosomes 1, 3 and 4. For
#'  zero doses, use 0.
#'  For more sophisticated simulations, we strongly recommend using PedigreeSim V2.0
#'  \url{https://github.com/PBR/pedigreeSim}
#'
#' @param parental.phases a list containing the linkage phase information for both parents
#' 
#' @param map.length the map length
#' 
#' @param n.ind number of individuals in the offspring
#' 
#' @param draw if \code{TRUE}, draws a graphical representation of the
#'     parental map, including the linkage phase configuration, in a
#'     pdf output (default = FALSE)
#'     
#' @param file name of the output file. It is ignored if
#'     \code{draw = TRUE}
#'  
#' @param prefix prefix used in all marker names. 
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
#'     h.temp <- sim_homologous(ploidy = 6, n.mrk = 20)
#'     fake.poly.dat <- cross_simulate(h.temp, map.length = 100, n.ind = 200)
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
cross_simulate <- function(parental.phases, 
                           map.length,
                           n.ind, 
                           draw = FALSE,
                           file = "output.pdf",
                           prefix = NULL,
                           seed = NULL,
                           width = 12,
                           height = 6,
                           prob.P = NULL,
                           prob.Q = NULL)
{
  n.mrk <- length(parental.phases$p)
  map.temp <- seq(0, map.length, length.out  = n.mrk)
  rf.vector <- mf_h(diff(map.temp))
  ploidy <- parental.phases$ploidy
  x <- sim_cross_two_informative_parents(ploidy,
                                         n.mrk,
                                         rf.vector,
                                         n.ind,
                                         parental.phases$hom.allele.p,
                                         parental.phases$hom.allele.q,
                                         seed = seed,
                                         prob.P = NULL,
                                         prob.Q = NULL)
  if(draw == TRUE)
    draw_cross(ploidy, dist.vec = map.temp,
               parental.phases$hom.allele.p,
               parental.phases$hom.allele.q,
               file = file, width = width, height = height)
  geno <- t(x[[1]])
  rownames(geno) <- paste0(prefix, rownames(geno))
  ind.names <- colnames(geno)
  mrk.names <- rownames(geno)
  res<-structure(list(ploidy = ploidy,
                      n.ind = n.ind,
                      n.mrk = n.mrk,
                      ind.names = ind.names,
                      mrk.names = mrk.names,
                      dosage.p1 = parental.phases$p,
                      dosage.p2 = parental.phases$q,
                      chrom = NA,
                      genome.pos = NA,
                      geno.dose = geno,
                      nphen = 0,
                      phen = NULL),
                 class = "mappoly.data")
  Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
  for(i in 0:ploidy)
    for(j in 0:ploidy)
      Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, dP = i, dQ = j)
  Dpop <- cbind(res$dosage.p1, res$dosage.p2)
  M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
  dimnames(M) <- list(res$mrk.names, c(0:ploidy))
  M <- cbind(M, res$geno.dose)
  res$chisq.pval <- apply(M, 1, mrk_chisq_test, ploidy = ploidy)
  return(res)
}
