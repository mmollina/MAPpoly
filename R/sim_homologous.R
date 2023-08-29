#' Simulate homology groups
#'
#' Simulate two homology groups (one for each parent) and their
#' linkage phase configuration.
#'
#' This function prevents the simulation of linkage phase
#' configurations which are impossible to estimate via two point
#' methods
#'
#' @param ploidy ploidy level. Must be an even number
#' 
#' @param n.mrk number of markers
#' 
#' @param prob.dose a vector indicating the proportion of markers for
#'    different dosage to be simulated (default = NULL)
#'     
#' @param seed random number generator seed
#'
#' @return a list containing the following components:
#'   \item{hom.allele.p}{ a list of vectors
#'    containing linkage phase configurations. Each vector contains the
#'    numbers of the homologous chromosomes in which the alleles are
#'    located. For instance, a vector containing \eqn{(1,3,4)} means that
#'    the marker has three doses located in the chromosomes 1, 3 and 4. For
#'    zero doses, use 0}
#'   \item{p}{contains the indices of the starting positions of the
#'     dosages, considering that the vectors contained in \code{p} are
#'     concatenated. Markers with no doses (zero doses are also
#'     considered)}
#'   \item{hom.allele.q}{Analogously to \code{hom.allele.p}}
#'   \item{q}{Analogously to \code{p}}
#'    \item{ploidy}{ploidy level}
#'
#' @examples
#'     h.temp <- sim_homologous(ploidy = 6, n.mrk = 20)
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
sim_homologous <- function(ploidy, n.mrk, prob.dose = NULL, seed = NULL)
{
  #prob.dose <- NULL
  if(!is.null(seed)) 
    set.seed(seed)
  
  if(is.null(prob.dose))
    prob.dose <- rep(1/(ploidy + 1), (ploidy + 1))
  
  if (length(prob.dose) != ploidy + 1)
    stop("Length of 'prob.dose' must be ploidy + 1.")
  
  
  
  # Create a contingency table
  contingency_table <- outer(prob.dose, prob.dose)
  
  
  # Remove the unacceptable conditions
  contingency_table[1, ploidy+1] <- 0
  contingency_table[ploidy+1, 1] <- 0
  contingency_table[1, 1] <- 0
  contingency_table[ploidy+1, ploidy+1] <- 0
  
  # Normalize the contingency table
  contingency_table <- contingency_table / sum(contingency_table)
  
  # Initialize matrices P1 and P2 with zeros
  P1 <- matrix(0, n.mrk, ploidy)
  P2 <- matrix(0, n.mrk, ploidy)
  
  for (i in 1:n.mrk) {
    # Generate a pair of sums for P1 and P2 based on the normalized contingency table
    
    pair_idx <- sample(1:(ploidy + 1)^2, 1, prob = as.vector(contingency_table))
    pair <- arrayInd(pair_idx, .dim = c(ploidy + 1, ploidy + 1)) - 1
    
    p1_sum <- pair[1]
    p2_sum <- pair[2]
    
    # Generate rows for P1 and P2
    p1_row <- c(rep(1, p1_sum), rep(0, ploidy - p1_sum))
    p2_row <- c(rep(1, p2_sum), rep(0, ploidy - p2_sum))
    
    P1[i,] <- sample(p1_row)
    P2[i,] <- sample(p2_row)
  }
  
      

  P1 <- ph_matrix_to_list(P1)
  P2 <- ph_matrix_to_list(P2)
  p <- unlist(lapply(P1, function(x) sum(as.logical(x))))
  q <- unlist(lapply(P2, function(x) sum(as.logical(x))))
  return(list(hom.allele.p = P1,
              hom.allele.q = P2,
              p = p, q = q, 
              ploidy = ploidy))
}
