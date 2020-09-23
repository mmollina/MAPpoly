#' Simulate homology groups
#'
#' Simulate two homology groups (one for each parent) and their
#' linkage phase configuration.
#'
#' This function prevents the simulation of linkage phase
#' configurations which are impossible to estimate via two point
#' methods
#'
#' @param m ploidy level. Must be an even number
#' 
#' @param n.mrk number of markers
#' 
#' @param min.d minimum dosage to be simulated (default = 0)
#' 
#' @param max.d maximum dosage to be simulated (default = m + 1)
#' 
#' @param prob.dose a vector indicating the proportion of markers for
#'    different dosage to be simulated (default = NULL)
#'    
#' @param max.ph maximum phase difference
#'  
#' @param restriction if TRUE (default), avoid cases where it is impossible to
#'     estimate recombination fraction and/or linkage phases via
#'     two-point analysis
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
#'
#' @examples
#'     h.temp<-sim_homologous(m=6, n.mrk=20, max.d=3, max.ph=3,
#'                            seed=123)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'
#' @export
sim_homologous<-function(m, n.mrk, min.d = 0, 
                         max.d = m+1, prob.dose = NULL,
                         max.ph, restriction=TRUE, 
                         seed=NULL)
{
  #prob.dose<-NULL
  if(!is.null(seed)) set.seed(seed)
  
  hom.allele.q<-hom.allele.p<-vector("list", n.mrk)
  count<-1
  while(count <= n.mrk)
  {
    hom.p.temp<-hom.q.temp<-0
    if(any(is.null(prob.dose)))
      p.temp<-sample(min.d:max.d,1)
    else
      p.temp<-sample(min.d:max.d,1, prob = prob.dose)
    if(all(p.temp!=0))
      hom.p.temp<-sample(1:m, p.temp)
    if(any(is.null(prob.dose)))
      q.temp<-sample(min.d:max.d,1)
    else
      q.temp<-sample(min.d:max.d,1, prob = prob.dose)
    if(all(q.temp!=0))
      hom.q.temp<-sample(1:m, q.temp)
    p <- sum(as.logical(hom.p.temp))
    q <- sum(as.logical(hom.q.temp))
    if(restriction && count > 1)
    {
      if(!any((p+q)==0,
              (p+q)==2*m,
              sum(as.logical(hom.allele.p[[count-1]])) - q == 0,
              sum(as.logical(hom.allele.q[[count-1]])) - p == 0,
              (p==m & q==0),
              (p==0 & q==m)))
      {
        hom.allele.p[[count]] <- hom.p.temp
        hom.allele.q[[count]] <- hom.q.temp
        count<-count+1
      }
    }
    else
    {
      if(!any((p+q)==0,
              (p+q)==2*m,
              (p==m & q==0),
              (p==0 & q==m)))
      {
        hom.allele.p[[count]] <- hom.p.temp
        hom.allele.q[[count]] <- hom.q.temp
        count<-count+1
      }
    }
    p<-unlist(lapply(hom.allele.p, function(x) sum(as.logical(x))))
    q<-unlist(lapply(hom.allele.q, function(x) sum(as.logical(x))))
  }
  return(list(hom.allele.p=hom.allele.p,
              hom.allele.q=hom.allele.q,
              p=p, q=q))
}
