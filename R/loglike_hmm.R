#' Multipoint log-likelihood computation
#' 
#' Update the multipoint log-likelihood of a given map 
#' using the method proposed by \cite{Mollinari and Garcia (2019)}.
#' 
#' @param input.map An object of class \code{mappoly.map}
#' 
#' @param input.data  An object of class \code{mappoly.data}, which was used to generate \code{input.map}
#' 
#' @param verbose If \code{TRUE}, map information is shown; if
#'     \code{FALSE}(default), no output is produced
#'     
#' @examples
#'  \donttest{
#'   hexa.map<-loglike_hmm(maps.hexafake[[1]])
#'   hexa.map
#'  }
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
#' @export loglike_hmm
#' 
loglike_hmm<-function(input.map, input.data = NULL, verbose = FALSE)
{
  ## Checking class of arguments
  if(!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  if(is.null(input.data))
    D<-get(input.map$info$data.name, pos=1)$geno.dose[input.map$info$mrk.names,]
  else
    D<-input.data$geno.dose[input.map$info$mrk.names,]
  dp <- input.map$info$seq.dose.p
  dq <- input.map$info$seq.dose.q
  for (j in 1:nrow(D))
    D[j, D[j, ] == input.map$info$m + 1] <- dp[j] + dq[j] + 1 + as.numeric(dp[j]==0 || dq[j]==0)
  for(i in 1:length(input.map$maps)){
    rf.temp<-input.map$maps[[i]]$seq.rf
    res.temp<-.Call("loglike_hmm",
                    input.map$info$m,
                    t(D),
                    lapply(input.map$maps[[i]]$seq.ph$P, function(x) x-1),
                    lapply(input.map$maps[[i]]$seq.ph$Q, function(x) x-1),
                    rf.temp,
                    verbose,
                    PACKAGE = "mappoly")
    input.map$maps[[i]]$loglike <- res.temp[[1]]
  }
  return(input.map)
}
