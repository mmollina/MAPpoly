#' Compute conditional probabilities of the genotypes using global error
#'
#' Conditional genotype probabilities are calculated for each marker
#' position and each individual given a map. 
#'
#' @param input.map An object of class \code{mappoly.map}
#' 
#' @param step 	Maximum distance (in cM) between positions at which 
#'              the genotype probabilities are calculated, though for 
#'              step = 0, probabilities are calculated only at the 
#'              marker locations.
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#'    
#' @param error the assumed global error rate (default = 0.01)
#' 
#' @param th.prob the threshold for using global error or genotype 
#'     probability distribution contained in the dataset (default = 0.95)
#'     
#' @param restricted if \code{TRUE} (default), restricts the prior to the 
#'                   possible classes under Mendelian non double-reduced 
#'                   segregation given the parental dosages 
#'                   
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @return An object of class 'mappoly.genoprob' which has two elements: a tridimensional
#' array containing the probabilities of all possible genotypes for each individual
#' in each marker position; and the marker sequence with it's recombination frequencies
#' 
#' @examples
#'  \donttest{
#'      probs.error<-calc_genoprob_error(input.map = solcap.err.map[[1]],
#'                                 error = 0.05,
#'                                 verbose = TRUE)
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
#' @export calc_genoprob_error
calc_genoprob_error<-function(input.map,  step = 0, phase.config = "best", error = 0.01, 
                              th.prob = 0.95, restricted = TRUE, verbose = TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  if (verbose && !capabilities("long.double")){
    cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  m<-input.map$info$m
  n.ind<-get(input.map$info$data.name, pos=1)$n.ind
  Dtemp<-get(input.map$info$data.name, pos=1)$geno.dose[input.map$maps[[1]]$seq.num,]
  map.pseudo <- create_map(input.map, step, phase.config = i.lpc)
  mrknames<-names(map.pseudo)
  n.mrk<-length(map.pseudo)
  indnames<-colnames(Dtemp)
  D<-matrix(m+1, nrow = length(map.pseudo), ncol = ncol(Dtemp), 
            dimnames = list(mrknames, indnames))
  D[rownames(Dtemp), ] <- as.matrix(Dtemp)
  dptemp <- get(input.map$info$data.name)$dosage.p[input.map$maps[[i.lpc]]$seq.num]
  dqtemp <- get(input.map$info$data.name)$dosage.q[input.map$maps[[i.lpc]]$seq.num]
  dq<-dp<-rep(m/2, length(mrknames))
  names(dp)<-names(dq)<-mrknames
  dp[names(dptemp)]<-dptemp
  dq[names(dqtemp)]<-dqtemp
 
  phP <- phQ <- vector("list", n.mrk)
  for(i in 1:length(phP)){
    phP[[i]] <- phQ[[i]] <- c(0:(m/2 - 1))
  }
  names(phP) <- names(phQ) <- mrknames
  phP[rownames(Dtemp)] <- input.map$maps[[i.lpc]]$seq.ph$P
  phQ[rownames(Dtemp)] <- input.map$maps[[i.lpc]]$seq.ph$Q

  ## Including error  
  gen<-vector("list", length(indnames))
  names(gen)<-indnames
  d.pq<-data.frame(dp = dp, 
                   dq = dq)
  d.pq$mrk<-rownames(d.pq)
  for(i in names(gen))
  {
    a<-matrix(0, nrow(D), input.map$info$m+1, dimnames = list(mrknames, 0:input.map$info$m))
    for(j in rownames(a)){
      if(D[j,i] == input.map$info$m+1){
        a[j,]<-segreg_poly(m = input.map$info$m, dP = dp[j], dQ = dq[j])
      } else {
        a[j,D[j,i]+1]<-1          
      }
    }
    a<-as.data.frame(a)
    a$mrk<-rownames(a)
    a.temp<-t(merge(a, d.pq, sort = FALSE)[,-c(1)])
    if(!is.null(error))
      a.temp<-apply(a.temp, 2, genotyping_global_error, 
                    m = input.map$info$m, error=error, 
                    th.prob = th.prob, 
                    restricted = restricted)
    else
      a.temp <- a.temp[1:(input.map$info$m+1), ]
    colnames(a.temp)<-a[,1]
    gen[[i]]<-a.temp
  }
  g = as.double(unlist(gen))
  p = as.numeric(unlist(phP))
  dp = as.numeric(cumsum(c(0, sapply(phP, function(x) sum(length(x))))))
  q = as.numeric(unlist(phQ))
  dq = as.numeric(cumsum(c(0, sapply(phQ, function(x) sum(length(x))))))
  rf = as.double(mf_h(diff(map.pseudo)))
  res.temp <-
    .Call(
      "calc_genoprob_prior",
      as.numeric(m),
      as.numeric(n.mrk),
      as.numeric(n.ind),
      as.numeric(p),
      as.numeric(dp),
      as.numeric(q),
      as.numeric(dq),
      as.double(g),
      as.double(rf),
      as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
      as.double(0),
      as.numeric(verbose),
      PACKAGE = "mappoly"
    )
  if (verbose) cat("\n")
  dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                          apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"),
                                mrknames, indnames)
  structure(list(probs = res.temp[[1]], map = map.pseudo), class="mappoly.genoprob")
}
