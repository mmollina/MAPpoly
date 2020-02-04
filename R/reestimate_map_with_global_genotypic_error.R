#' Introduces genotyping error into a full informative vector
#'
#' If \code{restricted = TRUE}, it restricts the prior to the 
#' possible classes under mendelian non double-reduced segregation 
#' given dosage of the parents
#'
#' @param void interfunction to be documented
#' @keywords internal
#'
#' @export genotyping_global_error
#'
genotyping_global_error<-function(x, m, restricted = TRUE,  error=0.01, th.prob=0.95)
{
  if(restricted){
    x1 <- x[1:(m+1)]
    if(sum(x1 > th.prob) == 1){
      x2 <- x[m + 2:3]
      id<-segreg_poly(m, dP = x2[1], dQ = x2[2]) > 0
      x3 <- x1[id]
      o<-which.max(x3)
      x3[o]<-1-error
      x3[-o]<-error/(sum(id)-1)
      x1[match(names(x3), names(x1))] <- x3
      return(x1)
    }
    return(x1)
  } else {
    x1 <- x[1:(m+1)]
    if(sum(x1 > th.prob)==1){
      o<-which.max(x1)
      x1[o]<-1-error
      x1[-o]<-error/(length(x1)-1)
      return(x1)
    }
    return(x1)
  }
}

#' Reestimate genetic map given a global genotyping error
#'
#' This function considers a global error when reestimating
#' a genetic map using Hidden Markov models
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @param error the assumed global error rate (default = NULL)
#' 
#' @param tol the desired accuracy (default = 10e-04)
#' 
#' @param restricted if \code{TRUE} (default), restricts the prior to the 
#'                   possible classes under mendelian non double-reduced segregation 
#'                   given dosage of the parents
#'                    
#' @param th.prob the threshold for using global error or genotype 
#'     probability distribution contained in the data set (default = 0.95)
#'      
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE} (default), no output is produced
#'
#' @return An object of class 'mappoly.map' with the following structure:
#' \item{m}{the ploidy level}
#' \item{mrk.names}{the names of markers present in the sequence}
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' \item{maps}{a list containing the sequence of markers, their recombination fractions,
#' the linkage phase configuration for all markers in both parents P and Q and the 
#' map's joint likelihood}
#'
#' @examples
#'   \dontrun{
#'     seq1.20<-make_seq_mappoly(hexafake, 1:20)
#'     counts<-cache_counts_twopt(seq1.20,
#'                                get.from.web=TRUE)
#'     subset.pairs<-est_pairwise_rf(seq1.20, counts,
#'                                n.clusters=1)
#'     subset.map <- est_rf_hmm_sequential(input.seq  = seq1.20,
#'                                         thres.twopt = 5,
#'                                         thres.hmm = 10,
#'                                         extend.tail = 10,
#'                                         tol = 0.1,
#'                                         tol.final = 10e-4,
#'                                         twopt = subset.pairs,
#'                                         verbose = TRUE)
#'     subset.map
#'     plot(subset.map)      
#'     subset.map.reest<-est_full_hmm_with_global_error(subset.map, 
#'                                                      error=0.01, 
#'                                                      tol=10e-4, 
#'                                                      verbose = TRUE)
#'     subset.map.reest
#'     plot(subset.map.reest)
#'     }
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
#' @export est_full_hmm_with_global_error
#'
est_full_hmm_with_global_error <- function(input.map, error=NULL, tol=10e-4, 
                                           restricted = TRUE, 
                                           th.prob=0.95, verbose = FALSE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  output.seq<-input.map
  mrknames<-get(input.map$info$data.name, pos=1)$mrk.names[input.map$maps[[1]]$seq.num]
  ## 
  if(nrow(get(input.map$info$data.name, pos=1)$geno)==get(input.map$info$data.name, pos=1)$n.mrk){
    geno.temp<-get(input.map$info$data.name, pos=1)$geno.dose[mrknames,]
    indnames<-get(input.map$info$data.name, pos=1)$ind.names
    gen<-vector("list", length(indnames))
    names(gen)<-indnames
    mrk<-ind<-NULL
    dp<-get(input.map$info$data.name, pos=1)$dosage.p[input.map$maps[[1]]$seq.num]
    dq<-get(input.map$info$data.name, pos=1)$dosage.q[input.map$maps[[1]]$seq.num]
    names(dp)<-names(dq)<-mrknames
    d.pq<-data.frame(dp = dp, 
                     dq = dq)
    d.pq$mrk<-mrknames
    for(i in names(gen))
    {
      a<-matrix(0, nrow(geno.temp), input.map$info$m+1, dimnames = list(mrknames, 0:input.map$info$m))
      for(j in rownames(a)){
        if(geno.temp[j,i] == input.map$info$m+1){
          a[j,]<-segreg_poly(m = input.map$info$m, dP = dp[j], dQ = dq[j])
        } else {
          a[j,geno.temp[j,i]+1]<-1          
        }
      }
      a<-as.data.frame(a)
      a$mrk<-rownames(a)
      a.temp<-t(merge(a, d.pq, sort = FALSE)[,-c(1)])
      if(!is.null(error))
        a.temp<-apply(a.temp, 2, genotyping_global_error, m = input.map$info$m, 
                      restricted = restricted, error=error, th.prob = th.prob)
      else
        a.temp <- a.temp[1:(input.map$info$m+1), ]
      colnames(a.temp)<-a[,1]
      gen[[i]]<-a.temp
    }
  } 
  else {##
    geno.temp<-subset(get(input.map$info$data.name, pos=1)$geno, mrk%in%mrknames)
    indnames<-get(input.map$info$data.name, pos=1)$ind.names
    gen<-vector("list", length(indnames))
    names(gen)<-indnames
    mrk<-ind<-NULL
    d.pq<-data.frame(dp = get(input.map$info$data.name, pos=1)$dosage.p[input.map$maps[[1]]$seq.num], 
                     dq = get(input.map$info$data.name, pos=1)$dosage.q[input.map$maps[[1]]$seq.num])
    d.pq$mrk<-mrknames
    for(i in names(gen))
    {
      a<-subset(geno.temp, ind%in%i)
      a<-a[match(mrknames, a$mrk),]
      a.temp<-t(merge(a, d.pq, sort = FALSE)[,-c(1:2)])
      if(!is.null(error))
        a.temp<-apply(a.temp, 2, genotyping_global_error, m = input.map$info$m, 
                      restricted = restricted, error=error, th.prob = th.prob)
      else
        a.temp <- a.temp[1:(input.map$info$m+1), ]
      colnames(a.temp)<-a[,1]
      gen[[i]]<-a.temp
    }
  }
  cat("
 ----------------------------------------------
 INFO: running HMM using full transition space:
       this operation may take a while.
-----------------------------------------------\n")
  for(i in 1:length(input.map$maps))
  {
    YP<-input.map$maps[[i]]$seq.ph$P
    YQ<-input.map$maps[[i]]$seq.ph$Q
    map<-poly_hmm_est(m = as.numeric(input.map$info$m),
                      n.mrk = as.numeric(input.map$info$n.mrk),
                      n.ind = as.numeric(length(gen)),
                      p = as.numeric(unlist(YP)),
                      dp = as.numeric(cumsum(c(0, sapply(YP, function(x) sum(length(x)))))),
                      q = as.numeric(unlist(YQ)),
                      dq = as.numeric(cumsum(c(0, sapply(YQ, function(x) sum(length(x)))))),
                      g = as.double(unlist(gen)),
                      rf = as.double(input.map$maps[[i]]$seq.rf),
                      verbose = verbose,
                      tol = tol)
    output.seq$maps[[i]]$seq.rf<-map$rf
    output.seq$maps[[i]]$loglike<-map$loglike
  }
  return(output.seq)
}
