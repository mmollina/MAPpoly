#' Introduces genotyping error into a full informative vector
#'
#' @param void interfunction to be documented
#' @keywords internal
#'
#' @export genotyping_global_error
#'
genotyping_global_error<-function(x, error=0.01, th.prob=0.999)
{
  if(sum(x > th.prob)==1){
    o<-which.max(x)
    x[o]<-1-error
    x[-o]<-error/(length(x)-1)
    return(x)
  }
  return(x)
}

#' Reestimate gentic map given a global genotyping error
#'
#' This function considers a global error when reestimating
#' a genetic map using hdden Markov models
#'
#' @param input.map an object of class \code{mappoly.map}.
#' @param error global error rate
#' @param tol the desired accuracy.
#' @param th.prob the threshold for using global error or genotype 
#'     probability distribution contained in the data set 
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#'
#' @return An object of class 'mappoly.map'
#'
#' @examples
#'   \dontrun{
#'     hexa.file<-system.file('extdata', 'hexafake', package = 'mappoly')
#'     hexa.dat<-read_geno(file.in = hexa.file)
#'     seq1.20<-make_seq_mappoly(hexa.dat, 1:20)
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
#'      
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
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @export est_full_hmm_with_global_error
#'
est_full_hmm_with_global_error <- function(input.map, error=NULL, tol=10e-4, 
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
      a.temp<-t(a)
      if(!is.null(error))
        a.temp<-apply(a.temp, 2, genotyping_global_error, error=error, th.prob = 0.9)
      gen[[i]]<-a.temp
    }
  } 
  else {##
    geno.temp<-subset(get(input.map$info$data.name, pos=1)$geno, mrk%in%mrknames)
    indnames<-get(input.map$info$data.name, pos=1)$ind.names
    gen<-vector("list", length(indnames))
    names(gen)<-indnames
    mrk<-ind<-NULL
    for(i in names(gen))
    {
      a<-subset(geno.temp, ind%in%i)
      a<-a[match(mrknames, a$mrk),]
      a.temp<-t(a[,-c(1:2)])
      if(!is.null(error))
        a.temp<-apply(a.temp, 2, genotyping_global_error, error=error, th.prob = th.prob)
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
