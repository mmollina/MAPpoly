#' Introduces genotyping error into a full informative vector
#'
#' @param void interfunction to be documented
#' @keywords internal
#'
#' @export genotyping_global_error
#'
genotyping_global_error<-function(x, error=0.01, th.num=0.999)
{
  if(sum(x > th.num)==1){
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
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#'
#' @return An object of class 'mappoly.map'
#'
#' @examples
#'   \dontrun{
#'     hexa_file<-system.file('extdata', 'hexa_fake', package = 'mappoly2')
#'     hexa_dat<-read_geno(file_in = hexa_file)
#'     all_mrk<-make_seq_mappoly(hexa_dat, 'all')
#'     counts_all_mrk_from_web<-cache_counts_twopt(all_mrk,
#'                                                 get.from.web=TRUE)
#'     all_pairs<-est_pairwise_rf(all_mrk, counts_all_mrk_from_web,
#'                                n.clusters=2)
#'
#'     seq1.10<-make_seq_mappoly(hexa_dat, 1:10)
#'
#'     l_seq1_3.0 <- ls_linkage_phases(input.seq = seq1.10, thres = 1.3,
#'                                     twopt = all_pairs)
#'     l_seq1_3.0
#'     plot(l_seq1_3.0)
#'     res<-est_rf_hmm(input.seq=seq1.10, input.ph=l_seq1_3.0,
#'                     verbose=TRUE, tol=10e-3)
#'     res
#'     plot(res)
#'
#'     res_final<-est_full_hmm_with_global_error(res, error=0.01, tol=10e-4)
#'     plot(res_final)
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
est_full_hmm_with_global_error <- function(input.map, error=NULL, tol=10e-4, verbose = FALSE)
  {
  output.seq<-input.map
  mrknames<-get(input.map$info$data.name, pos=1)$mrk.names[input.map$maps[[1]]$seq.num]
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
      a.temp<-apply(a.temp, 2, genotyping_global_error, error=error)
    colnames(a.temp)<-a[,1]
    gen[[i]]<-a.temp
  }
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
