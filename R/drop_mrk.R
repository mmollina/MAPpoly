#' Drop markers in a genetic map
#'
#' @param input.map An object of class \code{mappoly.map}
#' @param drop.mark Markers to be removed from the map
#' @param reest.map if \code{TRUE}, reestimate the map using the
#'     method indicated in the argument 'method'
#' @param input.mat An object of class \code{mappoly.rf.matrix}
#' @param tol tolerance for determining convergence.
#' @param phase.config should be a string \code{'best'} or the position of the
#'     configuration to be plotted. If \code{'best'}, plot the configuration
#'     with the highest likelihood.
#' @param method indicates whether to use Hidden Markov Models or Ordinary
#'     Least Squares to reestimate the recombination fraction
#' @param weight if \code{TRUE}, it uses the LOD scores to perform a weighted
#'    regression when the Ordinary Least Saquares is chosen
#' @param verbose If \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#' @examples
#'  \dontrun{
#'     data(hexafake)
#'     mrk.subset<-make_seq_mappoly(hexafake, 1:100)
#'     red.mrk<-elim_redundant(mrk.subset)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   count.cache = counts.web,
#'                                   n.clusters = 16,
#'                                   verbose=TRUE)
#'     system.time(
#'     subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
#'                                         thres.twopt = 5,
#'                                         thres.hmm = 3,
#'                                         extend.tail = 50,
#'                                         tol = 0.1,
#'                                         tol.final = 10e-3,
#'                                         twopt = subset.pairs,
#'                                         verbose = TRUE))
#'      plot(subset.map)
#'      #removing markers 'M10', 'M15' and 'M24'
#'      new.map<-drop_mrk(input.map = subset.map,
#'                        drop.mrk = c(10,15,24),
#'                        reest.map = TRUE,
#'                        tol = 10e-3,
#'                        method = "hmm")
#'      plot(new.map)
#'  }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2017) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_
#'
#' @export drop_mrk
#'
drop_mrk<-function(input.map, drop.mrk, reest.map = TRUE, input.mat = NULL,
                   tol = 10e-3,  phase.config = "best", method = c("hmm", "ols"),
                   weight = TRUE, verbose = TRUE)
{
  if (!inherits(input.map, "mappoly.map"))
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  output.map<-input.map
  ## indexing markers that should  be removed
  rem.rf <- rem <- match(drop.mrk, input.map$maps[[i.lpc]]$seq.num)
  if(any(rem == input.map$info$n.mrk))
    rem.rf <- rem.rf-1
  output.map$info$n.mrk<-output.map$info$n.mrk-length(drop.mrk)
  output.map$maps[[i.lpc]]$seq.num <- output.map$maps[[i.lpc]]$seq.num[-rem]
  output.map$maps[[i.lpc]]$seq.rf <- output.map$maps[[i.lpc]]$seq.rf[-rem.rf]
  output.map$maps[[i.lpc]]$seq.ph$P <- output.map$maps[[i.lpc]]$seq.ph$P[-rem]
  output.map$maps[[i.lpc]]$seq.ph$Q <- output.map$maps[[i.lpc]]$seq.ph$Q[-rem]
  output.map$maps[[i.lpc]]$loglike <- 0
  if(reest.map)
    output.map <- reest_map(input.map = output.map,
                            input.mat = input.mat,
                            tol = tol,
                            phase.config = phase.config,
                            method = method,
                            weight = weight,
                            verbose = verbose)
  return(output.map)
}

