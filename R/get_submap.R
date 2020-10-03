#' Extract sub-map from map
#'
#' Given a pre-constructed map, it extracts a sub-map for a provided 
#' sequence of marker positions. Optionally, it can update the linkage phase
#' configurations and respective recombination fractions.
#'
#' @param input.map An object of class \code{mappoly.map}
#' 
#' @param mrk.pos positions of the markers that should be considered in the new map. 
#' This can be in any order
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the configuration associated with the maximum 
#'                     likelihood
#'                     
#' @param reestimate.rf logical. If \code{TRUE} (default) the recombination fractions 
#' between markers are re-estimated
#' 
#' @param reestimate.phase logical. If \code{TRUE}, the linkage phase configurations are 
#' re-estimated (default = FALSE)
#' 
#' @param thres.twopt the LOD threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered (default = 5)
#'     
#' @param thres.hmm the threshold used to determine if the linkage
#'     phases compared via hmm analysis should be considered (default = 3)
#'     
#' @param extend.tail the length of the tail of the chain that should
#'     be used to calculate the likelihood of the linkage phases. If
#'     \code{info.tail = TRUE}, the function uses at least \code{extend.tail}
#'     as the length of the tail (default = 50)
#'     
#' @param tol the desired accuracy during the sequential phase (default = 0.1)
#' 
#' @param tol.final the desired accuracy for the final map (default = 10e-04)
#' 
#' @param use.high.precision logical. If \code{TRUE} uses high precision 
#' (long double) numbers in the HMM procedure implemented in C++,
#' which can take a long time to perform (default = FALSE)
#'     
#' @param verbose If \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @examples
#'  \donttest{
#'     ## selecting the six first markers in linkage group 1
#'     ## re-estimating the recombination fractions and linkage phases
#'     submap1.lg1<-get_submap(input.map = maps.hexafake[[1]], 
#'                            mrk.pos = 1:6, verbose = TRUE,
#'                            reestimate.phase = TRUE,  
#'                            tol.final = 10e-3)
#'    ## no recombination fraction re-estimation: first 20 markers
#'    submap2.lg1<-get_submap(input.map = maps.hexafake[[1]], 
#'                            mrk.pos = 1:20, reestimate.rf = FALSE,
#'                            verbose = TRUE, 
#'                            tol.final = 10e-3)
#'   plot(maps.hexafake[[1]])
#'   plot(submap1.lg1, mrk.names = TRUE, cex = .8)
#'   plot(submap2.lg1, mrk.names = TRUE, cex = .8)
#'   }
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
#' @export get_submap
#'
get_submap<-function(input.map, mrk.pos,  phase.config = "best", reestimate.rf = TRUE,
                     reestimate.phase = FALSE, thres.twopt = 5, thres.hmm = 3,
                     extend.tail = 50, tol = 0.1, tol.final = 10e-4,
                     use.high.precision = FALSE, verbose=TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  high.prec <- use.high.precision
  if(!use.high.precision & length(mrk.pos) > 1000)
  {
    high.prec <- TRUE
    if (verbose) message("Number of markers is greater than 1000: \nusing high precision estimation")
  }
 ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } 
  else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } 
  else i.lpc <- phase.config
  input.obj = get(input.map$info$data.name, pos = 1)
  if(reestimate.rf & !reestimate.phase)
  {
    seq.num<-input.map$maps[[i.lpc]]$seq.num[mrk.pos]
    ph<-list(P = input.map$maps[[i.lpc]]$seq.ph$P[as.character(seq.num)],
             Q = input.map$maps[[i.lpc]]$seq.ph$Q[as.character(seq.num)])
    s<-make_seq_mappoly(input.obj = input.obj,
                        arg = seq.num, 
                        data.name = input.map$info$data.name)
    rf.temp <- NULL#input.map$maps[[i.lpc]]$seq.rf[mrk.pos[-length(mrk.pos)]]
    res<-est_rf_hmm_single(input.seq = s, input.ph.single = ph, tol = tol.final,
                           verbose = verbose, rf.temp = rf.temp,
                           high.prec = high.prec)
    output.map<-input.map
    output.map$maps[[i.lpc]] <- res
    output.map$info$n.mrk <- length(mrk.pos)
    output.map$info$mrk.names <- input.map$info$mrk.names[mrk.pos]
    return(output.map)
  }
  else if(reestimate.phase)
  {
    if(verbose && !reestimate.rf)
      message("
    The recombination fraction will be reestimated 
    since 'reestimate.phase = TRUE'")
    s<-make_seq_mappoly(get(input.map$info$data.name, pos = 1),
                        input.map$maps[[i.lpc]]$seq.num[mrk.pos],
                        input.map$info$data.name)
    if(verbose)
      cat("\nEstimating pairwise recombination fraction for marker sequence ...")
    p<-est_pairwise_rf(input.seq = s, verbose = FALSE)
    if(verbose){
      cat("done\n")      
      cat("\nEstimating sequential map: \n----------------------------------------\n")
    }
    output.map<-est_rf_hmm_sequential(input.seq = s,
                                      thres.twopt = thres.twopt,
                                      thres.hmm = thres.hmm,
                                      extend.tail = extend.tail,
                                      twopt = p,
                                      verbose = verbose,
                                      tol = tol,
                                      tol.final = tol.final, 
                                      high.prec = high.prec)
    return(output.map)
  }
  output.map <- input.map
  z<-cumsum(c(0, imf_h(input.map$maps[[i.lpc]]$seq.rf)))[mrk.pos]
  rf.vec<-mf_h(abs(diff(z)))
  if (verbose) message("
    You selected: reestimate.rf = FALSE
    -----------------------------------------
    The recombination fractions provided were
    obtained using the marker positions in the 
    input map; For accurate values, plese 
    reestimate the map using functions 'reest_rf', 
    'est_full_hmm_with_global_error' or 
    'est_full_hmm_with_prior_prob'")
  ##phase info
  output.map$maps[[i.lpc]]$seq.rf <- rf.vec
  output.map$maps[[i.lpc]]$seq.num <- input.map$maps[[i.lpc]]$seq.num[mrk.pos]
  output.map$maps[[i.lpc]]$seq.ph$P <- input.map$maps[[i.lpc]]$seq.ph$P[mrk.pos]
  output.map$maps[[i.lpc]]$seq.ph$Q <- input.map$maps[[i.lpc]]$seq.ph$Q[mrk.pos]
  output.map$maps[[i.lpc]]$loglike <- 0
  ##map info
  output.map$info$n.mrk <- length(mrk.pos)
  output.map$info$mrk.names <- input.map$info$mrk.names[mrk.pos]
  output.map$info$seq.num <- input.map$info$seq.num [mrk.pos]
  output.map$info$mrk.names <- input.map$info$mrk.names[mrk.pos]
  output.map$info$seq.dose.p <- input.map$info$seq.dose.p[mrk.pos]
  output.map$info$seq.dose.q <- input.map$info$seq.dose.q[mrk.pos]
  output.map$info$sequence <- input.map$info$sequence[mrk.pos] 
  output.map$info$sequence.pos <- input.map$info$sequence.pos[mrk.pos]  
  output.map$info$seq.ref <- input.map$info$seq.ref[mrk.pos] 
  output.map$info$seq.alt <- input.map$info$seq.alt[mrk.pos] 
  output.map$info$chisq.pval <- input.map$info$chisq.pval[mrk.pos] 
  
  return(output.map)
}
