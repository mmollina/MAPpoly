#' Extract sub-map from map
#'
#' Given a pre-constructed map, it extracts a sub-map for a sequence of marker positions.
#'
#' @param input.map An object of class \code{mappoly.map}
#' @param mrk.pos positions of the markers that should be considered in the new map. This can be in any order.
#' @param phase.config which phase configuration should be used
#'    "best" will choose the one with highest likelihood
#' @param reestimate.rf logical. If \code{TRUE} (default) the recombination fractions between markers are reestimated.
#' @param reestimate.phase logical. If \code{TRUE}, the linkage phase configuration is reestimated. The default is \code{FALSE}
#' @param thres.twopt the threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered
#' @param thres.hmm the threshold used to determine if the linkage
#'     phases compared via hmm analysis should be considered
#' @param extend.tail the length of the tail of the chain that should
#'     be used to calculate the likelihood of the linkage phases. If
#'     \code{info.tail = TRUE}, the function uses at least \code{extend.tail}
#'     as the lengthof the tail.
#' @param count.cache an object of class \code{cache.info} containing
#'     pre-computed genotype frequencies, obtained with
#'     \code{\link[mappoly]{cache_counts_twopt}}.
#' @param tol the desired accuracy during the sequential phase.
#' @param tol.final the desired accuracy for the final map.
#' @param use.high.precision If \code{TRUE}, uses a high precision double
#'     variable in C++ routine. This is advised for maps with a large number
#'     of markers.
#' @param verbose If \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#' @examples
#'  \dontrun{
#'     data(hexafake)
#'     ##100 first markers as an example
#'     s1<-make_seq_mappoly(hexafake, 1:100)
#'     red.mrk<-elim_redundant(s1)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     unique.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                   count.cache = counts.web,
#'                                   n.clusters = 1,
#'                                   verbose=TRUE)
#'     map <- est_rf_hmm_sequential(input.seq = unique.mrks,
#'                                         thres.twopt = 5,
#'                                         thres.hmm = 3,
#'                                         extend.tail = 50,
#'                                         tol = 0.1,
#'                                         tol.final = 10e-2,
#'                                         twopt = subset.pairs,
#'                                         verbose = TRUE)
#'   #MDS whithout reestimating phase
#'   m1<-rf_list_to_matrix(input.twopt = unique.pairs)
#'   mds.ord <- mds_mappoly(input.mat = m1, p = NULL, n = NULL, ndim = 2)
#'   mds.seq.ord <- make_seq_mappoly(mds.ord)
#'   mds.map<-get_submap(input.map = map, mrk.pos = order(mds.seq.ord$seq.num),
#'                       verbose = TRUE, tol.final = 10e-2)
#'   plot(map)
#'   plot(mds.map)
#'   map$maps[[1]]$loglike
#'   mds.map$maps[[1]]$loglike
#'   plot(map$maps[[1]]$seq.num, mds.map$maps[[1]]$seq.num)
#'   }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @export get_submap
#'
get_submap<-function(input.map, mrk.pos,  phase.config = "best", reestimate.rf = TRUE,
                     reestimate.phase = FALSE, thres.twopt = 5, thres.hmm = 3,
                     extend.tail = 50, count.cache = NULL, tol = 0.1, tol.final = 10e-4,
                     use.high.precision = TRUE, verbose=TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  if(reestimate.rf)
  {
    seq.num<-input.map$maps[[i.lpc]]$seq.num[mrk.pos]
    ph<-list(P = input.map$maps[[i.lpc]]$seq.ph$P[as.character(seq.num)],
             Q = input.map$maps[[i.lpc]]$seq.ph$Q[as.character(seq.num)])
    s<-make_seq_mappoly(input.obj = get(input.map$info$data.name, pos = 1),
                        arg = seq.num, data.name = input.map$info$data.name)
    rf.temp <- NULL#input.map$maps[[i.lpc]]$seq.rf[mrk.pos[-length(mrk.pos)]]
    res<-est_rf_hmm_single(input.seq = s, input.ph.single = ph, tol = tol.final,
                           verbose = verbose, rf.temp = rf.temp, high.prec = use.high.precision)
    output.map<-input.map
    output.map$maps[[i.lpc]]<-res
    output.map$info$n.mrk<-length(mrk.pos)
    return(output.map)
  }
  else if(reestimate.phase)
  {
    s<-make_seq_mappoly(get(input.map$info$data.name, pos = 1),
                        input.map$maps[[i.lpc]]$seq.num[mrk.pos],
                        input.map$info$data.name)
    p<-est_pairwise_rf(input.seq = s, count.cache = count.cache)
    output.map<-est_rf_hmm_sequential(input.seq = s,
                                      thres.twopt = thres.twopt,
                                      thres.hmm = thres.hmm,
                                      extend.tail = extend.tail,
                                      twopt = p,
                                      verbose = verbose,
                                      tol = tol,
                                      tol.final = tol.final)
    return(output.map)
  }
  output.map<-input.map
  output.map$maps[[i.lpc]]$seq.num <- input.map$maps[[i.lpc]]$seq.num[mrk.pos]
  output.map$maps[[i.lpc]]$seq.rf <- input.map$maps[[i.lpc]]$seq.rf[mrk.pos[-length(mrk.pos)]]
  output.map$maps[[i.lpc]]$seq.ph$P <- input.map$maps[[i.lpc]]$seq.ph$P[mrk.pos]
  output.map$maps[[i.lpc]]$seq.ph$Q <- input.map$maps[[i.lpc]]$seq.ph$Q[mrk.pos]
  output.map$info$n.mrk<-length(mrk.pos)
  return(output.map)
}
