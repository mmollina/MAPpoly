#' Multilocus analysis using Hidden Markov Models (single parent, single phase)
#'
#' @param void internal function to be documented
#' @keywords internal
#' @examples
#'   \donttest{
#'     s <- make_seq_mappoly(solcap.dose.map[[1]])
#'     full.phase <- solcap.dose.map[[1]]$maps[[1]]$seq.ph
#'     P1.map <- est_rf_hmm_single_one_parent(s, full.phase, info.parent = 1, 
#'                                            uninfo.parent = 2)
#'     plot(P1.map)
#'     P2.map <- est_rf_hmm_single_one_parent(s, full.phase, info.parent = 2, 
#'                                            uninfo.parent = 1)
#'     plot(P2.map)
#'     plot_map_list(list(Atlantic = P1.map, B1829 = P2.map), , horiz = FALSE)
#'  }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export est_rf_hmm_single_one_parent
#' 
est_rf_hmm_single_one_parent <- function(input.seq,
                                         input.ph.single,
                                         info.parent = 1,
                                         uninfo.parent = 2,
                                         rf.vec = NULL,
                                         global.err = 0.0,
                                         tol = 10e-4,
                                         verbose = FALSE,
                                         ret.map.no.rf.estimation = FALSE)
{
  input_classes <- c("mappoly.sequence")
  if (!inherits(input.seq, input_classes[1])) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  if(length(input.seq$seq.num)  ==  1)
    stop("Input sequence contains only one marker.", call. = FALSE)
  if(global.err != 0)
    stop("Incorporation of global error not available yet")
  ploidy <- input.seq$ploidy
  P <- do.call(cbind, input.seq[grep(pattern = "seq.dose.p", names(input.seq))])
  id <- P[,info.parent] != 0 & apply(P[,uninfo.parent, drop = FALSE], 1, function(x) all(x==0))
  D <- get(input.seq$data.name, pos = 1)$geno.dose[input.seq$seq.num[id],]
  seq.num <- input.seq$seq.num[id]
  n.mrk <- nrow(D)
  n.ind <- ncol(D)
  ph <- input.ph.single[[info.parent]][id]
  h <- get_states_and_emission_one_parent(ploidy, ph, global.err, D)
  if(is.null(rf.vec))
    rf.vec <- rep(0.001, n.mrk-1)
  res.temp  <- 
    .Call("est_hmm_map_one_parent",
          ploidy,
          n.mrk,
          n.ind,
          h$states,
          h$emission,
          rf.vec,
          verbose,
          tol,
          ret.map.no.rf.estimation,
          PACKAGE = "mappoly")
  return(structure(list(info = list(ploidy = ploidy,
                                    n.mrk = sum(id),
                                    seq.num = input.seq$seq.num[id],
                                    mrk.names = input.seq$seq.mrk.names[id],
                                    seq.dose.p1 = input.seq$seq.dose.p1[input.seq$seq.mrk.names[id]],
                                    seq.dose.p2 = input.seq$seq.dose.p2[input.seq$seq.mrk.names[id]],
                                    chrom = input.seq$chrom[input.seq$seq.mrk.names[id]],
                                    genome.pos = input.seq$genome.pos[input.seq$seq.mrk.names[id]],
                                    seq.ref = input.seq$seq.ref[input.seq$seq.mrk.names[id]],
                                    seq.alt = input.seq$seq.alt[input.seq$seq.mrk.names[id]],
                                    chisq.pval = input.seq$chisq.pval[input.seq$seq.mrk.names[id]],
                                    data.name = input.seq$data.name,
                                    ph.thresh = input.seq$chisq.pval.thres),
                        maps = list(list(seq.num = seq.num,
                                         seq.rf = res.temp[[2]],
                                         seq.ph = list(P = input.ph.single[[1]][id],
                                                       Q = input.ph.single[[2]][id]),
                                         loglike = res.temp[[1]]))),
                   class = "mappoly.map"))
}
