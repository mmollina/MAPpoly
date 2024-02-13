#' Multilocus analysis using Hidden Markov Models (single parent, single phase)
#' @keywords internal
est_rf_hmm_single_phase_single_parent <- function(input.seq,
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
    stop("Incorporation of global error not available yet", call. = FALSE)
  ploidy <- input.seq$ploidy
  P <- do.call(cbind, input.seq[grep(pattern = "seq.dose.p", names(input.seq))])
  if(!all(unique(P[,uninfo.parent]) %in% c(0, input.seq$ploidy)))
    stop("wrong uninformative parent", call. = FALSE)
  D <- get(input.seq$data.name, pos = 1)$geno.dose[input.seq$seq.num,]
  seq.num <- input.seq$seq.num
  n.mrk <- nrow(D)
  n.ind <- ncol(D)
  ph <- input.ph.single[[info.parent]]
  h <- get_states_and_emission_single_parent(ploidy, ph, global.err, D, P[,uninfo.parent])
  if(is.null(rf.vec))
    rf.vec <- rep(0.001, n.mrk-1)
  res.temp  <- 
    .Call("est_hmm_map_single_parent",
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
  # return(structure(list(info = list(ploidy = ploidy,
  #                                   n.mrk = sum(id),
  #                                   seq.num = input.seq$seq.num,
  #                                   mrk.names = input.seq$seq.mrk.names,
  #                                   seq.dose.p1 = input.seq$seq.dose.p1[input.seq$seq.mrk.names],
  #                                   seq.dose.p2 = input.seq$seq.dose.p2[input.seq$seq.mrk.names],
  #                                   chrom = input.seq$chrom[input.seq$seq.mrk.names],
  #                                   genome.pos = input.seq$genome.pos[input.seq$seq.mrk.names],
  #                                   seq.ref = input.seq$seq.ref[input.seq$seq.mrk.names],
  #                                   seq.alt = input.seq$seq.alt[input.seq$seq.mrk.names],
  #                                   chisq.pval = input.seq$chisq.pval[input.seq$seq.mrk.names],
  #                                   data.name = input.seq$data.name,
  #                                   ph.thresh = input.seq$chisq.pval.thres),
  #                       maps = list(list(seq.num = seq.num,
  #                                        seq.rf = res.temp[[2]],
  #                                        seq.ph = list(P = input.ph.single[[1]],
  #                                                      Q = input.ph.single[[2]]),
  #                                        loglike = res.temp[[1]]))),
  #                  class = "mappoly.map"))
  map <- list(seq.num = seq.num,
               seq.rf = res.temp[[2]],
               seq.ph = list(P = input.ph.single[[1]],
                             Q = input.ph.single[[2]]),
               loglike = res.temp[[1]])
  map
}
