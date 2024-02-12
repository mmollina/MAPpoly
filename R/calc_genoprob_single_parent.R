#' Compute conditional probabilities of the genotype (one informative parent)
#'
#' Conditional genotype probabilities are calculated for each marker
#' position and each individual given a map 
#'
#' @param input.map An object of class \code{mappoly.map} (with exceptions)
#' 
#' @param step 	Maximum distance (in cM) between positions at which
#'              the genotype probabilities are calculated, though for
#'              step = 0, probabilities are calculated only at the
#'              marker locations.
#' @param info.parent index for informative parent
#' 
#' @param uninfo.parent index for uninformative parent
#' 
#' @param global.err the assumed global error rate (default = 0.0)
#' 
#' @param phase.config which phase configuration should be used. "best" (default)
#'                     will choose the phase configuration associated with the
#'                     maximum likelihood
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @return An object of class 'mappoly.genoprob' which has two elements: a tridimensional
#' array containing the probabilities of all possible genotypes for each individual
#' in each marker position; and the marker sequence with it's recombination frequencies
#'
#' @examples
#'  ## tetraploid example
#'  s <- make_seq_mappoly(tetra.solcap, 'seq12', info.parent = "p1")
#'  tpt <- est_pairwise_rf(s)
#'  map <- est_rf_hmm_sequential(input.seq = s,
#'                               twopt = tpt,
#'                                start.set = 10,
#'                                thres.twopt = 10, 
#'                                thres.hmm = 10,
#'                                extend.tail = 4,
#'                                info.tail = TRUE, 
#'                                sub.map.size.diff.limit = 8, 
#'                                phase.number.limit = 4,
#'                                reestimate.single.ph.configuration = TRUE,
#'                                tol = 10e-2,
#'                                tol.final = 10e-3)
#'  plot(map)                                     
#'  probs <- calc_genoprob_single_parent(input.map = map, 
#'                                    info.parent = 1, 
#'                                    uninfo.parent = 2, 
#'                                    step = 1)
#'  probs
#'  ## displaying individual 1, 6 genotypic states
#'  ## (rows) across linkage group 1 (columns)                          
#'  image(t(probs$probs[,,2]))
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378} 
#'
#' @export calc_genoprob_single_parent
#' @importFrom dplyr left_join
calc_genoprob_single_parent <- function(input.map,
                                        step = 0,
                                        info.parent = 1,
                                        uninfo.parent = 2,
                                        global.err = 0.0,
                                        phase.config = "best", 
                                        verbose = TRUE)
{
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  if (verbose && !capabilities("long.double")){
    cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config  ==  "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  ploidy <- input.map$info$ploidy
  g <- choose(ploidy, ploidy/2)
  P <- do.call(cbind, input.map$info[grep(pattern = "seq.dose.p", names(input.map$info))])
  if(!all(unique(P[,uninfo.parent]) %in% c(0, ploidy)))
    stop("wrong uninformative parent")
  new.map <- create_map(input.map, step = step, phase.config = i.lpc)
  D <- get(input.map$info$data.name, pos = 1)$geno.dose
  D <- D[names(new.map),]
  D[is.na(D)] <- ploidy + 1
  rownames(D) <- names(new.map)
  n.mrk <- nrow(D)
  n.ind <- ncol(D)
  ph <- input.map$maps[[i.lpc]]$seq.ph[[info.parent]]
  names(ph) <- input.map$info$mrk.names
  ph <- ph[names(new.map)]
  names(ph) <- names(new.map)
  rf.vec <- mf_h(diff(new.map))
  mp <- as.data.frame(new.map)
  mp$mrk <- rownames(mp)
  P <- data.frame(P)
  P$mrk <- rownames(P) 
  P <- left_join(mp, P, by = "mrk")
  X <- P[,3:4]
  X[is.na(X)] <- 0
  h <- get_states_and_emission_single_parent(ploidy, ph, global.err, D, X[,uninfo.parent])
  res.temp <- .Call("calc_genprob_single_parent",
                    ploidy,
                    n.mrk,
                    n.ind,
                    h$states,
                    h$emission,
                    rf.vec,
                    as.numeric(rep(0, g * n.mrk * n.ind)),
                    verbose,
                    PACKAGE = "mappoly")
  if(verbose) cat("\n")
  dim(res.temp[[1]]) <- c(g,n.mrk,n.ind)
  dimnames(res.temp[[1]]) <- list(apply(combn(letters[1:ploidy],ploidy/2),2, paste, collapse = ""),
                                  rownames(D), colnames(D))
  structure(list(probs = res.temp[[1]], map = new.map), class = "mappoly.genoprob")
}
