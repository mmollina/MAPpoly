#' Compute genotype conditional probabilities
#'
#' Conditional genotype probabilities are calculated for each marker
#' position and each individual given a map. In this version,
#' the probabilities are not calculated between markers.
#'
#' @param input.map An object of class \code{mappoly.map}
#'
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the phase configuration associated with the
#'                     maximum likelihood
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @param ... currently ignored
#'
#' @return An object of class 'mappoly.genoprob' which has two elements: a tridimensional
#' array containing the probabilities of all possible genotypes for each individual
#' in each marker position; and the marker sequence with it's recombination frequencies
#'
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
#'    probs<-calc_genoprob(input.map = subset.map,
#'                         verbose = TRUE)
#'    image(t(probs$probs[,,1]), col = rev(heat.colors(100)))
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
#' @export calc_genoprob
#'
calc_genoprob<-function(input.map,  phase.config = "best", verbose = TRUE)
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
  m<-input.map$info$m
  n.ind<-get(input.map$info$data.name, pos=1)$n.ind
  n.mrk<-input.map$info$n.mrk
  D<-get(input.map$info$data.name, pos=1)$geno.dose[input.map$maps[[1]]$seq.num,]
  mrknames<-rownames(D)
  indnames<-colnames(D)
  dp <- get(input.map$info$data.name)$dosage.p[input.map$maps[[i.lpc]]$seq.num]
  dq <- get(input.map$info$data.name)$dosage.q[input.map$maps[[i.lpc]]$seq.num]
  for (j in 1:nrow(D))
    D[j, D[j, ] == input.map$info$m + 1] <- dp[j] + dq[j] + 1 + as.numeric(dp[j]==0 || dq[j]==0)
  res.temp<-.Call("calc_genoprob",
                  m,
                  t(D),
                  lapply(input.map$maps[[i.lpc]]$seq.ph$P, function(x) x-1),
                  lapply(input.map$maps[[i.lpc]]$seq.ph$Q, function(x) x-1),
                  input.map$maps[[i.lpc]]$seq.rf,
                  as.numeric(rep(0, choose(m, m/2)^2 * n.mrk * n.ind)),
                  verbose=verbose,
                  PACKAGE = "mappoly")
  if(verbose) cat("\n")
  dim(res.temp[[1]])<-c(choose(m,m/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]])<-list(kronecker(apply(combn(letters[1:m],m/2),2, paste, collapse=""),
                                          apply(combn(letters[(m+1):(2*m)],m/2),2, paste, collapse=""), paste, sep=":"),
                                  mrknames, indnames)
  structure(list(probs = res.temp[[1]], map = create_map(input.map)), class="mappoly.genoprob")
}

#' @export
print.mappoly.genoprob <- function(x, ...) {
  cat("  This is an object of class 'mappoly.genoprob'")
  cat("\n  -----------------------------------------------------")
  ## printing summary
cat("\n  No. genotypic classes:                    ", dim(x$probs)[1], "\n")
  cat("  No. positions:                            ", dim(x$probs)[2], "\n")
  cat("  No. individuals:                          ", dim(x$probs)[3], "\n")
  cat("  -----------------------------------------------------\n")
}

#' Create a marker with pseudomarkers
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
create_map<-function(input.map, step = Inf,
                     phase.config = "best")
{
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  } else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  mrknames <- get(input.map$info$data.name, pos=1)$mrk.names[input.map$maps[[i.lpc]]$seq.num]
  if(length(unique(input.map$maps[[1]]$seq.num)) == 1){
    a<-rep(0,input.map$info$n.mrk)
    names(a)<-mrknames
    return(a)
  }
  map <- c(0, cumsum(imf_h(input.map$maps[[i.lpc]]$seq.rf)))
  names(map)<-mrknames
  if(is.null(step))
    return(map)
  minloc <- min(map)
  map <- map-minloc
  a <- seq(floor(min(map)), max(map), by = step)
  a <- a[is.na(match(a,map))]
  names(a) <- paste("loc",a,sep = "_")
  return(sort(c(a,map))+minloc)
}

