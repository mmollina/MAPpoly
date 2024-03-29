#' Compute conditional probabilities of the genotypes using probability distribution of dosages 
#'
#' Conditional genotype probabilities are calculated for each marker
#' position and each individual given a map. In this function,
#' the probabilities are not calculated between markers.
#'
#' @param input.map An object of class \code{mappoly.map}
#'
#' @param dat.prob an object of class \code{mappoly.data} containing the
#'                 probability distribution of the genotypes
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the phase configuration with the
#'                     maximum likelihood
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#'
#' @return An object of class 'mappoly.genoprob' which has two elements: a tridimensional
#' array containing the probabilities of all possible genotypes for each individual
#' in each marker position; and the marker sequence with it's recombination frequencies
#' 
#' @examples
#'  ## tetraploid example
#'  probs.t <- calc_genoprob_dist(input.map = solcap.prior.map[[1]],
#'                            dat.prob = tetra.solcap.geno.dist,
#'                            verbose = TRUE)
#'  probs.t
#'  ## displaying individual 1, 36 genotypic states 
#'  ## (rows) across linkage group 1 (columns)                          
#'  image(t(probs.t$probs[,,1]))
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
#' @importFrom dplyr select
#' 
#' @export calc_genoprob_dist
calc_genoprob_dist <- function(input.map, dat.prob = NULL, phase.config = "best", verbose = TRUE)
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
  if(is.null(dat.prob)){
    if(nrow(get(input.map$info$data.name, pos = 1)$geno) == get(input.map$info$data.name, pos = 1)$n.mrk){
      stop("
            The dataset associated to 'input.map'
            contains no genotypic probability distribution.
            Provide dataset using the argument 
           'dat.prob'.
           ")
    } else {
      dat.prob <- get(input.map$info$data.name, pos = 1)
    }
  }
  ind <- mrk <- NULL
  original.map.mrk <- input.map$info$mrk.names
  dat.prob.pos <- match(original.map.mrk, dat.prob$mrk.names)
  which.is.na <- which(is.na(dat.prob.pos))
  if(length(which.is.na) > 0)
    stop("Markers", original.map.mrk[which.is.na], "are not present in the 'dat.prob' object")
  temp.map <- input.map
  temp.map$info$seq.num <- temp.map$maps[[i.lpc]]$seq.num <- dat.prob.pos
  names(temp.map$maps[[i.lpc]]$seq.ph$P) <- names(temp.map$maps[[i.lpc]]$seq.ph$Q) <- dat.prob.pos
  if(!all(sort(get(temp.map$info$data.name, pos = 1)$ind.names) %in% sort(get(input.map$info$data.name, pos = 1)$ind.names)))
    stop("Individuals are different between dataset")
  g <- t(dat.prob$geno %>% filter(mrk%in%original.map.mrk) %>% select(!(mrk:ind)))
  ploidy = as.numeric(temp.map$info$ploidy)
  n.mrk = as.numeric(temp.map$info$n.mrk)
  n.ind = get(temp.map$info$data.name, pos = 1)$n.ind
  p1 = as.numeric(unlist(temp.map$maps[[1]]$seq.ph$P))
  dp1 = as.numeric(cumsum(c(0, sapply(temp.map$maps[[1]]$seq.ph$P, function(x) sum(length(x))))))
  p2 = as.numeric(unlist(temp.map$maps[[1]]$seq.ph$Q))
  dp2 = as.numeric(cumsum(c(0, sapply(temp.map$maps[[1]]$seq.ph$Q, function(x) sum(length(x))))))
  rf = temp.map$maps[[1]]$seq.rf
  ind.names <- dat.prob$ind.names
  res.temp  <- 
    .Call(
      "calc_genoprob_prior",
      as.numeric(ploidy),
      as.numeric(n.mrk),
      as.numeric(n.ind),
      as.numeric(p1),
      as.numeric(dp1),
      as.numeric(p2),
      as.numeric(dp2),
      as.double(g),
      as.double(rf),
      as.numeric(rep(0, choose(ploidy, ploidy/2)^2 * n.mrk * n.ind)),
      as.double(0),
      as.numeric(verbose),
      PACKAGE = "mappoly"
    )
  if (verbose) cat("\n")
  dim(res.temp[[1]]) <- c(choose(ploidy,ploidy/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]]) <- list(kronecker(apply(combn(letters[1:ploidy],ploidy/2),2, paste, collapse = ""),
                                          apply(combn(letters[(ploidy+1):(2*ploidy)],ploidy/2),2, paste, collapse = ""), paste, sep = ":"),
                                original.map.mrk, ind.names)
  structure(list(probs = res.temp[[1]], map = create_map(input.map)), class = "mappoly.genoprob")
}
