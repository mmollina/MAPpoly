#' Frequency of genotypes for two-point recombination fraction estimation
#'
#' Returns the frequency of each genotype for two-point reduction
#' of dimensionality. The frequency is calculated for all pairwise
#' combinations and for all possible linkage phase configurations.
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#'
#' @param get.from.web If \code{TRUE}, access the counts for all
#'     linkage phase configurations in a remote server (default = FALSE)
#'
#' @param cache.prev an object of class \code{cache.info} containing
#'     pre-computed genotype frequencies, obtained with
#'     \code{\link[mappoly]{cache_counts_twopt}} (optional, default = NULL)
#'
#' @param n.clusters Number of parallel processes to spawn (default = 1)
#'
#' @param verbose If \code{TRUE} (default), print the linkage phase
#'     confgurations. If \code{get.from.web = TRUE}, nothing is
#'     printed, since all linkage phase configurations will be cached.
#'
#' @param joint.prob If \code{FALSE} (default), returns the frequency of
#'     genotypes for transition probabilities (conditional
#'     probabilities). If \code{TRUE} returns the frequency for joint
#'     probabilities.
#'
#' @return An object of class \code{cache.info} which contains one (conditional probabilities)
#'     or two (both conditional and joint probabilities) lists. Each list
#'     contains all pairs of dosages between parents for all markers
#'     in the sequence. The names in each list are of the form 'A-B-C-D', where: A
#'     represents the dosage in parent 1, marker k; B represents the dosage in parent
#'     1, marker k+1; C represents the dosage in parent 2, marker k;
#'     and D represents the dosage in parent 2, marker k+1.  For each
#'     list, the frequencies were computed for all possible linkage
#'     phase configurations. The frequencies for each linkage phase
#'     configuration are distributed in matrices whose names
#'     represents the number of homologous chromosomes that share
#'     alleles. The rows on these matrices represents the dosages in markers k
#'     and k+1 for an individual in the offspring. See Table 3 of
#'     S3 Appendix in Mollinari and Garcia (2019) for an example.
#'
#' @examples
#'   \dontrun{
#'     data(hexafake)
#'     all.mrk<-make_seq_mappoly(hexafake, 'all')
#'     ## local computation
#'     counts<-cache_counts_twopt(all.mrk, n.clusters = 8)
#'     ## download from web (specially important for high ploidy levels)
#'     counts.web<-cache_counts_twopt(all.mrk, get.from.web = TRUE)
#'     }
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
#' @export
#' @import parallel Rcpp RCurl
#' @importFrom stats na.omit
#' @importFrom utils combn read.table
#'
cache_counts_twopt <- function(input.seq, get.from.web = FALSE, cache.prev = NULL, n.clusters = 1, verbose = TRUE, joint.prob = FALSE) {
    ## checking for correct object
    start <- proc.time()
    input_classes <- c("mappoly.sequence")
    if (!inherits(input.seq, input_classes)) {
        stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
    }
    cache.prev = NULL
    if(input.seq$m==2)
    {
      get.from.web <- FALSE
      cat("INFO: Computing genotype frequencies ...\n")
    }
    if (get.from.web)
        return(get_cache_two_pts_from_web(input.seq$m))
    temp.count <- NULL
    if (joint.prob) {
        temp.count <- cache_counts_twopt(input.seq, get.from.web = FALSE, cache.prev = cache.prev, n.clusters = n.clusters, verbose = verbose, joint.prob = FALSE)$cond
    }
    if (input.seq$m >= 8)
        message("\nploidy level ", input.seq$m, ": this operation could take a very long time.
                 \ntry to use function 'get_cache_two_pts.from.web' instead.\n")
    dose.names <- sort(unique(apply(rbind(combn(input.seq$seq.dose.p, 2), combn(input.seq$seq.dose.q, 2)), 2, paste, collapse = "-")))
    if (!is.null(cache.prev)) {
        if (!inherits(cache.prev, "cache.info")) {
            stop(deparse(substitute(cache.prev)), " is not an object of class 'cache.info'")
        }
        ## Number of distinct genotypic combinations for differennt ploidy levels
        x <- c(3, 6, 10, 15, 21, 28, 36)
        names(x) <- c("2", "4", "6", "8", "10", "12", "14")
        pl.class <- choose(1 + input.seq$m/2, 2) + 1 + input.seq$m/2
        first.na <- sapply(unlist(cache.prev, recursive = FALSE), function(x) !all(is.na(x)))

        ## check if the ploidy levels match
        if (ncol(cache.prev[[min(which(first.na))]][[1]]) != pl.class) {
            warning(deparse(substitute(cache.prev)), " contains counts for ploidy level ", names(which(x == ncol(cache.prev[[1]][[1]]))), "\n  Obtaining new counts for ploidy ",
                input.seq$m)
            return(cache_counts_twopt(input.seq = input.seq, cache.prev = NULL, n.clusters = n.clusters, verbose = verbose, joint.prob = joint.prob))
        }
        remaining <- which(is.na(pmatch(dose.names, names(cache.prev))))
        if (length(remaining) == 0) {
            cat("\n Nothing to add to 'cache.prev'. Returning original 'cache prev'.\n")
            return(cache.prev)
        }
        dose.names <- dose.names[remaining]
    }
    aux.mat <- matrix(unlist(lapply(strsplit(dose.names, split = "-"), as.numeric)), ncol = 4, byrow = TRUE)
    dimnames(aux.mat) <- list(paste("Conf.", 1:nrow(aux.mat), sep = ""), c("P.k", "P.k+1", "Q.k", "Q.k+1"))
    if (verbose) {
        cat("\n   Caching the following dosage combination: \n")
        print(aux.mat)
    }
    if (n.clusters > 1) {
        if (verbose)
            cat("INFO: Using ", n.clusters, " CPU's for calculation.\n")
        cl <- makeCluster(n.clusters)
        on.exit(stopCluster(cl))
        y <- parApply(cl, aux.mat, 1, get_counts_all_phases, m = input.seq$m, joint.prob = joint.prob)
        end <- proc.time()
    } else {
        if (verbose)
            cat("INFO: Going singlemode. Using one Core/CPU/PC for calculation.\n")
        y <- apply(aux.mat, 1, get_counts_all_phases, input.seq$m, joint.prob = joint.prob)
        end <- proc.time()
    }
    if (verbose) {
        cat("INFO: Done with", nrow(aux.mat), "phase configurations\n")
        cat("INFO: Calculation took:", round((end - start)[3], digits = 3), "seconds\n")
    }
    names(y) <- dose.names
    w <- c(cache.prev, y)
    joint <- NULL
    cond <- temp.count
    if (joint.prob)
        joint <- w else cond <- w
    all.dose.names <- apply(expand.grid(0:input.seq$m, 0:input.seq$m, 0:input.seq$m, 0:input.seq$m), 1, paste, collapse = "-")
    z1 <- z2 <- vector("list", length(all.dose.names))
    names(z1) <- names(z2) <- all.dose.names
    z1[names(cond)] <- cond
    if (joint.prob)
        z2[names(joint)] <- joint
    structure(list(cond = z1, joint = z2), class = "cache.info")
}

#' @export
print.cache.info <- function(x, ...) {
  x$cond[sapply(x$cond, is.null)]<-NA
  NonNAindex <- which(!is.na(x$cond))
  firstNonNA <- min(NonNAindex)
  y <- c(3, 6, 10, 15, 21, 28, 36)
  names(y) <- c("2", "4", "6", "8", "10", "12", "14")
  cat("  This is an object of class 'cache.info'")
  cat("\n  -----------------------------------------------------")
  ## printing summary
  cat("\n  Ploidy level:                               ", names(which(ncol(x$cond[[firstNonNA]][[1]])==y)), "\n")
  cat("  No. marker combinations:                    ", length(x$cond), "\n")
  cat("  -----------------------------------------------------\n")
}


