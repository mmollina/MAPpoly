#'  Allocate markers into blocks
#'
#' Function to allocate markers located in the same genomic region or
#' linkage disequilibrium blocks into bins.
#'
#' @param input.seq an object of class \code{mappoly.sequence}.
#'
#' @param search.type one of \code{'rf'} or \code{'seq'}. If
#'     \code{'rf'}, the search for bins is performed using the
#'     recombination fraction estimates at a certain level given by
#'     the argument \code{rf.limit}. If \code{'seq'}, the search for
#'     bins is performed using the sequence information contained in
#'     the input file. It uses as threshold to assume the markers are
#'     at the same bin the argument \code{seq.limit}
#'
#' @param rf.limit the maximum value to consider linked markers in
#'     case of \code{'search.type=rf'}
#'
#' @param seq.limit the distance limit to be considered when allocating
#'     markers into bins in case of \code{'search.type=seq'}
#'
#' @param reconstruct if \code{TRUE}, reconstructs the genetic map in
#'     each bin. If \code{FALSE}, assume all recombination fractions
#'     in the bin equal to zero. Only makes effect if
#'     \code{search.type = 'seq'}. When \code{search.type = 'rf'},
#'     this procedure is automatically performed.
#'
#' @param extend.tail trhe length of the tail of the chain that should
#'     be used to calculate the likelihood of the linakge phases
#'
#' @param ord.limit to be documented
#'
#' @param n.clusters Number of parallel processes to spawn
#'
#' @param ph.thres the threshold used to determine if the linkage
#'     phases compared via two-point analysis should be considered
#'
#' @param rf.mat an object of class \code{mappoly.rf.matrix}.
#'
#' @param count.cache an object of class \code{cache.info} containing
#'     pre-computed genotype frequencies, obtained with
#'     \code{\link[mappoly]{cache_counts_twopt}}.
#'
#' @param tol tolerance for the C routine, i.e., the value used to
#'     evaluate convergence in intermidiate procedures.
#'
#' @param tol.final tolerance for the C routine, i.e., the value used to
#'     evaluate convergence in the final maps that costitute marker blocks.
#'
#' @param error global error rate
#'
#' @param verbose if \code{FALSE} (default), simplified output is
#'     displayed.  if \code{TRUE}, detailed output is displayed.
#'
#' @param ask if \code{TRUE}, ask if the function should proceed
#'            with the phasing and recombination fraction estimation
#'            (if \code{reconstruct == TRUE}) after find marker blocks.
#'
#' @param x an object of class \code{mappoly.bins}
#'
#' @param ... curentlly ignored
#'
#' @param rf.thres recombination fraction threshold for plot. If \code{NULL},
#'        all recombination fractions are plotted.
#'
#' @return an object of class \code{mappoly.bins}
#'
#' @examples
#'   \dontrun{
#'     hexa_file<-system.file('extdata', 'hexa_fake', package = 'mappoly')
#'     hexa_dat<-read_geno(file_in = hexa_file)
#'     all_mrk<-make_seq_mappoly(hexa_dat, 'all')
#'     counts_all_mrk_from_web<-cache_counts_twopt(input.seq=all_mrk,
#'                                                 get.from.web=TRUE)
#'     all_pairs<-est_pairwise_rf(all_mrk, counts_all_mrk_from_web,
#'                                n.clusters=2)
#'
#'     mat<-rf_list_to_matrix(twopt.input=all_pairs,
#'                            thresh.LOD.ph=10,
#'                            thresh.LOD.rf=10,
#'                            thresh.rf = 0.0005)
#'    plot(mat)
#'
#'    blocks_rf<-find_marker_blocks(input.seq = all_mrk, search.type='rf',
#'                       rf.limit=0.0001, n.clusters=2, ph.thres=3,
#'                       rf.mat=mat, count.cache=counts_all_mrk_from_web,
#'                       tol=10e-3)
#'    blocks_rf
#'    plot(blocks_rf)
#'
#'    blocks_seq<-find_marker_blocks(input.seq = all_mrk, search.type='seq', ph.thres=5,
#'                        n.clusters=2, count.cache=counts_all_mrk_from_web)
#'    blocks_seq
#'    plot(blocks_seq)
#'
#'    ## Autotetraploid potato
#'
#'     data("potato_solcap")
#'     ch1<-make_seq_mappoly(potato_solcap, 'seq1')
#'     counts_all_mrk_from_web<-cache_counts_twopt(ch1,
#'                                                 get.from.web=TRUE)
#'     ch1_pairs<-est_pairwise_rf(ch1,
#'                                counts_all_mrk_from_web,
#'                                n.clusters=16)
#'
#'     ch1_mat<- rf_list_to_matrix(twopt.input = ch1_pairs,
#'                                  thresh.LOD.ph = 10,
#'                                  thresh.LOD.rf = 10,
#'                                  thresh.rf = 0.001)
#'
#'    ch1_blocks_rf<-find_marker_blocks(input.seq = ch1,
#'                                      search.type = 'rf',
#'                                      rf.limit = 0.0001,
#'                                      n.clusters = 16,
#'                                      ph.thres = 3,
#'                                      rf.mat = ch1_mat,
#'                                      count.cache = counts_all_mrk_from_web,
#'                                      tol = 10e-3,
#'                                      error = 0.2)
#'    ch1_blocks_rf
#'    plot(ch1_blocks_rf)
#'
#'    ch1_blocks_seq<-find_marker_blocks(input.seq = ch1,
#'                                       search.type = 'seq',
#'                                       seq.limit = 10000,
#'                                       reconstruct = TRUE,
#'                                       ph.thres = 5,
#'                                       n.clusters = 16,
#'                                       count.cache = counts_all_mrk_from_web,
#'                                       tol = 10e-3,
#'                                       error = 0.2)
#'    ch1_blocks_seq
#'    plot(ch1_blocks_seq)
#'    }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#'
#' @importFrom graphics abline
#' @importFrom stats as.dist cutree hclust
#' @export find_marker_blocks
#'
find_marker_blocks <- function(input.seq,
                               search.type = c("rf", "seq", "orig.ord"),
                               rf.limit = 1e-04,
                               seq.limit = 1000,
                               ord.limit = 10,
                               reconstruct = FALSE,
                               extend.tail = 10,
                               n.clusters = 1,
                               ph.thres = 3,
                               rf.mat = NULL,
                               tol = 1e-03,
                               tol.final = 1e-04,
                               error = NULL,
                               verbose = TRUE,
                               count.cache,
                               ask = FALSE)
{
  search.type <- match.arg(search.type)
  mrk.names <- input.seq$seq.mrk.names
  sub.map.size.limit = 10
  ###Constructing marker blocks - rf###################
  ## 
  ## based on recombination fraction
  if (search.type == "rf") {
    if (all(is.null(rf.mat)))
      stop("Please provide the filtered recobination fraction matrix.")
    pair <- rf.mat$id
    i <- 1
    bin <- NULL
    while (length(pair) >= 2) {
      if (verbose) {
        if (i%%100 == 0)
          cat("bin No.",
              i,
              "remaining mrks to process: ",
              length(unique(as.vector(pair))),
              "\n")
      }
      id <- c(pair[1, 1], pair[pair[, 1] == pair[1, 1], 2])
      M <- rf.mat$rec.mat[as.character(id), as.character(id)]
      M[is.na(M)] <- 0.5
      a <- hclust(as.dist(M), method = "complete")
      ct <- cutree(a, h = rf.limit)
      bin[[i]] <- id[ct == which.max(table(ct))]
      rmv <- which(pair[, 1] %in% bin[[i]])
      if (length(rmv) > 0) {
        pair <- pair[-rmv, ]
        pair <- matrix(pair, ncol = 2)
      }
      rmv <- which(pair[, 2] %in% bin[[i]])
      if (length(rmv) > 0) {
        pair <- pair[-rmv, ]
        pair <- matrix(pair, ncol = 2)
      }
      i <- i + 1
    }
    id <- which(is.na(match(unique(as.vector(rf.mat$id)), unlist(bin))))
    if (length(id) > 0) {
      id <- unique(as.vector(rf.mat$id))[id]
      M <- rf.mat$rec.mat[as.character(id), as.character(id)]
      M[is.na(M)] <- 0.5
      a <- hclust(as.dist(M), method = "complete")
      ct <- cutree(a, h = rf.limit)
      for (j in unique(ct)) {
        bin[[i]] <- id[ct == j]
        i <- i + 1
      }
    }
    id <- which(is.na(match(input.seq$seq.num, unlist(bin))))
    if (length(id) > 0) {
      for (j in id) {
        bin[[i]] <- j
        i <- i + 1
      }
    }
    one.snp <- which(sapply(bin, length) == 1)
    if (length(one.snp) == 0)
      one.snp <- NA
    multiple.snps <- which(sapply(bin, length) > 1)
  }
  ###Constructing marker blocks - sequence###################
  ## 
  ## based on sequence information
  else if (search.type == "seq")
  {
    if (is.null(input.seq$sequence))
      stop("There is no sequence information, please use another search method")
    so <- cbind(input.seq$sequence, input.seq$sequence.pos)
    so[is.na(so)]<-0
    dimnames(so) <- list(input.seq$seq.mrk.names, c("sequence", "sequencepos"))
    uso <- unique(so[, 1])
    p <- vector("list", length(uso))
    names(p) <- uso
    for (i in uso) {
      o <- order(input.seq$sequence.pos[so[, 1] == i])
      sq <- c(which(diff(input.seq$sequence.pos[so[, 1] == i][o]) > seq.limit),
              length(input.seq$sequence.pos[so[, 1] == i]))
      arg <- vector("list", length(sq))
      sqi <- c(1, sq + 1)
      for (j in 1:length(arg))
        arg[[j]] <- rownames(so)[so[, 1] == i &
                                   !is.na(match(so[, 2],
                                                input.seq$sequence.pos[so[, 1] == i][o][sqi[j]:sq[j]])
                                   )
                                 ]
      p[[as.character(i)]] <- arg
    }
    pn <- unlist(p, recursive = FALSE)
    nb <- sapply(pn, length)
    one.snp <- pn[nb == 1]
    if (length(one.snp) == 0)
      one.snp <- NA
    multiple.snps <- pn[nb > 1]
    if (verbose) {
      cat("INFO: ", length(p), "sequence found\n      ",
          sum(nb > 1), "marker blocks found\n      ",
          sum(nb == 1), "SNPs found\n----------------\n")
      cat("Distribution of the number of snps per block:\n")
      stem(sapply(multiple.snps, length))
    }
  }
  else if (search.type == "orig.ord")
  {
    if (is.null(input.seq$sequence))
      stop("There is no sequence information, please use another search method")
    so <- cbind(0, 1:length(input.seq$seq.num))
    dimnames(so) <- list(input.seq$seq.mrk.names, c("sequence", "sequencepos"))
    uso <- unique(so[, 1])
    p <- vector("list", length(uso))
    names(p) <- uso
    sq <- unique(c(seq(1, length(input.seq$seq.num), ord.limit), length(input.seq$seq.num)))
    arg <- vector("list", (length(sq)-1))
    for (j in 1:length(arg))
    arg[[j]] <- rownames(so)[sq[j]:(sq[j+1]-1)]
    p[[1]] <- arg
    pn <- unlist(p, recursive = FALSE)
    nb <- sapply(pn, length)
    one.snp <- pn[nb == 1]
    if (length(one.snp) == 0)
      one.snp <- NA
    multiple.snps <- pn[nb > 1]
    if (verbose) {
      cat("INFO: ", length(p), "sequence found\n      ",
          sum(nb > 1), "marker blocks found\n      ",
          sum(nb == 1), "SNPs found\n----------------\n")
    }
  }
  else stop("Invalid search type")
  if(interactive() && ask)
  {
    ANSWER <- readline("Do you want to proceed with analysis? ")
    if (substr(ANSWER, 1, 1) == "n")
      return(NULL)
  }
  if (verbose){
    cat("INFO: ",
        length(multiple.snps),
        "bins found\nGathering info to procede\n")
    cat("with recombinnation fraction computation\n")
  }
  ###Gathering two-point info###################
  index.seq <- vector("list", length(multiple.snps))
  for (i in 1:length(multiple.snps)) {
    if (verbose) {
      cat(".")
      if (i%%50 == 0)
        cat(" --> ", i, "/", length(multiple.snps), "\n", sep = "")
    }
    if (search.type == "rf")
      s <-  make_seq_mappoly(get(input.seq$data.name, pos = 1),
                             bin[[multiple.snps[i]]],
                             data.name = input.seq$data.name)
    else
      s <-  make_seq_mappoly(get(input.seq$data.name, pos = 1),
                             multiple.snps[[i]],
                             data.name = input.seq$data.name)
    twopt <- est_pairwise_rf(input.seq = s,
                             count.cache,
                             n.clusters = 1,
                             verbose = FALSE)
    index.seq[[i]] <- list(i = i,
                           seq.num = s,
                           twopt = twopt)
  }
  ### Estimating marker blocks #####
  ## parallel version
  if (n.clusters > 1) {
    start <- proc.time()
    if (verbose)
      cat("\n\nINFO: Using ", n.clusters, " Cores for calculation.\n")
    cl <- parallel::makeCluster(n.clusters)
    parallel::clusterEvalQ(cl, require(mappoly))
    parallel::clusterExport(cl,
                            varlist = c(input.seq$data.name,
                                        "ph.thres",
                                        "extend.tail",
                                        "tol",
                                        "tol.final",
                                        "sub.map.size.limit",
                                        "error"),
                            envir = environment())
    on.exit(parallel::stopCluster(cl))
    ## Call C++ routine
    bins <- parallel::clusterApplyLB(cl,
                                     index.seq,
                                     parallel_hmm_bins,
                                     thres.twopt = ph.thres,
                                     thres.hmm = 3,
                                     extend.tail = extend.tail,
                                     tol = tol,
                                     tol.final = tol.final,
                                     sub.map.size.limit = sub.map.size.limit,
                                     error = error)
    end <- proc.time()
    if (verbose) {
      cat("INFO: Done with", length(bins), " bins \n")
      cat("INFO: Calculation took:",
          round((end - start)[3],
                digits = 3), "seconds\n")
    }
  } else {
    bins <- lapply(index.seq,
                   parallel_hmm_bins,
                   thres.twopt = ph.thres,
                   thres.hmm = 3,
                   extend.tail = extend.tail,
                   tol = tol,
                   tol.final = tol.final,
                   sub.map.size.limit = sub.map.size.limit,
                   error = error)
  }
  names(bins) <- 1:length(bins)
  if (search.type == "rf")
    one.snp <- unlist(bin[one.snp])
  return(structure(list(bins = bins,
                        elim.bins = NULL,
                        snps = one.snp,
                        rf.thres = NULL,
                        seq.bins = names(bins),
                        maps = NULL,
                        error = error),
                   class = "mappoly.bins"))
}


#' Wrap function to estimate a map using parallel processing (for find_marker_blocks)
#'
#' @param index.seq sequence used as argument in function make_seq_temp
#' @param dat an object of class \code{mappoly.data}
#' @param thres.twopt the threshold used to determine if the linakge
#'     phases compared via two-point analysis should be considered
#' @param thres.hmm the threshold used to determine if the linakge
#'     phases compared via hmm analysis should be considered
#' @param extend.tail trhe length of the tail of the chain that should
#'     be used to calculate the likelihood of the linakge phases
#' @param dat.name name of the object of class \code{mappoly.data}
#' @param tol the desired accuracy during the sequential phase.
#' @param tol.final the desired accuracy for the final map.
#' @param block.estimate logical. If TRUE returns a map forcing all
#'     recombination fractions equal to 0 (1e-5)
#' @param error global error rate
#' @param verbose print the progress
#' @return a list of maps
#' @keywords internal
#' @export parallel_hmm_bins
parallel_hmm_bins <- function(index.seq,
                              thres.twopt = 10,
                              thres.hmm = 3,
                              extend.tail = 10,
                              tol = 10e-3,
                              tol.final = 10e-4,
                              sub.map.size.limit = 7,
                              phase.number.limit = 4,
                              error = NULL,
                              verbose = FALSE)
{
  w<-tryCatch(
    bla<-est_rf_hmm_sequential(input.seq = index.seq$s,
                               start.set = 4,
                               thres.twopt = thres.twopt,
                               thres.hmm = thres.hmm,
                               extend.tail = extend.tail,
                               twopt = index.seq$twopt,
                               verbose = verbose,
                               tol = tol,
                               tol.final = tol,
                               sub.map.size.diff.limit = sub.map.size.limit,
                               phase.number.limit = phase.number.limit),
    error = function(e) NA)
  if(!is.na(w) && !is.null(error))
    w<-est_full_hmm_with_global_error(input.map = w,
                                      error = error,
                                      tol = tol.final,
                                      verbose = FALSE)
  return(w)
}

#' Filter bins given a recombination fraction threshold
#'
#' Given a object of class \code{mappoly.bins}, returns a new object
#' of class \code{mappoly.bins} eliminating the bins whose
#' recombination fractions exceed a given threshold
#'
#' @param input.obj an object of class \code{mappoly.bins}
#'
#' @param rf.thres recombination fraction threshold used to eliminate
#'     bins
#'
#' @export filter_marker_blocks
filter_marker_blocks <- function(input.obj, rf.thres) {
  if (!(class(input.obj) == "mappoly.bins"))
    stop(deparse(substitute(input.obj)), " is not an object of class 'mappoly.bins'")
  single.bins <- lapply(input.obj$bins, function(x) x$maps[[which.max(sapply(x$maps, function(x) x$loglike))]])
  elim <- which(sapply(single.bins, function(x) any(x$seq.rf > rf.thres)))
  maps <- input.obj$maps
  if (length(elim) > 0 && !is.null(maps)) {
    warning("Number of bins changed, please reestimate the map.")
    maps <- NULL
  }
  structure(list(bins = input.obj$bins[-elim], elim.bins = input.obj$bins[elim], snps = input.obj$snps, rf.thres = rf.thres, seq.bins = input.obj$seq.bins,
                 maps = maps), class = "mappoly.bins")
}

#' @rdname find_marker_blocks
#' @export
print.mappoly.bins <- function(x, ...) {
  if (!any(class(x) == "mappoly.bins"))
    stop(deparse(substitute(x)), " is not an object of class 'mappoly.bins'")

  cat("  This is an object of class 'mappoly.bins'\n  ------------------------------------------\n")
  m <- x$bins[[1]]$info$m
  n.bins <- length(x$bins)
  n.mrk <- sapply(x$bins, function(x) x$info$n.mrk)
  tot.n.mrk <- sum(n.mrk)
  hap.len <- table(n.mrk)
  single.bins <- lapply(x$bins, function(x) x$maps[[which.max(sapply(x$maps, function(x) x$loglike))]])
  P.hap <- lapply(single.bins, function(x) ph_list_to_matrix(x$seq.ph$P, m))
  Q.hap <- lapply(single.bins, function(x) ph_list_to_matrix(x$seq.ph$Q, m))
  cod <- numeric(n.bins)
  for (i in 1:n.bins) cod[i] <- length(unique(c(apply(P.hap[[i]], 2, paste, collapse = ""), apply(Q.hap[[i]], 2, paste, collapse = ""))))
  elim.mrk <- 0
  if (!is.null(x$elim.bins))
    elim.mrk <- sum(sapply(x$elim.bins, function(x) x$info$n.mrk))
  bins.rf <- sapply(single.bins, function(x) x$seq.rf)
  cat("\n  Total No. of markers:                ", tot.n.mrk + length(x$snps) + elim.mrk, "\n")
  cat("  No. of bins:                         ", length(cod), "\n")
  cat("  No. markers into used bins           ", tot.n.mrk, "\n")
  cat("  No. single mrks:                     ", length(x$snps), "\n")
  cat("  No. eliminated bins:                 ", length(x$elim.bins), "\n")
  cat("  No. of markers into eliminated bins: ", elim.mrk, "\n")
  if (!is.null(x$rf.thres))
    cat("  rf threshold used to eliminate bins: ", x$rf.thres, "\n")
  cat("  Rec. Frac. range:                     ")
  cat(round(range(unlist(bins.rf)), 2), sep = "--")
  cat("\n")
  cat("\n  Number of markers per haplotypic alleles: ")
  print(table(cod))
}

#' @rdname find_marker_blocks
#' @export
plot.mappoly.bins <- function(x, rf.thres = NULL, ...) {
  m <- x$bins[[1]]$info$m
  n.bins <- length(x$bins)
  n.mrk <- sapply(x$bins, function(x) x$info$n.mrk)
  tot.n.mrk <- sum(n.mrk)
  hap.len <- table(n.mrk)
  single.bins <- lapply(x$bins, function(x) x$maps[[which.max(sapply(x$maps, function(x) x$loglike))]])
  P.hap <- lapply(single.bins, function(x) ph_list_to_matrix(x$seq.ph$P, m))
  Q.hap <- lapply(single.bins, function(x) ph_list_to_matrix(x$seq.ph$Q, m))
  cod <- numeric(n.bins)
  for (i in 1:n.bins) cod[i] <- length(unique(c(apply(P.hap[[i]], 2, paste, collapse = ""), apply(Q.hap[[i]], 2, paste, collapse = ""))))
  bins.rf <- sapply(single.bins, function(x) x$seq.rf)
  cod.old <- cod
  cod <- cod * max(unlist(bins.rf))/max(cod)
  op <- par(bty = "n")
  plot(x = -10:n.bins, y = rep(0, n.bins + 11), ylim = c(0, max(unlist(bins.rf))), type = "n", xlab = "bins", ylab = "recombintion fraction")
  if (is.null(rf.thres))
    rf.thres <- 0.5
  w <- lapply(bins.rf, function(x) x > rf.thres)
  for (i in 1:length(w)) {
    for (j in 1:length(w[[i]])) {
      if (w[[i]][j])
        w[[i]][j] <- as.character(i) else w[[i]][j] <- "*"
    }
  }
  ct <- 1
  d <- !duplicated(cod)
  for (i in rev(order(cod))) {
    text(x = rep(ct, length(bins.rf[[i]])), y = bins.rf[[i]], labels = as.character(w[[i]]))
    points(x = ct, y = cod[i], pch = 17, col = 2)
    if (d[i]) {
      text(x = -10, y = cod[i], labels = as.character(cod.old[i]), col = 2)
      abline(h = cod[i], lwd = 0.5, lty = 2)
    }
    ct <- ct + 1
  }
  par(op)
  invisible(as.numeric(unlist(w)[unlist(w) != "*"]))
}
