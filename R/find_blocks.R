#'  Allocate markers into linkage blocks
#'
#' Function to allocate markers into linkage blocks.  This is an 
#' EXPERIMENTAL FUNCTION and should be used with caution.
#'
#' @param input.seq an object of class \code{mappoly.sequence}.
#'
#' @param clustering.type if \code{'rf'}, it uses UPGMA clusterization based on 
#' the recombination fraction matrix to assemble blocks. Linkage blocks are 
#' assembled by cutting the clusterization tree at \code{rf.limit}. 
#' If \code{'genome'}, it splits the marker sequence at neighbor markers morre than 
#' \code{'genome.block.threshold'} apart.
#'
#' @param rf.limit the maximum value to consider linked markers in
#'     case of \code{'clustering.type = rf'}
#'
#' @param genome.block.threshold the threshold to assume markers are in the same linkage block.
#' to be considered when allocating markers into blocks in case of \code{'clustering.type = genomee'}
#' 
#' @param rf.mat an object of class \code{mappoly.rf.matrix}.
#' 
#' @param ncpus Number of parallel processes to spawn 
#'
#' @param ph.thres the threshold used to sequentially phase markers. 
#' Used in \code{thres.twopt} and \code{thres.hmm}. See \code{\link[mappoly]{est_rf_hmm_sequential}} 
#' for details.
#'
#' @param phase.number.limit the maximum number of linkage phases of the sub-maps.
#'  The default is 10. See \code{\link[mappoly]{est_rf_hmm_sequential}} for details.
#'     
#' @param error the assumed global genotyping error rate. If \code{NULL} (default) it does 
#' not include an error in the block estimation.
#' 
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced.
#'     
#' @param tol tolerance for the C routine, i.e., the value used to
#'     evaluate convergence.
#'
#' @param tol.err tolerance for the C routine, i.e., the value used to
#'     evaluate convergence, including the global genotyping error in the model.
#'
#' @return a list containing 1: a list of blocks in form of \code{mappoly.map} objects;
#'    2: a vector containing markers that were not included into blocks.
#'
#' @examples
#'   \dontrun{
#'   ## Selecting 50 markers in chromosome 5
#'   s5 <- make_seq_mappoly(tetra.solcap, "seq5")
#'   s5 <- make_seq_mappoly(tetra.solcap, s5$seq.mrk.names[1:50])
#'   tpt5 <- est_pairwise_rf(s5)
#'   m5 <- rf_list_to_matrix(tpt5, 3, 3)
#'   fb.rf <- find_blocks(s5, rf.mat = m5, verbose = FALSE, ncpus = 2)
#'   bl.rf <- fb.rf$blocks
#'   plot_map_list(bl.rf)
#'   
#'   ## Merging resulting maps
#'   map.merge <- merge_maps(bl.rf, tpt5)
#'   plot(map.merge, mrk.names = T)
#'   
#'   ## Comparing linkage phases with pre assembled map
#'   id <- na.omit(match(map.merge$info$mrk.names, solcap.err.map[[5]]$info$mrk.names))
#'   map.orig <- get_submap(solcap.err.map[[5]], mrk.pos = id)
#'   p1.m<-map.merge$maps[[1]]$seq.ph$P
#'   p2.m<-map.merge$maps[[1]]$seq.ph$Q
#'   names(p1.m) <- names(p2.m) <- map.merge$info$mrk.names
#'   p1.o<-map.orig$maps[[1]]$seq.ph$P
#'   p2.o<-map.orig$maps[[1]]$seq.ph$Q
#'   names(p1.o) <- names(p2.o) <- map.orig$info$mrk.names
#'   n <- intersect(names(p1.m), names(p1.o))
#'   plot_compare_haplotypes(4, p1.o[n], p2.o[n], p1.m[n], p2.m[n])
#'   
#'   ### Using genome
#'   fb.geno <- find_blocks(s5, clustering.type = "genome", genome.block.threshold = 10^4)
#'   plot_map_list(fb.geno$blocks)
#'   splt <- lapply(fb.geno$blocks, split_mappoly, 1)
#'   plot_map_list(splt)
#' }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @importFrom graphics abline
#' @importFrom stats as.dist cutree hclust
#' @export find_blocks
find_blocks <- function(input.seq,
                        clustering.type = c("rf", "genome"),
                        rf.limit = 1e-04,
                        genome.block.threshold = 10000,
                        rf.mat = NULL,
                        ncpus = 1,
                        ph.thres = 3,
                        phase.number.limit = 10,
                        error = 0.05,
                        verbose = TRUE, 
                        tol = 1e-2,
                        tol.err = 1e-3)
{
  dat.name <- input.seq$data.name
  clustering.type <- match.arg(clustering.type)
  if (clustering.type == "rf") {
    if (all(is.null(rf.mat)))
      stop("Please provide 'rf.mat'")
    i <- 1
    blocks <- NULL
    M <- rf.mat$rec.mat
    M[is.na(M)] <- 0.5
    while (length(M) > 1) {
      if (verbose) {
        if (i%%100 == 0)
          cat("block No.", i, " ---> remaining mrks: ", nrow(M),"\n")
      }
      a <- hclust(as.dist(M), method = "complete")
      ct <- cutree(a, h = 1e-4)
      blocks[[i]] <- rownames(M)[ct == which.max(table(ct))]
      id <- setdiff(rownames(M), unlist(blocks[1:i]))
      M <- M[id, id]
      i <- i + 1
    }
    blocks <- blocks[sapply(blocks, length) != 1]
  } else if (clustering.type == "genome") {
    if (is.null(input.seq$genome.pos))
      stop("There is no genome information.")
    so <- get_genomic_order(input.seq)$ord
    uso <- unique(so[, 1])
    p <- vector("list", length(uso))
    names(p) <- uso
    for (i in uso) {
      o <- order(input.seq$genome.pos[so[, 1] == i])
      sq <- c(which(diff(input.seq$genome.pos[so[, 1] == i][o]) > genome.block.threshold),
              length(input.seq$genome.pos[so[, 1] == i]))
      arg <- vector("list", length(sq))
      sqi <- c(1, sq + 1)
      for (j in 1:length(arg))
        arg[[j]] <- rownames(so)[so[, 1] == i & !is.na(match(so[, 2], input.seq$genome.pos[so[, 1] == i][o][sqi[j]:sq[j]]))]
      p[[as.character(i)]] <- arg
    }
    pn <- unlist(p, recursive = FALSE)
    nb <- sapply(pn, length)
    if (verbose) {
      cat("INFO: ", length(p), "chromosome \n      ",
          sum(nb > 1), "marker blocks\n      ",
          sum(nb == 1), "SNPs \n")
    }
    blocks <- pn[nb > 1]
  }
  if (ncpus > 1) {
    if (verbose)
      cat("INFO: Using ", ncpus, " CPUs for calculation.\n")
    cl = parallel::makeCluster(ncpus)
    parallel::clusterExport(cl, "parallel_block")
    parallel::clusterExport(cl, dat.name)
    on.exit(parallel::stopCluster(cl))
    res <- parallel::parLapply(cl,
                               blocks,
                               parallel_block,
                               dat.name = dat.name, 
                               ph.thres = ph.thres, 
                               tol = tol,
                               phase.number.limit = phase.number.limit, 
                               error = error, 
                               tol.err = tol.err)
  } else {
    res <- vector("list", length(blocks))
    for(i in 1:length(res)){
      res[[i]] <- parallel_block(mrk.vec = blocks[[i]] ,
                                 dat.name = dat.name, 
                                 ph.thres = ph.thres, 
                                 tol = tol,
                                 phase.number.limit = phase.number.limit, 
                                 error = error, 
                                 tol.err = tol.err)
    }
  } 
  res <- res[!is.na(res)]
  one.snp <- setdiff(input.seq$seq.mrk.names, unlist(lapply(res, function(x) x$info$mrk.names)))
  return(structure(list(blocks = res,
                        snps = one.snp),
                   class = "mappoly.block"))
}

#' Auxiliary function to estimate a map in a block of markers using parallel 
#' processing 
#' 
#' @keywords internal
#' @export parallel_block
parallel_block <- function(mrk.vec, 
                           dat.name,
                           ph.thres = 3, 
                           tol = 1e-2,
                           phase.number.limit = 20,
                           error = 0.05,
                           tol.err = 1e-3, 
                           verbose = FALSE) {
  func <- function(){
    dat <- get(dat.name, pos = 1)
    st <- make_seq_mappoly(dat, mrk.vec, data.name = dat.name)
    tp <- est_pairwise_rf(st, verbose = verbose)
    mp <- est_rf_hmm_sequential(input.seq = st,
                                start.set = 4,
                                thres.twopt = ph.thres,
                                thres.hmm = ph.thres,
                                info.tail = TRUE,
                                twopt = tp,
                                phase.number.limit = phase.number.limit,
                                tol = tol,
                                tol.final = tol.err, 
                                verbose = verbose)
    mp <- filter_map_at_hmm_thres(mp, 0.01)
    if(mp$info$n.mrk > 3) {
      mp.up <- split_and_rephase(mp, tp, gap.threshold = 1, size.rem.cluster = 3, verbose = verbose)
    } else if(mp$info$n.mrk > 2) {
      mp.up <- split_and_rephase(mp, tp, gap.threshold = 1, size.rem.cluster = 2, verbose = verbose)      
    } else {
      mp.up <- split_and_rephase(mp, tp, gap.threshold = 1, size.rem.cluster = 1, verbose = verbose)      
    } 
    mp.up <- est_full_hmm_with_global_error(mp.up, error = error, tol = tol.err, verbose = FALSE)
    mp.up    
  }
  return(tryCatch(func(), error = function(e) NA))
}
