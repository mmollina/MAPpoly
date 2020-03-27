#' Pairwise two-point analysis
#'
#' Performs the two-point pairwise analysis, accounting for the recombination fractions
#' between all pairwise markers given one or more linkage phase configurations.
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#'
#' @param count.cache an object of class \code{cache.info} containing
#'     pre-computed genotype frequencies, obtained with
#'     \code{\link[mappoly]{cache_counts_twopt}}
#'
#' @param n.clusters Number of parallel processes (cores) to spawn (default = 1)
#'
#' @param tol the desired accuracy. See \code{optimize()} for details
#'
#' @param mrk.pairs a matrix of dimensions 2*N, containing N
#'    pairs of markers to be analyzed. If \code{NULL} (default), all pairs are
#'    considered
#'
#' @param n.batches The number of batches of marker pairs that should be analyzed 
#'    in parallel.
#'
#' @param verbose If \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @param memory.warning if \code{TRUE}, prints a memory warning if the 
#' number of markers is greater than 2000. Mostly for internal use 
#' 
#' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
#' 
#' @return An object of class \code{poly.est.two.pts.pairwise} which
#'     is a list  containing the following components:
#'     \item{data.name}{name of the object of class
#'         \code{mappoly.data} with the raw data}
#'     \item{n.mrk}{number of markers in the sequence}
#'     \item{seq.num}{a \code{vector} containing the (ordered) indices
#'         of markers inthe sequence, according to the input file}
#'     \item{pairwise}{a list of size
#'     \code{choose(length(input.seq$seq.num), 2)}, each of them containing: a matrix
#'     with the linkage phase configuration i.e. the number of
#'     homologous that share alelels; the LOD Score relative to each
#'     one of these configurations; the recombination fraction
#'     estimated for each configuration; and their associated LOD
#'     score}
#'     \code{chisq.pval.thres}{threshold used to perform the segregation tests}
#'     \code{chisq.pval}{p-values associated with the performed segregation tests}
#'
#' @examples
#'   \dontrun{
#'   ## Tetraploid example: 
#'   all.mrk <- make_seq_mappoly(tetra.solcap, 'all')
#'   red.mrk <- elim_redundant(all.mrk)
#'   unique.mrks <- make_seq_mappoly(red.mrk)
#'   counts <- cache_counts_twopt(unique.mrks, cached = TRUE)
#'   # will take ~ 13 min
#'   all.pairs <- est_pairwise_rf(input.seq = unique.mrks,
#'                                count.cache = counts,
#'                                n.clusters = 7, 
#'                                verbose=TRUE)
#'    all.pairs
#'    plot(all.pairs, 90, 91)
#'    mat <- rf_list_to_matrix(all.pairs)
#'    plot(mat)
#'    
#'    ## Hexaploid example
#'    fl = "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/BT/sweetpotato_chr1.vcf.gz"
#'    tempfl <- tempfile(pattern = 'chr1_', fileext = '.vcf.gz')
#'    download.file(fl, destfile = tempfl)
#'    dat.dose.vcf = read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2")
#'    
#'    ## Filtering dataset by marker
#'    dat.filt.mrk <- filter_missing(input.data = dat.dose.vcf, 
#'                                   type = "marker", 
#'                                   filter.thres = 0.10, 
#'                                   inter = FALSE)
#'    ## Filtering dataset by individual
#'    dat.filt.ind <- filter_missing(input.data = dat.filt.mrk, 
#'                                   type = "individual", 
#'                                   filter.thres = 0.10, 
#'                                   inter = FALSE)
#'   ## Segregation test
#'   pval.bonf <- 0.05/dat.filt.ind$n.mrk
#'   mrks.chi.filt <- filter_segregation(dat.filt.ind,
#'                                       chisq.pval.thres =  pval.bonf, 
#'                                       inter = FALSE)
#'   seq.ch1<-make_seq_mappoly(mrks.chi.filt)
#'   plot(seq.ch1)
#'   counts <- cache_counts_twopt(seq.ch1, cached = TRUE)
#'   ## will take ~  19 min / peak of memory usage ~ 10GB
#'   all.pairs.1 <- est_pairwise_rf(input.seq = seq.ch1,
#'                                  count.cache = counts,
#'                                  n.clusters = 7, 
#'                                  verbose=TRUE)
#'   ## same thing, but it will take ~  21 min / peak of memory usage ~ 6GB
#'   all.pairs.2 <- est_pairwise_rf(input.seq = seq.ch1,
#'                                  count.cache = counts,
#'                                  n.clusters = 7, 
#'                                  n.batch = 10,
#'                                  verbose=TRUE)   
#'    plot(all.pairs, 90, 91)
#'    mat <- rf_list_to_matrix(all.pairs.1)
#'    plot(mat)
#'    }
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
#' @export est_pairwise_rf
#' @importFrom parallel makeCluster clusterEvalQ stopCluster parLapply
est_pairwise_rf <- function(input.seq, count.cache,
                            n.clusters = 1,
                            tol = .Machine$double.eps^0.25,
                            mrk.pairs = NULL,
                            n.batches = 1,
                            verbose = TRUE,
                            memory.warning = FALSE, 
                            parallelization.type = c("PSOCK", "FORK"))
{
  ## checking for correct object
  if (!class(input.seq) == "mappoly.sequence")
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  parallelization.type <- match.arg(parallelization.type)
  dpl <- duplicated(input.seq$seq.num)
  ## checking for duplicated markers
  if (any(dpl))
    stop("There are duplicated markers in the sequence:\n Check markers: ", unique(input.seq$seq.num[dpl]), " at position(s) ", which(dpl))
  # Memory warning
  ANSWER = "flag"

  if(input.seq$m < 6){
    if (length(input.seq$seq.num) > 10000 && interactive() && n.batches == 1 && !memory.warning){
      while (substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !=""){
        cat("  Ploidy level:", input.seq$m, "\n~~~~~~~~~~\n")
        message("
  The sequence contains more than 10000 markers. 
  This requires high-performance computing resources.
  Do you want to proceed? (Y/n): ")
        ANSWER <- readline("")
        if (substr(ANSWER, 1, 1) == "n" || substr(ANSWER, 1, 1) == "no") 
          stop("  You decided to stop 'est_pairwise_rf'.")
      }
    } 
  } else {
    if (length(input.seq$seq.num) > 3000 && interactive() && n.batches == 1 && !memory.warning){
      while (substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !=""){
        cat("  Ploidy level:", input.seq$m, "\n~~~~~~~~~~\n")
        message("
  The sequence contains more than 3000 markers. 
  This requires high-performance computing resources.
  Do you want to proceed? (Y/n): ")
        ANSWER <- readline("")
        if (substr(ANSWER, 1, 1) == "n" || substr(ANSWER, 1, 1) == "no") 
          stop("  You decided to stop 'est_pairwise_rf'.")
      }
    } 
  }
  geno <- get(input.seq$data.name, pos=1)$geno.dose
  ## all possible pairs
  if (is.null(mrk.pairs)) {
    mrk.pairs <- combn(sort(input.seq$seq.num), 2) - 1
  }
  batch.size <- NULL
  if(n.batches > 1)
    batch.size <- ceiling(ncol(mrk.pairs)/n.batches)
  if (is.null(batch.size)) {
    ## splitting pairs in chunks
    if (length(input.seq$seq.num) < 10)
      n.clusters <- 1
    id <- ceiling(seq(1, (ncol(mrk.pairs) + 1), length.out = n.clusters + 1))
    input.list <- vector("list", n.clusters)
    for (i in 1:n.clusters) input.list[[i]] <- mrk.pairs[, id[i]:(id[i + 1] - 1)]
    ## parallel version
    if (n.clusters > 1) {
      start <- proc.time()
      if (verbose)
        cat("INFO: Using ", n.clusters, " CPUs for calculation.\n")
      cl = parallel::makeCluster(n.clusters, type = parallelization.type)
      parallel::clusterEvalQ(cl, require(mappoly))
      on.exit(parallel::stopCluster(cl))
      res <- parallel::parLapply(cl,
                                 input.list,
                                 paralell_pairwise,
                                 input.seq = input.seq,
                                 geno = geno,
                                 dP = get(input.seq$data.name)$dosage.p,
                                 dQ = get(input.seq$data.name)$dosage.q,
                                 count.cache = count.cache,
                                 tol = tol)
      end <- proc.time()
      if (verbose) {
        cat("INFO: Done with",
            ncol(mrk.pairs),
            " pairs of markers \n")
        cat("INFO: Calculation took:",
            round((end - start)[3],
                  digits = 3),
            "seconds\n")
      }
    } else {
      if (verbose) {
        cat("INFO: Going singlemode. Using one CPU for calculation.\n")
        if (length(input.seq$seq.num) < 10)
          cat("Also, number of markers is too small to perform parallel computation.\n")
      }
      res <- lapply(input.list,
                    paralell_pairwise,
                    input.seq = input.seq,
                    geno = geno,
                    dP = get(input.seq$data.name)$dosage.p,
                    dQ = get(input.seq$data.name)$dosage.q,
                    count.cache = count.cache,
                    tol = tol)
    }
    res <- unlist(res,
                  recursive = FALSE)
    names(res) <- apply(mrk.pairs + 1,
                        2,
                        paste,
                        collapse = "-")
    return(structure(list(data.name = input.seq$data.name,
                          n.mrk = length(input.seq$seq.num),
                          seq.num = input.seq$seq.num,
                          pairwise = res,
                          chisq.pval.thres = input.seq$chisq.pval.thres,
                          chisq.pval = input.seq$chisq.pval,
                          nas  = sapply(res, function(x) any(is.na(x)))),
                     class = "poly.est.two.pts.pairwise"))
  } else {
    if (verbose)
      cat("INFO: There are ",
          ncol(mrk.pairs),
          " recombination fractions to be estimated.\n")
    id.batch <- c(seq(batch.size + 1,
                      ncol(mrk.pairs),
                      batch.size),
                  + 1)
    
    id.batch <- unique(c(seq(batch.size, ncol(mrk.pairs), batch.size), ncol(mrk.pairs)))
    id.batch <- cbind(c(1, 1+id.batch[1:(length(id.batch)-1)]), id.batch)
    res <- vector("list", ncol(mrk.pairs))
    if (verbose)
      cat("INFO: Estimating the first batch of ",
          batch.size,
          " recombination fractions.\n")
    z <- system.time(res[id.batch[1,1]:id.batch[1,2]] <- est_pairwise_rf(input.seq = input.seq,
                                                                         count.cache = count.cache,
                                                                         n.clusters = n.clusters,
                                                                         tol = tol,
                                                                         parallelization.type = parallelization.type,
                                                                         mrk.pairs = mrk.pairs[,id.batch[1,1]:id.batch[1,2]],
                                                                         verbose = FALSE,
                                                                         memory.warning = TRUE)$pairwise)
    if (verbose) {
      cat("INFO:",
          batch.size,
          " recombination fractions estimated in",
          round(z[3]/60, digits = 3), "minutes.\n")
    }
    if (verbose) {
      cat("INFO: Estimated time to complete the job (",
          nrow(id.batch),
          " batches):",
          round(z[3] * (nrow(id.batch) - 1)/60,
                digits = 3),
          " minutes.\n")
    }
    for (i in 2:nrow(id.batch)) {
      if (verbose)
        cat("batch ", i, " of ", nrow(id.batch), "\n")
      res[id.batch[i,1]:id.batch[i,2]] <- est_pairwise_rf(input.seq = input.seq,
                                                          count.cache = count.cache,
                                                          n.clusters = n.clusters,
                                                          tol = tol,
                                                          mrk.pairs = mrk.pairs[,id.batch[i,1]:id.batch[i,2]],
                                                          verbose = FALSE, 
                                                          memory.warning = TRUE, 
                                                          parallelization.type = parallelization.type)$pairwise
    }
  }
  return(structure(list(data.name = input.seq$data.name,
                        n.mrk = length(input.seq$seq.num),
                        seq.num = input.seq$seq.num,
                        pairwise = res,
                        chisq.pval.thres = input.seq$chisq.pval.thres,
                        chisq.pval = input.seq$chisq.pval,
                        nas  = sapply(res, function(x) any(is.na(x)))),
                   class = "poly.est.two.pts.pairwise"))
}

#' Wrapper function to pairwise two-point estimation in C++
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
paralell_pairwise <- function(mrk.pairs,
                              input.seq,
                              geno,
                              dP,
                              dQ,
                              count.cache,
                              tol = .Machine$double.eps^0.25)
{
  res <- .Call("pairwise_rf_estimation",
               input.seq$m,
               as.matrix(mrk.pairs),
               as.matrix(geno),
               as.vector(dP),
               as.vector(dQ),
               count.cache$cond,
               tol = tol,
               PACKAGE = "mappoly")
  return(lapply(res, format_rf))
}

#' Format results from pairwise two-point estimation in C++
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export format_rf
format_rf <- function(res) {
  x <- res
  if (length(x) != 4) {
    LOD_ph <- min(x[2, ]) - x[2, ]
    rf <- x[1, order(LOD_ph, decreasing = TRUE)]
    LOD_rf <- abs(x[2, ] - x[3, ])[order(LOD_ph, decreasing = TRUE)]
    LOD_ph <- LOD_ph[order(LOD_ph, decreasing = TRUE)]
    return(cbind(LOD_ph, rf, LOD_rf))
  } else {
    nm <- paste(min(x[1], x[2]), min(x[3], x[4]), sep = "-")
    x <- cbind(0, rf = NA, LOD = NA)
    rownames(x) <- nm
    colnames(x) <- c("ph_LOD", "rf", "rf_LOD")
    return(x)
  }
}

#' @rdname est_pairwise_rf
#' @keywords internal
#' @export
print.poly.est.two.pts.pairwise <- function(x, ...) {
  cat("  This is an object of class 'poly.est.two.pts.pairwise'")
  cat("\n  -----------------------------------------------------")
  ## printing summary
  cat("\n  No. markers:                            ", x$n.mrk, "\n")
  cat("  No. estimated recombination fractions:  ", choose(x$n.mrk, 2) - sum(x$nas))
  cat(" (", round(100 - 100 * sum(x$nas)/length(x$nas), 1), "%)\n", sep = "")
  cat("  -----------------------------------------------------\n")
}

#' @rdname est_pairwise_rf
#' @keywords internal
#' @export
plot.poly.est.two.pts.pairwise <- function(x, first.mrk, second.mrk, ...) {
  i<-which(names(x$pairwise)%in%paste(sort(c(first.mrk, second.mrk)), collapse = "-"))
  if(length(i)==0)
    stop("The requested combination of markers is not included in the two-point object")
  data.name <- x$data.name
  x<-x$pairwise[[i]]
  variable<-value<-sh<-NULL
  x<-reshape::melt(x, id=rownames(x))
  colnames(x)<-c("sh", "variable", "value")
  rfs<-as.character(format(round(t((subset(x, variable=="rf", select=value))), digits=2), digits=2))
  x.temp<-subset(x, variable!="rf")
  mLOD<-max(x.temp$value)
  x.temp <- transform(x.temp, sh = factor(sh, levels = unique(x.temp$sh)))
  ds<-paste0(paste0(get(data.name, pos=1)$dosage.p[first.mrk] , "---", get(data.name, pos=1)$dosage.p[second.mrk]), "  x  ",
             paste0(get(data.name, pos=1)$dosage.q[first.mrk] , "---", get(data.name, pos=1)$dosage.q[second.mrk]))
  ggplot2::ggplot(x.temp, ggplot2::aes(sh,    value, label = sh    ,fill=variable))  +
    ggplot2::geom_bar(stat = "identity", position = "identity") +
    ggplot2::annotate("text", x = c(1:length(rfs)), y = rep(mLOD + mLOD/20,length(rfs)), label = rfs, size=4) +
    ggplot2::annotate("text", x = 1, y = mLOD + mLOD/7, label = "recombination fractions", size=4, hjust=0, vjust=0, fontface=3)+
    ggplot2::scale_x_discrete(name=paste("Homologous sharing alleles   \n dosage in parents ", ds)) +
    ggplot2::scale_y_continuous(name="LOD") +
    ggplot2::scale_fill_manual(values=c("#E69F00", "#56B4E9"),
                               name="",
                               breaks=c("LOD_ph", "LOD_rf"),
                               labels=c("LOD_phase", "LOD_rf"))
}
