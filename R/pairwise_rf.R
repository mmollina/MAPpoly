#' Pairwise two-point analysis
#'
#' Performs the two-point pairwise analysis between all markers in a sequence. 
#' For each pair, the function estimates the recombination fraction for all 
#' possible linkage phase configurations and associated LOD Scores. 
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#'
#' @param count.cache an object of class \code{cache.info} containing
#'     pre-computed genotype frequencies, obtained with
#'     \code{\link[mappoly]{cache_counts_twopt}}. If \code{NULL} (default),
#'     genotype frequencies are internally loaded.
#'        
#' @param count.matrix similar to \code{count.cache}, but in matrix format. 
#' Mostly for internal use.
#'
#' @param ncpus Number of parallel processes (cores) to spawn (default = 1)
#'
#' @param mrk.pairs a matrix of dimensions 2*N, containing N
#'    pairs of markers to be analyzed. If \code{NULL} (default), all pairs are
#'    considered
#'
#' @param n.batches deprecated. Not available on MAPpoly 0.3.0 or higher
#'    
#' @param est.type  Indicates whether to use the discrete ("disc") or the probabilistic ("prob") dosage scoring 
#'                  when estimating the two-point recombination fractions. 
#'
#' @param verbose If \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'     
#' @param memory.warning if \code{TRUE}, prints a memory warning if the 
#' number of markers is greater than 10000 for ploidy levels up to 4, and 
#' 3000 for ploidy levels > 4. 
#' 
#' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
#'
#' @param tol the desired accuracy. See \code{optimize()} for details
#' 
#' @param ll will return log-likelihood instead of LOD scores. 
#' (for internal use) 
#' 
#' @return An object of class \code{mappoly.twopt} which
#'     is a list  containing the following components:
#'     \item{data.name}{name of the object of class
#'         \code{mappoly.data} with the raw data}
#'     \item{n.mrk}{number of markers in the sequence}
#'     \item{seq.num}{a \code{vector} containing the (ordered) indices
#'         of markers in the sequence, according to the input file}
#'     \item{pairwise}{a list of size
#'     \code{choose(length(input.seq$seq.num), 2)}, each of them containing a 
#'     matrix where the name of the rows have the form x-y, where x and y indicate 
#'     how many homologues share the same allelic variant in parents P and Q, 
#'     respectively (see Mollinari and Garcia, 2019 for notation). The first 
#'     column indicates the LOD Score in relation to the most likely linkage 
#'     phase configuration. The second column shows the estimated recombination 
#'     fraction for each configuration, and the third indicates the LOD Score 
#'     comparing the likelihood under no linkage (r = 0.5) with the estimated 
#'     recombination fraction (evidence of linkage).}
#'     \code{chisq.pval.thres}{threshold used to perform the segregation tests}
#'     \code{chisq.pval}{p-values associated with the performed segregation tests}
#'
#' @examples
#'   ## Tetraploid example (first 50 markers) 
#'   all.mrk <- make_seq_mappoly(tetra.solcap, 1:50)
#'   red.mrk <- elim_redundant(all.mrk)
#'   unique.mrks <- make_seq_mappoly(red.mrk)
#'   all.pairs <- est_pairwise_rf(input.seq = unique.mrks,
#'                                ncpus = 1, 
#'                                verbose = TRUE)
#'    all.pairs
#'    plot(all.pairs, 20, 21)
#'    mat <- rf_list_to_matrix(all.pairs)
#'    plot(mat)
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378}
#'     
#' @export est_pairwise_rf
#' @importFrom parallel makeCluster clusterEvalQ stopCluster parLapply
#' @importFrom Rcpp sourceCpp
#' @importFrom reshape2 melt acast
#' @importFrom dplyr filter arrange
est_pairwise_rf <- function(input.seq, count.cache = NULL, 
                            count.matrix = NULL, ncpus = 1L,
                            mrk.pairs = NULL, n.batches = 1L,
                            est.type = c("disc","prob"),
                            verbose = TRUE, memory.warning = TRUE, 
                            parallelization.type = c("PSOCK", "FORK"), 
                            tol = .Machine$double.eps^0.25, ll = FALSE)
{
  ## checking for correct object
  if (!inherits(input.seq, "mappoly.sequence"))
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  parallelization.type <- match.arg(parallelization.type)
  dpl <- duplicated(input.seq$seq.num)
  ## checking for duplicated markers
  if (any(dpl))
    stop("There are duplicated markers in the sequence:\n Check markers: ", unique(input.seq$seq.num[dpl]), " at position(s) ", which(dpl))
  if(is.null(count.cache))
    count.cache = cache_counts_twopt(input.seq, cached = TRUE)
  # Memory warning
  ANSWER = "flag"
  if(input.seq$ploidy < 6){
    if (length(input.seq$seq.num) > 10000 && interactive() && !memory.warning){
      while (substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != ""){
        cat("  Ploidy level:", input.seq$ploidy, "\n~~~~~~~~~~\n")
        message("
  The sequence contains more than 10000 markers. 
  This requires high-performance computing resources.
  Do you want to proceed? (Y/n): ")
        ANSWER <- readline("")
        if (substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N") 
          stop("  You decided to stop 'est_pairwise_rf'.")
      }
    } 
  } 
  else {
    if (length(input.seq$seq.num) > 3000 && interactive() && !memory.warning){
      while (substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != ""){
        cat("  Ploidy level:", input.seq$ploidy, "\n~~~~~~~~~~\n")
        message("
  The sequence contains more than 3000 markers. 
  This requires high-performance computing resources.
  Do you want to proceed? (Y/n): ")
        ANSWER <- readline("")
        if (substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N") 
          stop("  You decided to stop 'est_pairwise_rf'.")
      }
    } 
  }
  est.type = match.arg(est.type)
  ## Checking for genotype probability 
  if(!exists('geno', where = get(input.seq$data.name, pos = 1)) & est.type != "disc"){
    warning("There is no probabilistic dosage scoring in the dataset. Using est.type = 'disc'")
    est.type <- "disc"
  }
  ## get genotypes
  if(est.type  ==  "disc"){
    geno <- as.matrix(get(input.seq$data.name, pos = 1)$geno.dose)
  } 
  else {
    d1 <- get(input.seq$data.name, pos = 1)$geno 
    d2 <- reshape2::melt(d1, id.vars = c("mrk", "ind"))
    geno <- reshape2::acast(d2, mrk ~ variable ~ ind)
    geno <- geno[get(input.seq$data.name, pos = 1)$mrk.names,,get(input.seq$data.name, pos = 1)$ind.names]
  }
  ## all possible pairs
  if (is.null(mrk.pairs)) {
    mrk.pairs <- combn(sort(input.seq$seq.num), 2) - 1
  } else {
    mrk.pairs <- mrk.pairs - 1
  }
  if(n.batches > 1) 
    stop("n.batch is a deprecated argument. To avoid memory 
          overflow, use function 'est_pairwise_rf2' instead. 
          After grouping markers, you should use 
          'est_pairwise_rf'.")
  ## splitting pairs in chunks
  if (length(input.seq$seq.num) < 10)
    ncpus <- 1
  id <- ceiling(seq(1, (ncol(mrk.pairs) + 1), length.out = ncpus + 1))
  input.list <- vector("list", ncpus)
  for (i in 1:ncpus) input.list[[i]] <- mrk.pairs[, id[i]:(id[i + 1] - 1)]
  ## parallel version
  if (ncpus > 1) {
    start <- proc.time()
    if (verbose)
      cat("INFO: Using ", ncpus, " CPUs for calculation.\n")
    cl = parallel::makeCluster(ncpus, type = parallelization.type)
    if(est.type  ==  "disc")
      parallel::clusterExport(cl, "paralell_pairwise_discrete")
    if(est.type  ==  "prob")
      parallel::clusterExport(cl, "paralell_pairwise_probability")
    on.exit(parallel::stopCluster(cl))
    if(est.type  ==  "disc"){
      res <- parallel::parLapply(cl,
                                 input.list,
                                 paralell_pairwise_discrete,
                                 input.seq = input.seq,
                                 geno = geno,
                                 dP = get(input.seq$data.name)$dosage.p1,
                                 dQ = get(input.seq$data.name)$dosage.p2,
                                 count.cache = count.cache,
                                 tol = tol,
                                 ll = ll)
    } 
    else if(est.type  ==  "prob") {
      res <- parallel::parLapply(cl,
                                 input.list,
                                 paralell_pairwise_probability,
                                 input.seq = input.seq,
                                 geno = geno,
                                 dP = get(input.seq$data.name)$dosage.p1,
                                 dQ = get(input.seq$data.name)$dosage.p2,
                                 count.cache = count.cache,
                                 tol = tol,
                                 ll =ll)
    } 
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
  } 
  else {
    if (verbose) {
      cat("INFO: Going singlemode. Using one CPU for calculation.\n")
      if (length(input.seq$seq.num) < 10)
        cat("Also, number of markers is too small to perform parallel computation.\n")
    }
    if(est.type  ==  "disc"){
      res <- lapply(input.list,
                    paralell_pairwise_discrete,
                    input.seq = input.seq,
                    geno = geno,
                    dP = get(input.seq$data.name)$dosage.p1,
                    dQ = get(input.seq$data.name)$dosage.p2,
                    count.cache = count.cache,
                    tol = tol,
                    ll = ll)
    } 
    else if(est.type  ==  "prob") {
      res <- lapply(input.list,
                    paralell_pairwise_probability,
                    input.seq = input.seq,
                    geno = geno,
                    dP = get(input.seq$data.name)$dosage.p1,
                    dQ = get(input.seq$data.name)$dosage.p2,
                    count.cache = count.cache,
                    tol = tol,
                    ll = ll)
    } 
  }
  res <- unlist(res,
                recursive = FALSE)
  names(res) <- apply(mrk.pairs + 1,
                      2,
                      paste,
                      collapse = "-")
  nas <- sapply(res, function(x) any(is.na(x)))
  return(structure(list(data.name = input.seq$data.name,
                        n.mrk = length(input.seq$seq.num),
                        seq.num = input.seq$seq.num,
                        pairwise = res,
                        chisq.pval.thres = input.seq$chisq.pval.thres,
                        chisq.pval = input.seq$chisq.pval,
                        nas  = nas),
                   class = "mappoly.twopt"))
}

#' Wrapper function to discrete-based pairwise two-point estimation in C++
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
paralell_pairwise_discrete <- function(mrk.pairs,
                                       input.seq,
                                       geno,
                                       dP,
                                       dQ,
                                       count.cache,
                                       tol = .Machine$double.eps^0.25,
                                       ll = ll)
{
  res <- .Call("pairwise_rf_estimation_disc",
               input.seq$ploidy,
               as.matrix(mrk.pairs),
               as.matrix(geno),
               as.vector(dP),
               as.vector(dQ),
               count.cache$cond,
               tol = tol,
               PACKAGE = "mappoly")
  if(ll) return(res)
  return(lapply(res, format_rf))
}

#' Wrapper function to probability-based pairwise two-point estimation in C++
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
paralell_pairwise_probability <- function(mrk.pairs,
                                          input.seq,
                                          geno,
                                          dP,
                                          dQ,
                                          count.cache,
                                          tol = .Machine$double.eps^0.25,
                                          ll = ll)
{
  res <- .Call("pairwise_rf_estimation_prob",
               input.seq$ploidy,
               as.matrix(mrk.pairs),
               as.integer(dim(geno)),
               as.double(geno),
               as.vector(dP),
               as.vector(dQ),
               count.cache$cond,
               tol = tol,
               PACKAGE = "mappoly")
  if(ll) return(res)
  return(lapply(res, format_rf))
}

#' Format results from pairwise two-point estimation in C++
#'
#' @param void internal function to be documented
#' @keywords internal
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

#' @export
print.mappoly.twopt <- function(x, ...) {
  cat("  This is an object of class 'mappoly.twopt'")
  cat("\n  -----------------------------------------------------")
  ## printing summary
  cat("\n  No. markers:                            ", x$n.mrk, "\n")
  cat("  No. estimated recombination fractions:  ", choose(x$n.mrk, 2) - sum(x$nas))
  cat(" (", round(100 - 100 * sum(x$nas)/length(x$nas), 1), "%)\n", sep = "")
  cat("  -----------------------------------------------------\n")
}

#' @export
plot.mappoly.twopt <- function(x, first.mrk, second.mrk, ...) {
  i <- which(names(x$pairwise)%in%paste(sort(c(first.mrk, second.mrk)), collapse = "-"))
  if(length(i) == 0)
    stop("The requested combination of markers is not included in the two-point object")
  data.name <- x$data.name
  x <- x$pairwise[[i]]
  variable <- value <- sh <- NULL
  x <- reshape2::melt(x, id = rownames(x))
  colnames(x) <- c("sh", "variable", "value")
  rfs <- as.character(format(round(t((subset(x, variable == "rf", select = value))), digits = 2), digits = 2))
  x.temp <- subset(x, variable != "rf")
  mLOD <- max(x.temp$value)
  x.temp <- transform(x.temp, sh = factor(sh, levels = unique(x.temp$sh)))
  ds <- paste0(paste0(get(data.name, pos = 1)$dosage.p1[first.mrk] , "---", get(data.name, pos = 1)$dosage.p1[second.mrk]), "  x  ",
               paste0(get(data.name, pos = 1)$dosage.p2[first.mrk] , "---", get(data.name, pos = 1)$dosage.p2[second.mrk]))
  ggplot2::ggplot(x.temp, ggplot2::aes(sh,    value, label = sh    ,fill = variable))  +
    ggplot2::geom_bar(stat = "identity", position = "identity") +
    ggplot2::annotate("text", x = c(1:length(rfs)), y = rep(mLOD + mLOD/20,length(rfs)), label = rfs, size = 4) +
    ggplot2::annotate("text", x = 1, y = mLOD + mLOD/7, label = "recombination fractions", size = 4, hjust = 0, vjust = 0, fontface = 3)+
    ggplot2::scale_x_discrete(name = paste("Homologous sharing alleles   \n dosage in parents ", ds)) +
    ggplot2::scale_y_continuous(name = "LOD") +
    ggplot2::scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                               name = "",
                               breaks = c("LOD_ph", "LOD_rf"),
                               labels = c("LOD_phase", "LOD_rf"))
}


#' Pairwise two-point analysis - RcppParallel version
#'
#' Performs the two-point pairwise analysis between all markers in a sequence. 
#' For each pair, the function estimates the recombination fraction for all 
#' possible linkage phase configurations and associated LOD Scores.
#' 
#'  Differently from est_pairwise_rf this function returns only the values associated 
#'  to the best linkage phase configuration. 
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#'
#' @param ncpus Number of parallel processes (cores) to spawn (default = 1)
#'
#' @param mrk.pairs a matrix of dimensions 2*N, containing N
#'    pairs of markers to be analyzed. If \code{NULL} (default), all pairs are
#'    considered
#'
#' @param verbose If \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @param tol the desired accuracy. See \code{optimize()} for details
#' 
#' @return An object of class \code{mappoly.twopt2} 
#' 
#' @examples
#'   ## Tetraploid example  
#'   all.mrk <- make_seq_mappoly(tetra.solcap, 100:200)
#'   all.pairs <- est_pairwise_rf2(input.seq = all.mrk, ncpus = 2)
#'   m <- rf_list_to_matrix(all.pairs)
#'   plot(m, fact = 2)
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
#' @export est_pairwise_rf2
#' @importFrom Rcpp sourceCpp
#' @importFrom reshape2 melt acast
#' @importFrom dplyr filter arrange
est_pairwise_rf2 <- function(input.seq, 
                             ncpus = 1L,
                             mrk.pairs = NULL, 
                             verbose = TRUE, 
                             tol = .Machine$double.eps^0.25)
{
  ## checking for correct object
  if (!inherits(input.seq, "mappoly.sequence"))
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  dpl <- duplicated(input.seq$seq.num)
  ## checking for duplicated markers
  if (any(dpl))
    stop("There are duplicated markers in the sequence:\n Check markers: ", unique(input.seq$seq.num[dpl]), " at position(s) ", which(dpl))
  count.cache = cache_counts_twopt(input.seq, cached = TRUE)
  ## get genotypes
  geno <- as.matrix(get(input.seq$data.name, pos = 1)$geno.dose)
  ## all possible pairs
  if (is.null(mrk.pairs)) {
    mrk.pairs <- combn(input.seq$seq.num, 2) - 1
    flag <- 1
  } else {
    mrk.pairs <- mrk.pairs - 1
  }
  ## splitting pairs in chunks
  id <- ceiling(seq(1, (ncol(mrk.pairs) + 1), length.out = 1))
  input.list <- mrk.pairs
  ## parallel version 
  if (verbose) {
    cat("Going RcppParallel mode.\n...")
  }
  RcppParallel::setThreadOptions(numThreads = ncpus)
  count.vector = unlist(count.cache$cond)
  count.phases = unlist(lapply(count.cache$cond, function(x) paste0(names(x), collapse = '/')))
  count.matrix.rownames = unlist(lapply(count.cache$cond, function(x) paste0(rownames(x[[1]]), collapse = '/')))
  count.matrix.number = unlist(lapply(count.cache$cond, length))
  count.matrix.length = unlist(lapply(count.cache$cond, function(x) length(c(unlist(x)) )))
  count.matrix.pos = cumsum(c(1, count.matrix.length[-length(count.matrix.length)]))
  res <- paralell_pairwise_discrete_rcpp(mrk.pairs = as.matrix(input.list),
                                         m = input.seq$ploidy,
                                         geno = geno,
                                         dP = get(input.seq$data.name)$dosage.p1,
                                         dQ = get(input.seq$data.name)$dosage.p2,
                                         count.vector = count.vector,
                                         count.phases = count.phases,
                                         count.matrix.rownames = count.matrix.rownames,
                                         count.matrix.number = count.matrix.number,
                                         count.matrix.pos = count.matrix.pos,
                                         count.matrix.length = count.matrix.length,
                                         tol = tol, threads = ncpus)
  res[res == -1] = NA
  colnames(res) = c("Sh_P1","Sh_P2","rf","LOD_rf","LOD_ph")
  if (verbose) {
    cat("\nDone with pairwise computation.\n~~~~~~~\nReorganizing results\n...")
  }
  if(!exists(x = "flag"))
    return(res)
  seq.mrk.names <- input.seq$seq.mrk.names  
  n <- length(seq.mrk.names)
  Sh_P1 = v_2_m(res[,1], n)
  Sh_P2 = v_2_m(res[,2], n)
  rf = v_2_m(res[,3], n)
  LOD_rf = v_2_m(res[,4], n)
  LOD_ph = v_2_m(res[,5], n)
  invisible(structure(list(seq.mrk.names = seq.mrk.names,
                           data.name  = input.seq$data.name,
                           chisq.pval.thres = input.seq$chisq.pval.thres,
                           chisq.pval = input.seq$chisq.pval,
                           pairwise = list(rf = rf, LOD.rf = LOD_rf,
                                           LOD.ph = LOD_ph, Sh.P1 = Sh_P1, 
                                           Sh.P2 = Sh_P2)), class = "mappoly.twopt2"))
}

#' Wrapper function to discrete-based pairwise two-point estimation in C++
#'
#' @param void internal function to be documented
#' @keywords internal
#' @importFrom RcppParallel RcppParallelLibs
#' @export
paralell_pairwise_discrete_rcpp <- function(mrk.pairs,
                                            m,
                                            geno,
                                            dP,
                                            dQ,
                                            count.vector,
                                            count.phases,
                                            count.matrix.rownames,
                                            count.matrix.number,
                                            count.matrix.pos,
                                            count.matrix.length,
                                            tol = .Machine$double.eps^0.25,
                                            threads)
{
  res <- .Call("pairwise_rf_estimation_disc_rcpp",
               as.matrix(mrk.pairs),
               m,
               as.matrix(geno),
               as.vector(dP),
               as.vector(dQ),
               count.vector,
               count.phases,
               count.matrix.rownames,
               count.matrix.number,
               count.matrix.pos,
               count.matrix.length,
               tol = tol,
               threads = threads,
               PACKAGE = "mappoly")
  return(res)
  ## return(lapply(res, format_rf))
}
