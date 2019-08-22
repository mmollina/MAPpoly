#' Two-point analysis for marker blocks
#'
#' Performs the two-point pairwise analysis between two blcoks of phased SNPs.
#' 
#' When uThe \code{max.inc} argument 
#' 
#'
#' @param block1 a block of phased SNPs contained in a object of
#'               class\code{\link[mappoly]{mappoly.map}}.
#'
#' @param block2 a block of phased SNPs contained in a object of
#'               class\code{\link[mappoly]{mappoly.map}}.
#'
#' @param ph1 an integer indicating the position of
#'        the linkage phase configuration for the first marker block.
#'         If \code{NULL}, it uses the position with maximum likelihood.
#'
#' @param ph2 an integer indicating the position of
#'        the linkage phase configuration for the second marker block.
#'        If \code{NULL}, it uses the position with maximum likelihood.
#'
#' @param M a matrix containg the number of homologous chromosomes that
#'        share alleic variations between the SNPs contained in the phased
#'        marker blocks. The rows represent SNPs in the first marker block
#'        and the columns represent SNPs in the second marker block. This
#'        matrix can be obtained using the function class\code{\link[mappoly]{mat_share}}
#'
#' @param max.inc the maximum inconsitency acceptable when filtering the
#'        linkage phase configuration based on two-point information. 
#'        Pairs with higher numbers will not be taken in to account in 
#'        the filtering process.
#'
#' @param block1.tail the number of SNPs that should be used in the first
#'        phased marker block. If \code{NULL} (default) uses all SNPs
#'
#' @return an object of class \code{mappoly.mrkblock.est.twopt}
#'
#' @examples
#'   \dontrun{
#'   data(hexafake)
#'   mrk.subset<-make_seq_mappoly(hexafake, 21:30)
#'   red.mrk<-elim_redundant(mrk.subset)
#'   unique.mrks<-make_seq_mappoly(red.mrk)
#'   counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'   subset.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                 count.cache = counts.web)
#'   subset.map <- est_rf_hmm_sequential(input.seq = unique.mrks,
#'                                       thres.twopt = 5,
#'                                       thres.hmm = 10,
#'                                       extend.tail = 10,
#'                                       tol = 0.1,
#'                                       tol.final = 10e-3,
#'                                       twopt = subset.pairs,
#'                                       verbose = TRUE,
#'                                       high.prec = FALSE)
#'                                       
#'   s1 <- make_seq_mappoly(hexafake, subset.map$maps[[1]]$seq.num[1:5])                                     
#'   map1 <- get_submap(subset.map, 1:5)
#'   
#'   s2 <- make_seq_mappoly(hexafake, subset.map$maps[[1]]$seq.num[6:10])                                     
#'   map2 <- get_submap(subset.map, 6:10)
#'   
#'    twopt.sub <- make_pairs_mappoly(subset.pairs, 
#'                                    make_seq_mappoly(hexafake, 
#'                                                     c(map1$maps[[1]]$seq.num, 
#'                                                       map2$maps[[1]]$seq.num)))           
#'    M<-mat_share(map1,
#'                 map2,
#'                 twopt.sub,
#'                 count.cache = counts.web,
#'                 thres = 3)
#'    M
#'    rf_map1_map2<-est_rf_marker_blocks(block1 = map1,
#'                                       block2 = map2,
#'                                       ph1 = "best",
#'                                       ph2 = "best",
#'                                       M = M,
#'                                       max.inc = 0,
#'                                       block1.tail = NULL,
#'                                       tol = 0.01)
#'    rf_map1_map2$rf.stats                                   
#'                                
#'    new.map<-subset.map
#'    new.map$maps[[1]]$seq.ph <- list(P = c(rf_map1_map2$phM1[[1]]$config.to.test[[1]]$P,
#'                                           rf_map1_map2$phM2[[1]]$config.to.test[[1]]$P),
#'                                     Q = c(rf_map1_map2$phM1[[1]]$config.to.test[[1]]$Q,
#'                                           rf_map1_map2$phM2[[1]]$config.to.test[[1]]$Q))
#'   new.map$maps[[1]]$seq.rf[25]<-rf_map1_map2$rf.stats[1, "rf"]
#'   plot_compare_haplotypes(m = 6, 
#'                           hom.allele.p1 = subset.map$maps[[1]]$seq.ph$P, 
#'                           hom.allele.p2 = new.map$maps[[1]]$seq.ph$P, 
#'                           hom.allele.q1 = subset.map$maps[[1]]$seq.ph$Q,
#'                           hom.allele.q2 = new.map$maps[[1]]$seq.ph$Q)
#'    }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'
#' @export est_rf_marker_blocks

est_rf_marker_blocks <- function(block1,
                                 block2,
                                 ph1 = "best",
                                 ph2 = "best",
                                 M,
                                 max.inc = 0,
                                 block1.tail = NULL,
                                 tol = 0.01)
{
  ## Checking the input objects
  if (all(is.na(match(class(block1), c("mappoly.map")))))
    stop(deparse(substitute(block1)), " is not an object of class 'mappoly.map'")
  if (all(is.na(match(class(block2), c("mappoly.map")))))
    stop(deparse(substitute(block2)), " is not an object of class 'mappoly.map'")
  
  ## Select all SNPs in first marker block
  if (is.null(block1.tail))
    block1.tail <- block1$info$n.mrk
  
  ## Getting the desired linkage phase configuration for both marker blocks
  if (ph1 == "best")
    ph1 <- which.min(get_LOD(block1, sorted = FALSE))
  if (ph2 == "best")
    ph2 <- which.min(get_LOD(block2, sorted = FALSE))
  
  ## Checking if it is possible to estimate the recombination
  ## fraction given the information present in the blocks
  M1P <- !as.logical(sum(sapply(block1$maps[[1]]$seq.ph$P, function(x) sum(as.logical(x)))))
  M1Q <- !as.logical(sum(sapply(block1$maps[[1]]$seq.ph$Q, function(x) sum(as.logical(x)))))
  M2P <- !as.logical(sum(sapply(block2$maps[[1]]$seq.ph$P, function(x) sum(as.logical(x)))))
  M2Q <- !as.logical(sum(sapply(block2$maps[[1]]$seq.ph$Q, function(x) sum(as.logical(x)))))
  
  ## If it is impossible to compute the recombination fraction, return NA
  if ((M1P && M1Q) || (M2P && M2Q) || (M1P && M2Q) || (M2P && M1Q))
    return(structure(list(rf.stats = NA, phM1 = NA, phM2 = NA), class = c("poly.haplo.est.two.pts")))
  
  M1 <- block1$maps[[ph1]] #map from haplotype 1
  M2 <- block2$maps[[ph2]] #map from haplotype 2
  
  ## Based on a two-point elimination, enumerate all phase configurations
  ## that need to be tested
  w <- generate_all_link_phases_elim_equivalent_haplo(M1, M2, M, block1$info$m, max.inc)
  
  ## Allocating space for all possible maps originated from selected configurations
  res_H0 <- res <- vector("list", nrow(w$wp) * nrow(w$wq))
  phM1 <- phM2 <- vector("list", nrow(w$wp) * nrow(w$wq))
  ct <- 0
  
  # linkage configuration in P
  for (i in 1:nrow(w$wp)) {
    zp <- strsplit(w$wp[i, ], split = "-")
    p1.n <- sapply(zp, function(x) x[1])
    p2.n <- sapply(zp, function(x) x[2])
    # linkage configuration in Q
    for (j in 1:nrow(w$wq)) {
      ct <- ct + 1
      h <- list("vector", 2)
      zq <- strsplit(w$wq[j, ], split = "-")
      q1.n <- sapply(zq, function(x) x[1])
      q2.n <- sapply(zq, function(x) x[2])
      ## Assembling ph to obtain the haplotype index for block1
      ph <- NULL
      ph$m <- block1$info$m
      ph$seq.num <- tail(M1$seq.num, n = block1.tail)
      ph$config.to.test[[1]]$P <- ph_matrix_to_list(w$hP1[as.character(ph$seq.num), p1.n])
      names(ph$config.to.test[[1]]$P) <- ph$seq.num
      ph$config.to.test[[1]]$Q <- ph_matrix_to_list(w$hQ1[as.character(ph$seq.num), q1.n])
      names(ph$config.to.test[[1]]$Q) <- ph$seq.num
      ph$data.name <- block1$info$data.name
      ## getting haplotype indices
      h[[1]] <- haplotype_index(ph = ph)
      phM1[[ct]] <- ph
      ## Assembling ph to obtain the haplotype index for block2
      ph$seq.num <- M2$seq.num
      ph$config.to.test[[1]]$P <- ph_matrix_to_list(w$hP2[, p2.n])
      names(ph$config.to.test[[1]]$P) <- M2$seq.num
      ph$config.to.test[[1]]$Q <- ph_matrix_to_list(w$hQ2[, q2.n])
      names(ph$config.to.test[[1]]$Q) <- M2$seq.num
      ph$data.name <- block1$info$data.name
      h[[2]] <- haplotype_index(ph = ph)
      phM2[[ct]] <- ph
      res[[ct]] <- est_haplo_hmm(m = block1$info$m, 
                                 nmar = length(h),
                                 nind = length(h[[1]]),
                                 haplo = h,
                                 rf_vec = rep(0.01, (length(h) - 1)),
                                 verbose = FALSE,
                                 tol = tol)
      res_H0[[ct]] <- est_haplo_hmm(m = block1$info$m, nmar = length(h),
                                    nind = length(h[[1]]),
                                    haplo = h,
                                    rf_vec = rep(0.5, (length(h) - 1)),
                                    verbose = FALSE,
                                    use_H0 = TRUE,
                                    tol = tol)
    }
  }
  w <- matrix(unlist(res), ncol = 2, byrow = TRUE)
  w_H0 <- matrix(unlist(res_H0), ncol = 2, byrow = TRUE)
  index <- w[, 1] - max(w[, 1])
  res <- matrix(cbind(index, w[, 2], w[,1] - w_H0[,1])[order(index, decreasing = TRUE), ], ncol = 3)
  colnames(res) <- c("ph_LOD", "rf", "rf_LOD")
  res <- data.frame(res)
  return(structure(list(rf.stats = res, phM1 = phM1[order(index, decreasing = TRUE)],
                        phM2 = phM2[order(index, decreasing = TRUE)]),
                   class = c("mappoly.mrkblock.est.twopt")))
}