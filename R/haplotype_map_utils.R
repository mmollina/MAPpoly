#' Return the number of shared homologous chromosomes between all markers in two
#' marker blocks, represented by two objects of type \code{mappoly.map}
#'
#' @param x1 first object of class \code{mappoly.map}
#' @param x2 second object of class \code{mappoly.map}
#' @param ph1 a vector of integers indicating the position of the linkage phase
#'   configuration for \code{x1}. If \code{NULL}, it uses the position with
#'   maximum likelihood.
#' @param ph2 a vector of integers indicating the position of the linkage phase
#'   configuration for \code{x2}. If \code{NULL}, it uses the position with
#'   maximum likelihood.
#' @param count.cache an object of class \code{cache.info} containing
#'   pre-computed genotype frequencies, obtained with
#'   \code{\link[mappoly]{cache_counts_twopt}}.
#' @param thres LOD score threshold for linkage phase configuration.
#'
#' @examples
#'   \dontrun{
#'   data(hexafake)
#'   mrk.subset<-make_seq_mappoly(hexafake, 1:50)
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
#'   s1 <- make_seq_mappoly(hexafake, 
#'                          unique.mrks$seq.mrk.names[1:25])                                     
#'   map1 <- est_rf_hmm_sequential(input.seq = s1,
#'                                 thres.twopt = 5,
#'                                 thres.hmm = 10,
#'                                 extend.tail = 10,
#'                                 tol = 0.01,
#'                                 tol.final = 10e-3,
#'                                 twopt = subset.pairs,
#'                                 verbose = TRUE,
#'                                 high.prec = FALSE)
#'                                 
#'   s2 <- make_seq_mappoly(hexafake, 
#'                          unique.mrks$seq.mrk.names[26:47])                                     
#'   map2 <- est_rf_hmm_sequential(input.seq = s2,
#'                                 thres.twopt = 5,
#'                                 thres.hmm = 10,
#'                                 extend.tail = 10,
#'                                 tol = 0.01,
#'                                 tol.final = 10e-3,
#'                                 twopt = subset.pairs,
#'                                 verbose = TRUE,
#'                                 high.prec = FALSE)
#'    twopt.sub <- make_pairs_mappoly(subset.pairs, 
#'                                    make_seq_mappoly(hexafake, 
#'                                                     c(map1$maps[[1]]$seq.num, 
#'                                                       map2$maps[[1]]$seq.num)))           
#'    M<-mat_share(map1,
#'                 map2,
#'                 twopt,
#'                 count.cache = counts.web,
#'                 thres = -1)
#'    M
#'    }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @keywords internal
#' @export mat_share

mat_share <-
  function(x1,
           x2,
           twopt,
           ph1 = "best",
           ph2 = "best",
           count.cache = NULL,
           thres) {
    if (ph1 == "best")
      ph1 <- which.min(get_LOD(x1))
    if (ph2 == "best")
      ph2 <- which.min(get_LOD(x2))
    s1 <- x1$maps[[ph1]]$seq.num
    s2 <- x2$maps[[ph2]]$seq.num
    M <- list(P = matrix(NA, length(x1$maps[[ph1]]$seq.num), length(x2$maps[[ph2]]$seq.num)),
              Q = matrix(NA, length(x1$maps[[ph1]]$seq.num), length(x2$maps[[ph2]]$seq.num)))
    dimnames(M$P) <-
      list(x1$maps[[ph1]]$seq.num, x2$maps[[ph2]]$seq.num)
    dimnames(M$Q) <-
      list(x1$maps[[ph1]]$seq.num, x2$maps[[ph2]]$seq.num)
    for (i in x1$maps[[ph1]]$seq.num) {
      for (j in x2$maps[[ph2]]$seq.num) {
        w <- twopt$pairwise[paste(sort(c(i, j)), collapse = "-")][[1]]
        if (length(w) > 3) {
          if (abs(w[2, "LOD_ph"]) >= thres) {
            temp <- as.numeric(unlist(strsplit(rownames(w)[1], "-")))
            M$P[as.character(i), as.character(j)] <- temp[1]
            M$Q[as.character(i), as.character(j)] <- temp[2]
          }
        }
      }
    }
    return(M)
  }

#' Generates all possible linkage phases between two blocks of markers,
#' eliminating equivalent configurations, i.e., configurations with the
#' same likelihood
#'
#' @param block1 submap with markers of the first block
#' @param block2 submap with markers of the second block
#' @param rf.matrix matrix obtained with the function \code{rf_list_to_matrix}
#' using the parameter \code{shared.alleles = TRUE}
#' @param m ploidy level (i.e. 4, 6 and so on)
#' @param max.inc maximum number of allowed inconsistencies
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} and Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @export generate_all_link_phases_elim_equivalent_haplo
generate_all_link_phases_elim_equivalent_haplo <-
    function(block1, block2, rf.matrix, m, max.inc = NULL) {
        ## Getting M matrix
        M = list(P = rf.matrix$ShP[as.character(block1$seq.num),as.character(block2$seq.num)], Q = rf.matrix$ShQ[as.character(block1$seq.num),as.character(block2$seq.num)])
        
        ## Parent P: all permutations between blocks
        hP1 <- ph_list_to_matrix(L = block1$seq.ph$P, m = m)
        p1 <- apply(hP1, 2, paste, collapse = "")
        hP2 <- ph_list_to_matrix(L = block2$seq.ph$P, m = m)
        p2 <- apply(hP2, 2, paste, collapse = "")
        dimnames(hP1) <- list(block1$seq.num, p1)
        dimnames(hP2) <- list(block2$seq.num, p2)
        p2 <- perm_tot(p2)
        
        ## Parent Q: all permutations between blocks
        hQ1 <- ph_list_to_matrix(L = block1$seq.ph$Q, m = m)
        q1 <- apply(hQ1, 2, paste, collapse = "")
        hQ2 <- ph_list_to_matrix(L = block2$seq.ph$Q, m = m)
        q2 <- apply(hQ2, 2, paste, collapse = "")
        dimnames(hQ1) <- list(block1$seq.num, q1)
        dimnames(hQ2) <- list(block2$seq.num, q2)
        q2 <- perm_tot(q2)

        ## WP: removing redundancy and accounting for shared elleles
        wp <- NULL
        for (i in 1:nrow(p2))
            wp <- rbind(wp, paste(p1, p2[i, ], sep = "-"))
        wp.n <- unique(t(apply(unique(wp), 1, sort)))
        yp <- apply(t(apply(wp, 1, function(x) sort(x))), 1, paste, collapse = "|")
        yp.n <- apply(wp.n, 1, function(x) paste(sort(x), collapse = "|"))
        wp <- wp[match(yp.n, yp), , drop = FALSE]
        ct <- numeric(nrow(wp))
        for (i in 1:nrow(wp)){
            a = matrix(unlist(strsplit(wp[i, ], "-")), ncol = 2, byrow = TRUE)
            sharedP = hP1[,a[,1]]%*%t(hP2[,a[,2]])
            ct[i] = sum((M$P != sharedP), na.rm = T)
        }
        ## Checking inconsistency
        if(is.null(max.inc)){
            id <- rep(TRUE, length(ct))      
        } else
            id <- ct <= max.inc
        ## Maximum inconsistency
        if (sum(id) == 0)
            id <- which.min(ct)
        wp <- matrix(wp[id, ], ncol = m)

        ## WQ: removing redundancy and accounting for shared elleles
        wq <- NULL
        for (i in 1:nrow(q2))
            wq <- rbind(wq, paste(q1, q2[i, ], sep = "-"))
        wq.n <- unique(t(apply(unique(wq), 1, sort)))
        yq <- apply(t(apply(wq, 1, function(x) sort(x))), 1, paste, collapse = "|")
        yq.n <- apply(wq.n, 1, function(x) paste(sort(x), collapse = "|"))
        wq <- wq[match(yq.n, yq), , drop = FALSE]
        ct <- numeric(nrow(wq))
        for (i in 1:nrow(wq)){
            a = matrix(unlist(strsplit(wq[i, ], "-")), ncol = 2, byrow = TRUE)
            sharedQ = hQ1[,a[,1]]%*%t(hQ2[,a[,2]])
            ct[i] = sum((M$Q != sharedQ), na.rm = T)
        }
        ## Checking inconsistency
        if(is.null(max.inc)){
            id <- rep(TRUE, length(ct))      
        } else
            id <- ct <= max.inc
        ## Maximum inconsistency
        if (sum(id) == 0)
            id <- which.min(ct)
        wq <- matrix(wq[id, ], ncol = m)
        
        ## Re-arranging phases
        phase.to.test <- vector("list", nrow(wp) * nrow(wq))
        cte <- 1
        for(i in 1:nrow(wp)){
            for(j in 1:nrow(wq)){
                phase.to.test[[cte]] <- list(P = ph_matrix_to_list(hP2[,sapply(strsplit(wp[i,], "-"), function(x) x[2])]),
                                             Q = ph_matrix_to_list(hQ2[,sapply(strsplit(wq[j,], "-"), function(x) x[2])]))
                cte <- cte + 1
            }
        }
        phase.to.test <- unique(phase.to.test)
        names(phase.to.test) <- paste0("config.", 1:length(phase.to.test))
        return(phase.to.test)
    }

#' Generates all possible linkage phases between two blocks of markers,
#' eliminating equivalent configurations, i.e., configurations with the
#' same likelihood (old function)
#'
#' @param block1 submap with markers of the first block
#' @param block2 submap with markers of the second block
#' @param M matrix obtained with \code{mate_share} function
#' @param m ploidy level (i.e. 4, 6 and so on)
#' @param max.inc maximum number of allowed inconsistencies
#' @keywords internal
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export generate_all_link_phases_elim_equivalent_haplo_old
generate_all_link_phases_elim_equivalent_haplo_old <-
  function(block1, block2, M, m, max.inc = NULL) {
    ## Parent P
    hP1 <- ph_list_to_matrix(L = block1$seq.ph$P, m = m)
    p1 <- apply(hP1, 2, paste, collapse = "")
    hP2 <- ph_list_to_matrix(L = block2$seq.ph$P, m = m)
    p2 <- apply(hP2, 2, paste, collapse = "")
    dimnames(hP1) <- list(block1$seq.num, p1)
    dimnames(hP2) <- list(block2$seq.num, p2)
    ## Here I am anchoring the linkage phase configuration in block1 and listing all possible permutations in hap2
    p2 <- perm_tot(p2)
    ## Parent Q
    hQ1 <- ph_list_to_matrix(L = block1$seq.ph$Q, m = m)
    q1 <- apply(hQ1, 2, paste, collapse = "")
    hQ2 <- ph_list_to_matrix(L = block2$seq.ph$Q, m = m)
    q2 <- apply(hQ2, 2, paste, collapse = "")
    dimnames(hQ1) <- list(block1$seq.num, q1)
    dimnames(hQ2) <- list(block2$seq.num, q2)
    q2 <- perm_tot(q2)
    wp <- NULL
    for (i in 1:nrow(p2))
      wp <- rbind(wp, paste(p1, p2[i, ], sep = "-"))
    wp.n <- unique(t(apply(unique(wp), 1, sort)))
    yp <- apply(t(apply(wp, 1, function(x) sort(x))), 1, paste, collapse = "|")
    yp.n <- apply(wp.n, 1, function(x) paste(sort(x), collapse = "|"))
    wp <- wp[match(yp.n, yp), , drop = FALSE]
    ct <- numeric(nrow(wp))
    #####################
    ## Speed up this loop
    #####################
    for (i in 1:nrow(wp)) {
      a <- matrix(unlist(strsplit(wp[i, ], "-")), ncol = 2, byrow = TRUE)
      for (j in 1:nrow(M$P)) {
        for (k in 1:ncol(M$P)) {
          a1 <- M$P[as.character(block1$seq.num[j]), as.character(block2$seq.num[k])]
          a2 <- sum((hP1[as.character(block1$seq.num[j]), 
                         a[, 1]] + hP2[as.character(block2$seq.num[k]), a[, 2]]) == 2)
          if (is.na(a1))
            a1 <- a2
          if (a1 != a2)
            ct[i] <- ct[i] + 1
        }
      }
    }
    if(is.null(max.inc)){
      id <- rep(TRUE, length(ct))      
    } else
      id <- ct <= max.inc
    #maximum inconsistency
    if (sum(id) == 0)
      id <- which.min(ct)
    wp <- matrix(wp[id, ], ncol = m)
    ## Choosing unique configurations among all possible permutations (Q)
    wq <- NULL
    for (i in 1:nrow(q2))
      wq <- rbind(wq, paste(q1, q2[i, ], sep = "-"))
    wq.n <- unique(t(apply(unique(wq), 1, sort)))
    yq <-
      apply(t(apply(wq, 1, function(x)
        sort(x))), 1, paste, collapse = "|")
    yq.n <-
      apply(wq.n, 1, function(x)
        paste(sort(x), collapse = "|"))
    wq <- wq[match(yq.n, yq), , drop = FALSE]
    ct <- numeric(nrow(wq))
    for (i in 1:nrow(wq)) {
      a <- matrix(unlist(strsplit(wq[i, ], "-")), ncol = 2, byrow = TRUE)
      for (j in 1:nrow(M$Q)) {
        for (k in 1:ncol(M$Q)) {
          a1 <- M$Q[as.character(block1$seq.num[j]), as.character(block2$seq.num[k])]
          a2 <-
            sum((hQ1[as.character(block1$seq.num[j]), a[, 1]] + hQ2[as.character(block2$seq.num[k]), a[, 2]]) == 2)
          if (is.na(a1))
            a1 <- a2
          if (a1 != a2)
            ct[i] <- ct[i] + 1
        }
      }
    }
    if(is.null(max.inc)){
      id <- rep(TRUE, length(ct))      
    } else
      id <- ct <= max.inc
    if (sum(id) == 0)
      id <- which.min(ct)
    wq <- matrix(wq[id, ], ncol = m)
    ## Re-arraging phases
    phase.to.test <- vector("list", nrow(wp) * nrow(wq))
    cte <- 1
    for(i in 1:nrow(wp)){
      for(j in 1:nrow(wq)){
        phase.to.test[[cte]] <- list(P = ph_matrix_to_list(hP2[,sapply(strsplit(wp[i,], split = "-"), function(x) x[2])]),
                                     Q = ph_matrix_to_list(hQ2[,sapply(strsplit(wq[j,], split = "-"), function(x) x[2])]))
        cte <- cte + 1
      }
    }
    phase.to.test <- unique(phase.to.test)
    names(phase.to.test) <- paste0("config.", 1:length(phase.to.test))
    return(phase.to.test)
  }

#' Generates a list where each element represents one individual. Each
#' element of the list contains the states that should be visited
#' given a marker block
#'
#' @param void intern function to be documented
#' @keywords internal
#' @export haplotype_index
haplotype_index <- function(ph, conf = 1) {
  if (conf > length(ph$config.to.test))
    warning(
      "there is ",
      length(ph$config.to.test),
      " configuration(s) to test. Using the first one."
    )
  P <- sapply(ph$config.to.test[[conf]]$P, function(x, m) {
    v <- numeric(m)
    v[x] <- 1
    v
  }, m = ph$m)
  Q <- sapply(ph$config.to.test[[conf]]$Q, function(x, m) {
    v <- numeric(m)
    v[x] <- 1
    v
  }, m = ph$m)
  A <- map_haplo(ph$m, P, Q)
  B <- get(ph$data.name, pos=1)$geno.dose[ph$seq.num,]
  if(any(B == ph$m + 1))
    B[B == ph$m + 1] <- NA 
  ##
  obs <- apply(B, 2, paste, collapse = "")
  stat <- apply(A[, -c(1, 2)], 1, paste, collapse = "")
  H <- vector("list", length(obs))
  for (i in 1:length(H)){
    H[[i]] <-
      matrix(A[which(stat == obs[i]), c(1, 2)], ncol = 2) - 1
    #This would be to include probabilities associated with the genotypic states
    #H[[i]] <-  cbind(H[[i]],rep(1, nrow(H[[i]]))) 
  }
  for (i in grep("NA", obs))
  {
    H[[i]] <-
      as.matrix(expand.grid(0:(choose(ph$m, ph$m / 2) - 1), 0:(choose(ph$m, ph$m /
                                                                        2) - 1), stringsAsFactors = FALSE)[2:1])
    #H[[i]] <-  cbind(H[[i]],rep(1, nrow(H[[i]])))
  }
  ## FIXME: instead of using NA in case of c.o., use the minimum tail with no c.o. 
  for (i in   which(sapply(H, length)==0))
  {
    H[[i]] <-
      as.matrix(expand.grid(0:(choose(ph$m, ph$m / 2) - 1), 0:(choose(ph$m, ph$m /
                                                                        2) - 1), stringsAsFactors = FALSE)[2:1])
    #H[[i]] <-  cbind(H[[i]],rep(1, nrow(H[[i]])))
  }
  H
}

#' States that should be visited given a molecular phenotype of the
#' haplotype
#'
#' @param m ploidy level
#' @param P a list of vectors indicating the haplotype configuration
#'     for genitor P
#' @param Q a list of vectors indicating the haplotype configuration
#'     for genitor Q
#' @return a matrix whose the first two columns indicates all the
#'     possible states and the remaining columns indicates the
#'     observed molecular phenotype of the haplotype
#' @keywords internal
#' @export map_haplo
map_haplo <- function(m, P, Q) {
  ngen <- choose(m, m / 2)
  res <- matrix(NA, ncol = 2 + ncol(P), nrow = ngen ^ 2)
  idx <- combn(m, m / 2)
  
  if (ncol(P) == 1) {
    for (i in 1:ncol(idx)) {
      temp <- sum(P[idx[, i], ])
      for (j in 1:ncol(idx)) {
        res[(i - 1) * ngen + j, ] <- c(i, j, temp + sum(Q[idx[, j], ]))
      }
    }
    return(res)
  }
  for (i in 1:ncol(idx)) {
    temp <- apply(P[idx[, i], ], 2, sum)
    for (j in 1:ncol(idx)) {
      res[(i - 1) * ngen + j, ] <-
        c(i, j, temp + apply(Q[idx[, j], ], 2, sum))
    }
  }
  return(res)
}


#' Estimate a genetic map given a sequence of block markers
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
est_haplo_hmm <-
  function(m, nmar, nind, haplo, emit = NULL, 
           rf_vec, verbose, use_H0 = FALSE, tol) {
    ## In case no genotypic probabilities distrubutions are provided
    if(is.null(emit)){
      emit <- vector("list", length(haplo))
      for(i in  1:length(haplo)){
        tempemit <- vector("list", length(haplo[[i]]))
        for(j in 1:length(haplo[[i]])){
          tempemit[[j]] <- rep(1, nrow(haplo[[i]][[j]]))
        }
        emit[[i]] <- tempemit
      }
    }
    res.temp <-
      .Call("est_haplotype_map",
            m,
            nmar,
            nind,
            haplo,
            emit,
            rf_vec,
            verbose,
            tol,
            use_H0,
            PACKAGE = "mappoly")
    res.temp
  }

