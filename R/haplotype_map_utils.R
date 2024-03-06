#' Eliminate equivalent linkage phases
#' 
#' Generates all possible linkage phases between two blocks of markers
#' (or a block and a marker), eliminating equivalent configurations, 
#' i.e. configurations with the same likelihood and also considering
#' the two-point information (shared alleles)
#'
#' @param block1 submap with markers of the first block
#' 
#' @param block2 submap with markers of the second block, 
#' or just a single marker identified by its name
#'  
#' @param rf.matrix matrix obtained with the function \code{rf_list_to_matrix}
#' using the parameter \code{shared.alleles = TRUE}
#' 
#' @param ploidy ploidy level (i.e. 4, 6 and so on)
#' 
#' @param max.inc maximum number of allowed inconsistencies (default = NULL: don't check inconsistencies)
#' 
#' @keywords internal
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} and Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#' 
#' @export generate_all_link_phases_elim_equivalent_haplo
#' 
generate_all_link_phases_elim_equivalent_haplo <- 
  function(block1, block2, rf.matrix, ploidy, max.inc = NULL) {
    ## Check block2 class (block or single marker)
    if (length(block2)==1){
      dp <- get(rf.matrix$data.name)$dosage.p1[block2]
      dq <- get(rf.matrix$data.name)$dosage.p2[block2]
      if(dp != 0) dp <- 1:dp
      if(dq != 0) dq <- 1:dq
      seq.ph = list(P = list(dp), Q = list(dq))
      block2 = list(seq.ph = seq.ph, mrk.names = block2)
    }
    ## Getting M matrix
    M = list(P = rf.matrix$ShP[block1$mrk.names,block2$mrk.names], 
             Q = rf.matrix$ShQ[block1$mrk.names,block2$mrk.names])
    
    ## Parent P: all permutations between blocks
    hP1 <- ph_list_to_matrix(L = block1$seq.ph$P, ploidy = ploidy)
    p1 <- apply(hP1, 2, paste, collapse = "")
    hP2 <- ph_list_to_matrix(L = block2$seq.ph$P, ploidy = ploidy)
    p2 <- apply(hP2, 2, paste, collapse = "")
    dimnames(hP1) <- list(block1$mrk.names, p1)
    dimnames(hP2) <- list(block2$mrk.names, p2)
    p2 <- perm_tot(p2)
    
    ## Parent Q: all permutations between blocks
    hQ1 <- ph_list_to_matrix(L = block1$seq.ph$Q, ploidy = ploidy)
    q1 <- apply(hQ1, 2, paste, collapse = "")
    hQ2 <- ph_list_to_matrix(L = block2$seq.ph$Q, ploidy = ploidy)
    q2 <- apply(hQ2, 2, paste, collapse = "")
    dimnames(hQ1) <- list(block1$mrk.names, q1)
    dimnames(hQ2) <- list(block2$mrk.names, q2)
    q2 <- perm_tot(q2)
    
    ## WP: removing redundancy and accounting for shared alleles
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
      sharedP = tcrossprod(hP1[,a[,1], drop = FALSE], hP2[,a[,2], drop = FALSE])            
      ct[i] = sum((M$P != sharedP), na.rm = T)
    }
    ## Checking inconsistency
    if(is.null(max.inc)){
      id <- rep(TRUE, length(ct))      
    } else
      id <- ct <= max.inc
    ## Maximum inconsistency
    if (sum(id)  ==  0)
      id <- which.min(ct)
    wp <- matrix(wp[id, ], ncol = ploidy)
    
    ## WQ: removing redundancy and accounting for shared alleles
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
      sharedQ = tcrossprod(hQ1[,a[,1], drop = FALSE], hQ2[,a[,2], drop = FALSE])
      ct[i] = sum((M$Q != sharedQ), na.rm = T)
    }
    ## Checking inconsistency
    if(is.null(max.inc)){
      id <- rep(TRUE, length(ct))      
    } else
      id <- ct <= max.inc
    ## Maximum inconsistency
    if (sum(id)  ==  0)
      id <- which.min(ct)
    wq <- matrix(wq[id, ], ncol = ploidy)
    
    ## Re-arranging phases
    phase.to.test <- vector("list", nrow(wp) * nrow(wq))
    cte <- 1
    for(i in 1:nrow(wp)){
      for(j in 1:nrow(wq)){
        P <- ph_matrix_to_list(hP2[,sapply(strsplit(wp[i,], "-"), function(x) x[2]), drop = FALSE])
        Q <- ph_matrix_to_list(hQ2[,sapply(strsplit(wq[j,], "-"), function(x) x[2]), drop = FALSE])
        names(P) <- names(Q) <- block2$seq.num
        phase.to.test[[cte]] <- list(P = P, Q = Q)
        cte <- cte + 1
      }
    }
    phase.to.test <- unique(phase.to.test)
    names(phase.to.test) <- paste0("config.", 1:length(phase.to.test))
    return(phase.to.test)
  }

#' Estimate a genetic map given a sequence of block markers
#' @keywords internal
est_haplo_hmm  <- 
  function(ploidy, n.mrk, n.ind, haplo, emit = NULL, 
           rf_vec, verbose = TRUE, use_H0 = FALSE, 
           highprec = FALSE, tol = 10e-4) {
    ## Checking capabilities
    if (verbose && !capabilities("long.double") && highprec){
      cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
      highprec = FALSE
    }
    ## In case no genotypic probabilities distributions are provided
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
    if(highprec){
      res.temp  <- 
        .Call("est_haplotype_map_highprec",
              ploidy,
              n.mrk,
              n.ind,
              haplo,
              emit,
              rf_vec,
              verbose,
              tol,
              use_H0,
              PACKAGE = "mappoly")
      return(res.temp)
      
    } else {
      res.temp  <- 
        .Call("est_haplotype_map",
              ploidy,
              n.mrk,
              n.ind,
              haplo,
              emit,
              rf_vec,
              verbose,
              tol,
              use_H0,
              PACKAGE = "mappoly")
      return(res.temp)
    }
  }

#' Estimate a genetic map given a sequence of block markers 
#' given the conditional probabilities of the genotypes
#' @keywords internal
est_map_haplo_given_genoprob <- function(map.list,
                                       genoprob.list,
                                       tol = 10e-5){
  ## Checking capabilities
  if (!capabilities("long.double")){
    message("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
  }  
  ploidy <- map.list[[1]]$info$ploidy
  ## number of genotypic states
  ngam <- choose(ploidy, ploidy/2)
  ## Number of genotypes in the offspring
  n.gen <- ngam^2
  ## number of markers
  n.mrk <- sapply(map.list, function(x) length(x$info$seq.num))
  ## number of individuals
  n.ind <- dim(genoprob.list[[1]]$probs)[3]
  ## the threshold for visiting states: 1/n.gen
  thresh.cut.path <- 1/n.gen
  
  ## Hash table: homolog combination --> states to visit in both parents
  A <- as.matrix(expand.grid(0:(ngam-1), 
                           0:(ngam-1))[,2:1])
  rownames(A) <- dimnames(genoprob.list[[1]]$probs)[[1]]
  ## h: states to visit in both parents
  ## e: probability distribution 
  h <- e <- NULL
  for(j in 1:length(map.list)){
    e.temp <- h.temp <- vector("list", n.ind)
    for(i in 1:n.ind){
      a <- genoprob.list[[j]]$probs[,dim(genoprob.list[[j]]$probs)[2],i]  
      e.temp[[i]] <- a[a > thresh.cut.path]
      h.temp[[i]] <- A[names(e.temp[[i]]), , drop = FALSE]
    }
    h <- c(h, list(h.temp))
    e <- c(e, list(e.temp))    
  }
  map <- est_haplo_hmm(ploidy = ploidy, 
                     n.mrk = length(h), 
                     n.ind = n.ind, 
                     haplo = h, 
                     emit = e, 
                     rf_vec = rep(0.01, length(h)-1), 
                     verbose = FALSE, 
                     use_H0 = FALSE, 
                     tol = tol)
  genoprob <- calc_genoprob_haplo (ploidy = ploidy, 
                                 n.mrk = length(h), 
                                 n.ind = n.ind, 
                                 haplo = h, 
                                 emit = e,  
                                 rf_vec = map[[2]],
                                 ind.names = dimnames(genoprob.list[[1]]$probs)[[3]],
                                 verbose = FALSE)
  list(map = map, genoprob = genoprob)
}

#' Compute conditional probabilities of the genotypes given a sequence 
#' of block markers
#' @keywords internal
calc_genoprob_haplo <- function(ploidy, n.mrk, n.ind, haplo, emit = NULL, 
                                rf_vec, ind.names, verbose = TRUE, 
                                highprec = FALSE) {
  ## Checking capabilities
  if (verbose && !capabilities("long.double") && highprec){
    cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
    highprec = FALSE
  }
  ## In case no genotypic probabilities distributions are provided
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
  mrk.names <- names(haplo)
  if(highprec){
    res.temp <- calc_genprob_haplo_highprec_cpp(ploidy,
                                                n.mrk,
                                                n.ind,
                                                haplo,
                                                emit,
                                                rf_vec,
                                                as.numeric(rep(0, choose(ploidy, ploidy/2)^2 * n.mrk * n.ind)),
                                                verbose)
  } else{
    res.temp <- calc_genprob_haplo_cpp(ploidy,
                                       n.mrk,
                                       n.ind,
                                       haplo,
                                       emit,
                                       rf_vec,
                                       as.numeric(rep(0, choose(ploidy, ploidy/2)^2 * n.mrk * n.ind)),
                                       verbose)
  }
  if(verbose) cat("\n")
  dim(res.temp[[1]]) <- c(choose(ploidy,ploidy/2)^2,n.mrk,n.ind)
  dimnames(res.temp[[1]]) <- list(kronecker(apply(combn(letters[1:ploidy],ploidy/2),2, paste, collapse = ""),
                                          apply(combn(letters[(ploidy+1):(2*ploidy)],ploidy/2),2, paste, collapse = ""), paste, sep = ":"),
                                mrk.names, ind.names)
  structure(list(probs = res.temp[[1]], map = rf_vec), class = "mappoly.genoprob")
}
