#' Eliminates equivalent linkage phase configurations
#' 
#' Drop equivalent linkage phase configurations, i.e. the ones which
#' have permuted homologous chromosomes
#'
#' @param Z a list of matrices whose columns represent homologous
#'     chromosomes and the rows represent markers
#' 
#' @return a unique list of matrices
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @keywords internal
elim_equiv <- function(Z) {
  a <- sapply(Z, function(x) paste(sort(apply(x, 2, paste, collapse = "")), collapse = ""))
  Z[which(!duplicated(a))]
}

#' Given a pair of character indicating the numbers i and j : 'i-j',
#' returns a numeric pair c(i,j)
#'
#' @param w a pair of characters 'i-j'
#' @return a numeric pair c(i,j)
#' @keywords internal
get_ij <- function(w) as.numeric(unlist(strsplit(w, split = "-")))


#' Linkage phase format conversion: matrix to list
#' 
#' This function converts linkage phase configurations from matrix
#' form to list
#'
#' @param M matrix whose columns represent homologous chromosomes and
#'     the rows represent markers
#' 
#' @return a list of linkage phase configurations
#' 
#' @keywords internal
#' @export
ph_matrix_to_list <- function(M) {
  w <- lapply(split(M, seq(NROW(M))), function(x, M) which(x  ==  1))
  w[sapply(w, function(x) length(x)  ==  0)] <- 0
  w
}

#' Linkage phase format conversion: list to matrix
#' 
#' This function converts linkage phase configurations from list
#' to matrix form
#'
#' @param L a list of configuration phases
#' 
#' @param ploidy ploidy level
#' 
#' @return a matrix whose columns represent homologous chromosomes and
#'     the rows represent markers
#' 
#' @keywords internal
#' @export
ph_list_to_matrix <- function(L, ploidy) {
  M <- matrix(0, nrow = length(L), ncol = ploidy)
  for (i in 1:nrow(M)) if (all(L[[i]] != 0))
    M[i, L[[i]]] <- 1
  M
}

#' Get the indices of selected linkage phases given a threshold
#'
#' @param x a data frame containing information about two markers. In
#'     this data frame, the lines indicate the possible configuration
#'     phases and the columns indicate the LOD for configuration phase
#'     (ph_LOD), the recombination fraction (rf), and the LOD for
#'     recombination fraction (rf_LOD) 
#'     
#' @param thres a threshold from which the linkage phases can be
#'     discarded (if abs(ph_LOD) > thres)
#' @return a list of indices for both parents
#' @keywords internal
get_indices_from_selected_phases <- function(x, thres) {
  y <- rownames(x)[which(abs(x[, 1]) <= abs(thres))]
  y <- matrix(as.numeric(unlist(strsplit(y, split = "-"))), ncol = 2, byrow = TRUE)
  list(P = unique(y[, 1]), Q = unique(y[, 2]))
}

#' Generate all possible linkage phases in matrix form given the dose
#' and the number of shared alleles between a inserted marker and a
#' pre-computed linkage configuration.
#'
#' @param X a list of matrices whose columns represent homologous
#'     chromosomes and the rows represent markers. Each element of the
#'     list represents a linkage phase configuration.
#' @param d the dosage of the inserted marker
#' @param sh the number of shared alleles between k1 (marker already
#'     present on the sequence) and k2 (the inserted marker)
#' @param ploidy the ploidy level
#' @param k1 marker already present in the sequence
#' @param k2 inserted marker
#' @return a unique list of matrices representing linkage phases
#' @keywords internal
generate_all_link_phase_elim_equivalent <- function(X, d, sh, ploidy, k1, k2) {
  mat <- matrix(0, nrow = ploidy, ncol = choose(ploidy, d))
  i <- combn(ploidy, d)
  j <- cumsum(c(0, rep(ploidy, choose(ploidy, d) - 1)))
  mat[as.numeric(apply(i, 1, function(x, y) x + y, y = j))] <- 1
  ct <- NULL
  Y <- NULL
  for (i in 1:ncol(mat)) {
    Y[[i]] <- rbind(X, mat[, i])
    so <- sum(apply(Y[[i]][c(k1, k2), ], 2, sum)  ==  2)
    if (any(so  ==  sh))
      ct <- c(ct, i)
  }
  elim_equiv(Y[ct])
}

#' Concatenate new marker
#' 
#' Inserts a new marker at the end of the sequence, taking into account
#' the two-point information
#'
#' @param X a list of matrices whose columns represent homologous
#'     chromosomes and the rows represent markers
#' @param d the dosage of the inserted marker
#' @param sh a list of shared alleles between all markers in the sequence
#' @param seq.num a vector of integers containing the number of each marker in the raw data file
#' @param ploidy the ploidy level
#' @param mrk the marker to be inserted
#' @return a unique list of matrices representing linkage phases
#'
#' @keywords internal
concatenate_new_marker <- function(X = NULL, d, sh = NULL, seq.num = NULL, ploidy, mrk = 1) {
  if (is.null(X) & is.null(sh) & mrk  ==  1 & is.null(seq.num)) {
    Y <- numeric(ploidy)
    if (d != 0)
      Y[1:d] <- 1
    return(list(matrix(Y, ncol = ploidy)))
  }
  Y.final <- NULL
  for (i in 1:length(X)) {
    for (j in (mrk - 1):1) {
      id <- paste(sort(c(seq.num[j], seq.num[mrk])), collapse = "-")
      Y.temp <- generate_all_link_phase_elim_equivalent(X[[i]], d = d[mrk], sh = sh[[id]], ploidy = ploidy, k1 = j, k2 = mrk)
      if (length(Y.temp)  ==  1 && j  ==  (mrk - 1)) {
        Y <- Y.temp
        best.conf <- 1
        (break)()
      }
      if (j  ==  (mrk - 1)) {
        Y <- Y.temp
        best.conf <- 1:length(Y)
      } else {
        best.conf <- na.omit(match(Y, Y.temp))
        if (length(best.conf) != 0)
          Y <- Y.temp[best.conf]
        if (length(Y)  ==  1)
          (break)()
      }
    }
    Y.final <- c(Y.final, Y)
  }
  Y.final
}

#' Eliminate configurations using two-point information
#' 
#' Drops unlikely configuration phases given the two-point information
#' and a LOD threshold
#'
#' @param input.seq an object of class \code{mappoly.sequence}.
#' @param twopt an object of class \code{mappoly.twopt}
#' @param thres threshold from which the linkage phases can be
#'     discarded (if abs(ph_LOD) > thres)
#' @return a unique list of matrices representing linkage phases
#' @keywords internal
elim_conf_using_two_pts <- function(input.seq, twopt, thres) {
  if (!inherits(input.seq, "mappoly.sequence"))
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  check <- check_pairwise(input.seq, twopt)
  if (any(check != 0)){
    stop("There is no information for pairs: \n", paste(capture.output(print(check)), collapse = "\n"))
  }
  info.par <-detect_info_par(input.seq)
  index <- apply(apply(combn(input.seq$seq.num, 2), 2, sort), 2, paste, collapse = "-")
  w <- twopt$pairwise[index]
  sh <- lapply(w, get_indices_from_selected_phases, thres = thres)
  dp <- input.seq$seq.dose.p1
  dq <- input.seq$seq.dose.p2
  ploidy <- input.seq$ploidy
  if(info.par == "both")
  {
    dp <- input.seq$seq.dose.p1
    dq <- input.seq$seq.dose.p2
    sp <- lapply(sh, function(x) x$P)
    sq <- lapply(sh, function(x) x$Q)
    XP <- concatenate_new_marker(d = dp[1], ploidy = input.seq$ploidy)
    for (i in 2:length(dp)) {
      if (is.null(XP)) {
        warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
        return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
      }
      XP <- concatenate_new_marker(X = XP, d = dp, sh = sp, seq.num = input.seq$seq.num, ploidy = input.seq$ploidy, mrk = i)
    }
    if (is.null(XP)) {
      warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
      return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
    } else for (i in 1:length(XP)) dimnames(XP[[i]]) <- list(input.seq$seq.num, paste("h", 1:input.seq$ploidy, sep = ""))
    XQ <- concatenate_new_marker(d = dq[1], ploidy = input.seq$ploidy)
    for (i in 2:length(dq)) {
      if (is.null(XP)) {
        warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
        return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
      }
      XQ <- concatenate_new_marker(X = XQ, d = dq, sh = sq, seq.num = input.seq$seq.num, ploidy = input.seq$ploidy, mrk = i)
    }
    if (is.null(XQ)) {
      warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
      return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
    } else for (i in 1:length(XQ)) dimnames(XQ[[i]]) <- list(input.seq$seq.num, paste("h", 1:input.seq$ploidy, sep = ""))
    return(list(P = XP, Q = XQ))
  } else if (info.par == "p1") 
  {
    sp <- lapply(sh, function(x) x$P)
    XP <- concatenate_new_marker(d = dp[1], ploidy = input.seq$ploidy)
    for (i in 2:length(dp)) {
      if (is.null(XP)) {
        warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
        return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
      }
      XP <- concatenate_new_marker(X = XP, d = dp, sh = sp, seq.num = input.seq$seq.num, ploidy = input.seq$ploidy, mrk = i)
    }
    if (is.null(XP)) {
      warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
      return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
    } else for (i in 1:length(XP)) dimnames(XP[[i]]) <- list(input.seq$seq.num, paste("h", 1:input.seq$ploidy, sep = ""))
    XQ <- XP[[1]]
    for(i in 1:nrow(XQ)){
      if(dq[i] == ploidy)
        XQ[i, ] <- rep(1, ploidy)
      else
        XQ[i, ] <- rep(0, ploidy)
    }
    return(list(P = XP, Q = list(XQ)))
  } else if (info.par == "p2")
  {
    sq <- lapply(sh, function(x) x$Q)
    XQ <- concatenate_new_marker(d = dq[1], ploidy = input.seq$ploidy)
    for (i in 2:length(dq)) {
      if (is.null(XQ)) {
        warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
        return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
      }
      XQ <- concatenate_new_marker(X = XQ, d = dq, sh = sq, seq.num = input.seq$seq.num, ploidy = input.seq$ploidy, mrk = i)
    }
    if (is.null(XQ)) {
      warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
      return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
    } else for (i in 1:length(XQ)) dimnames(XQ[[i]]) <- list(input.seq$seq.num, paste("h", 1:input.seq$ploidy, sep = ""))
    XP <- XQ[[1]]
    for(i in 1:nrow(XP)){
      if(dp[i] == ploidy)
        XP[i, ] <- rep(1, ploidy)
      else
        XP[i, ] <- rep(0, ploidy)
    }
    return(list(P = list(XP), Q = XQ))
  } else stop("Should not get here.")
}

#' Given a homology group in matrix form, it returns the number shared
#' homologous for all pairs of markers in this group
#'
#' @param M matrix whose columns represent homologous chromosomes and
#'     the rows represent markers
#' @return a vector containing the number of shared homologous for all
#'     pairs of markers
#' @keywords internal
get_ph_conf_ret_sh <- function(M) {
  y <- combn(rownames(M), 2)
  y <- apply(y, 2, sort)
  x <- apply(apply(y, 2, function(x) apply(M[x, ], 2, sum)), 2, function(x) sum(x  ==  2))
  names(x) <- apply(y, 2, paste, collapse = "-")
  x
}

#' Check if all pairwise combinations of elements of \code{input.seq}
#' are contained in \code{twopt}
#'
#' @param input.seq An object of class \code{mappoly.sequence}
#' @param twopt An object of class \code{mappoly.twopt}
#' @return If all pairwise combinations of elements of
#'     \code{input.seq} are contained in \code{twopt}, the function
#'     returns 0. Otherwise, returns the missing pairs.
#' @keywords internal
check_pairwise <- function(input.seq, twopt) {
  if(!(inherits(input.seq, "mappoly.sequence") || is.integer(input.seq) || is.numeric(input.seq) || is.character(input.seq)))
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence', 'numeric' or 'integer'")
  if(!inherits(twopt, "mappoly.twopt"))
    stop(deparse(substitute(twopt)), " is not an object of class 'mappoly.twopt' or 'poly.haplo.est.two.pts.pairwise'")
  id.seq <- input.seq
  if(inherits(input.seq, "mappoly.sequence"))
    id.seq <- input.seq$seq.num
  dpl <- duplicated(id.seq)
  if(any(dpl))
    stop("There are duplicated markers in the sequence:\n Check markers: ", unique(id.seq[dpl]), " at position(s) ", which(dpl))
  index <- apply(apply(combn(id.seq, 2), 2, sort), 2, paste, collapse = "-")
  miss.pairs <- which(is.na(match(index, names(twopt$pairwise))))
  if(length(miss.pairs) > 0)
    return(index[miss.pairs])
  return(0)
}

#' Get the recombination fraction for a sequence of markers given an
#' object of class \code{mappoly.twopt} and a list
#' containing the linkage phase configuration. This list can be found
#' in any object of class \code{two.pts.linkage.phases}, in
#' x$config.to.test$'Conf-i', where x is the object of class
#' \code{two.pts.linkage.phases} and i is one of the possible
#' configurations.
#'
#' @param twopt an object of class \code{mappoly.twopt}
#' @param ph.list a list containing the linkage phase configuration. This
#'     list can be found in any object of class
#'     \code{two.pts.linkage.phases}, in x$config.to.test$'Conf-i',
#'     where x is the object of class \code{two.pts.linkage.phases}
#'     and i is one of the possible configurations.
#' @return a vector with the recombination fraction between markers
#'     present in ph.list, for that specific order.
#' @keywords internal
get_rf_from_list <- function(twopt, ph.list) {
  nm <- as.numeric(names(ph.list$P))
  rf <- numeric(length(nm) - 1)
  for (i in 1:(length(nm) - 1)) {
    id <- paste(sort(c(nm[i], nm[i + 1])), collapse = "-")
    if (all(ph.list$P[[i]]  ==  0) || all(ph.list$P[[i + 1]]  ==  0))
      id.ph.P <- "0" else id.ph.P <- sum(!is.na(match(ph.list$P[[i]], ph.list$P[[i + 1]])))
      if (all(ph.list$Q[[i]]  ==  0) || all(ph.list$Q[[i + 1]]  ==  0))
        id.ph.Q <- "0" else id.ph.Q <- sum(!is.na(match(ph.list$Q[[i]], ph.list$Q[[i + 1]])))
        rf[i] <- twopt$pairwise[[id]][paste(id.ph.P, id.ph.Q, sep = "-"), 2]
  }
  return(rf)
}


#' List of linkage phases
#'
#' Returns a list of possible linkage phase configurations using
#' the two-point information contained in the object \code{mappoly.twopt}
#' as elimination criteria
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#' 
#' @param thres the LOD threshold used to determine whether linkage phases
#'     compared via two-point analysis should be considered
#'     
#' @param twopt an object of class \code{mappoly.twopt}
#'     containing the two-point information
#'     
#' @param mrk.to.add marker to be added to the end of the linkage
#'     group. If \code{NULL} (default) adds all markers contained in
#'     \code{input.seq}. Mostly for internal usage
#'     
#' @param prev.info (optional) an object of class \code{two.pts.linkage.phases}
#'     containing the previous info about linkage phase configuration.
#'     Mostly for internal usage
#'     
#' @param x an object of the class \code{two.pts.linkage.phases}
#' 
#' @param ... currently ignored
#' 
#' @return An object of class \code{two.pts.linkage.phases} which
#'     contains the following structure: 
#'     \item{config.to.test}{a matrix with all possible linkage phase configurations
#'      for both parents, P and Q} 
#'     \item{rec.frac}{a matrix with all recombination fractions}
#'     \item{ploidy}{the ploidy level}
#'     \item{seq.num}{the sequence of markers}
#'     \item{thres}{the LOD threshold}
#'     \item{data.name}{the dataset name}
#'
#' @examples
#' seq.all.mrk <- make_seq_mappoly(hexafake, 'all')
#' id <- get_genomic_order(seq.all.mrk)
#' s <- make_seq_mappoly(id)
#' seq10 <- make_seq_mappoly(hexafake, s$seq.mrk.names[1:10])
#' twopt <- est_pairwise_rf(seq10)
#' 
#' ## Using the first 10 markers 
#' l10.seq.3.0 <- ls_linkage_phases(input.seq = seq10, thres = 3, twopt = twopt)
#' l10.seq.3.0
#' plot(l10.seq.3.0)
#' l10.seq.2.0 <- ls_linkage_phases(input.seq = seq10, thres = 2.0, twopt = twopt)
#' l10.seq.2.0
#' plot(l10.seq.2.0)
#' l10.seq.1.0 <- ls_linkage_phases(input.seq = seq10, thres = 1.0, twopt = twopt)
#' l10.seq.1.0
#' plot(l10.seq.1.0)
#' 
#' ## Using the first 5 markers 
#' seq5 <- make_seq_mappoly(hexafake, s$seq.mrk.names[1:5])
#' l5.seq.5.0 <- ls_linkage_phases(input.seq = seq5, thres = 5, twopt = twopt)
#' l5.seq.5.0
#' plot(l5.seq.5.0)
#' l5.seq.3.0 <- ls_linkage_phases(input.seq = seq5, thres = 3, twopt = twopt)
#' l5.seq.3.0
#' plot(l5.seq.3.0)
#' l5.seq.1.0 <- ls_linkage_phases(input.seq = seq5, thres = 1, twopt = twopt)
#' l5.seq.1.0
#' plot(l5.seq.1.0)
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \doi{10.1534/g3.119.400378} 
#'
#' @importFrom utils tail
#' @export ls_linkage_phases

ls_linkage_phases <- function(input.seq, 
                              thres, 
                              twopt, 
                              mrk.to.add = NULL, 
                              prev.info = NULL) {
  info.par <-detect_info_par(input.seq)
  if (is.null(mrk.to.add) || is.null(prev.info)) 
  {
    X <- elim_conf_using_two_pts(input.seq, twopt, thres)
    mrk.seq.temp <- input.seq$seq.num
    if (!is.list(X))
      return(structure(list(config.to.test = NULL, rec.frac = NULL, 
                            ploidy = input.seq$ploidy, 
                            seq.num = input.seq$seq.num, thres = thres, 
                            data.name = input.seq$data.name),
                       class = "two.pts.linkage.phases"))
  } 
  else 
  {
    seq.temp <- c(input.seq$seq.num, mrk.to.add)
    index <- apply(apply(combn(seq.temp, 2), 2, sort), 2, paste, collapse = "-")
    sh <- lapply(twopt$pairwise[index], get_indices_from_selected_phases, thres = thres)
    ploidy <- input.seq$ploidy
    data.name <- input.seq$data.name
    if(info.par == "both")
    {
      P <- all_phases_with_added_mrk(prev.info, seq.temp, sh, ploidy, 1, data.name)
      Q <- all_phases_with_added_mrk(prev.info, seq.temp, sh, ploidy, 2, data.name)
    }
    else if (info.par == "p1") 
    {
      P <- all_phases_with_added_mrk(prev.info, seq.temp, sh, ploidy, 1, data.name)
      dq <- get(input.seq$data.name, pos = 1)$dosage.p2
      Q <- P
      for(i in 1:length(Q)) 
      {
        dqt <- dq[as.numeric(rownames(Q[[i]]))]
        for(j in 1:nrow(Q[[i]])){
          if(dqt[j] == ploidy)
            Q[[i]][j, ] <- rep(1, ploidy)
          else
            Q[[i]][j, ] <- rep(0, ploidy)
        }
      }
    }
    else if (info.par == "p2")
    {
      Q <- all_phases_with_added_mrk(prev.info, seq.temp, sh, ploidy, 2, data.name)
      dp <- get(input.seq$data.name, pos = 1)$dosage.p1
      P <- Q
      for(i in 1:length(P)) 
      {
        dpt <- dp[as.numeric(rownames(P[[i]]))]
        for(j in 1:nrow(P[[i]])){
          if(dpt[j] == ploidy)
            P[[i]][j, ] <- rep(1, ploidy)
          else
            P[[i]][j, ] <- rep(0, ploidy)
        }
      }
    }
    else stop("Should not get here.")
    X <- list(P = P, Q = Q)
    input.seq <- make_seq_mappoly(get(input.seq$data.name), seq.temp, data.name = input.seq$data.name)
    mrk.seq.temp <- input.seq$seq.num
  }
  ## Insert marker to be added at the end of the linkage group.
  config.to.test <- vector("list", length(X$P) * length(X$Q))
  count <- 1
  for (i in 1:length(X$P)) {
    for (j in 1:length(X$Q)) {
      config.to.test[[count]] <- list(P = ph_matrix_to_list(X$P[[i]]), Q = ph_matrix_to_list(X$Q[[j]]))
      names(config.to.test[[count]]$P) <- names(config.to.test[[count]]$Q) <- rownames(X$P[[i]])
      count <- count + 1
    }
  }
  # Filtering for identical configurations
  id <- which(!duplicated(config.to.test))
  config.to.test <- config.to.test[id]
  structure(list(config.to.test = config.to.test, 
                 rec.frac = NULL, 
                 ploidy = input.seq$ploidy, 
                 seq.num = input.seq$seq.num, 
                 thres = thres, 
                 data.name = input.seq$data.name),
            class = "two.pts.linkage.phases")
}


#' @rdname ls_linkage_phases
#' @keywords internal
#' @export
print.two.pts.linkage.phases <- function(x, ...) {
  cat("\nThis object is too complex to print. Here is a summary:")
  cat("\n---------------------------------------------")
  cat("\nThere is (are) ", length(x$config.to.test), " possible linkage phase(s)\nbased on two-point analysis.")
  cat("\n---------------------------------------------")
  cat("\nThe threshold assumed to discard unlikely\nlinkage phases was ", x$thres, "\n")
}
#' @rdname ls_linkage_phases
#' @keywords internal
#' @export
plot.two.pts.linkage.phases <- function(x, ...) {
  if (length(x$config.to.test)  ==  1) {
    draw_phases(ploidy = get(x$data.name)$ploidy, hom.allele.p = x$config.to.test[1][[1]]$P, hom.allele.q = x$config.to.test[1][[1]]$Q)
  } else {
    n.col <- ceiling(sqrt(length(x$config.to.test)))
    n.row <- ceiling(length(x$config.to.test)/n.col)
    oldpar <- par(xaxt = "n", bty = "n", mar = c(2, 2, 2, 2), mfrow = c(max(n.row, n.col), max(n.row, n.col)))
    on.exit(par(oldpar))
    for (k in names(x$config.to.test)) {
      draw_phases(ploidy = get(x$data.name)$ploidy, hom.allele.p = x$config.to.test[k][[1]]$P, hom.allele.q = x$config.to.test[k][[1]]$Q)
    }
  }
  
}

#' Plot the linkage phase configuration given a list of homologous chromosomes
#'
#' @param ploidy ploidy level
#' @param hom.allele.p a \code{list} of vectors containing linkage
#'     phase configuration for parent P. Each vector contains the
#'     numbers of the homologous chromosomes in which the alleles are
#'     located.
#' @param hom.allele.q same for parent Q
#' @keywords internal
#' @importFrom graphics lines par plot points text
#'
draw_phases <- function(ploidy, hom.allele.p, hom.allele.q) {
  col1 <- "#e41a1c"
    col2 <- "#377eb8"
      n.mrk <- length(hom.allele.p)
      plot(c(0, 22), c(0, -(ploidy + 15)), type = "n", axes = FALSE, xlab = "", main = "", ylab = "")
      for (i in -(1:ploidy)) {
        lines(c(0, 10), c(i, i), lwd = 1, col = "darkgray", lty = 2)
        lines(c(12, 22), c(i, i), lwd = 1, col = "darkgray", lty = 2)
      }
      pos.p <- cumsum(c(0, rep(1, n.mrk - 1)/sum(rep(1, n.mrk - 1))) * 10)
      for (i in 1:n.mrk) {
        points(x = rep(pos.p[i], ploidy), y = -c(1:ploidy), pch = 15, col = col2, cex = 2)
        if (any(hom.allele.p[[i]] != 0))
          points(x = rep(pos.p[i], length(hom.allele.p[[i]])), y = -hom.allele.p[[i]], col = col1, pch = 15, cex = 2)
      }
      pos.q <- pos.p + 12
      for (i in 1:n.mrk) {
        points(x = rep(pos.q[i], ploidy), y = -c(1:ploidy), col = col2, pch = 15, cex = 2)
        if (any(hom.allele.q[[i]] != 0))
          points(x = rep(pos.q[i], length(hom.allele.q[[i]])), y = -hom.allele.q[[i]], col = col1, pch = 15, cex = 2)
      }
      text(x = 11, y = -(ploidy + 1)/2, labels = "X", cex = 1)
}

#' Enumerate all phases with an added marker
#'
#' @param void internal function to be documented
#' @keywords internal
all_phases_with_added_mrk <- function(prev.info, seq.temp, sh, ploidy, parent, data.name){
  L <- lapply(prev.info$config.to.test, function(x) x[[parent]])
  M <- lapply(L, function(x) ph_list_to_matrix(x, ploidy))
  d <- get(data.name)[[c("dosage.p1", "dosage.p2")[parent]]][seq.temp]
  s <- lapply(sh, function(x) x[[parent]])
  X <- concatenate_new_marker(X = M, 
                              d = d, 
                              sh = s, 
                              seq.num = seq.temp, 
                              ploidy = ploidy, 
                              mrk = length(seq.temp))
  for (i in 1:length(X)) 
    rownames(X[[i]]) <- seq.temp
  return(X)
}
