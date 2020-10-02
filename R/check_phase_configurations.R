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
    w <- lapply(split(M, seq(NROW(M))), function(x, M) which(x == 1))
    w[sapply(w, function(x) length(x) == 0)] <- 0
    w
}

#' Linkage phase format conversion: list to matrix
#' 
#' This function converts linkage phase configurations from list
#' to matrix form
#'
#' @param L a list of configuration phases
#' 
#' @param m ploidy level
#' 
#' @return a matrix whose columns represent homologous chromosomes and
#'     the rows represent markers
#' 
#' @keywords internal
#' @export
ph_list_to_matrix <- function(L, m) {
    M <- matrix(0, nrow = length(L), ncol = m)
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
    list(P = sort(unique(y[, 1])), Q = sort(unique(y[, 2])))
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
#' @param m the ploidy level
#' @param k1 marker already present on the sequence
#' @param k2 inserted marker
#' @return a unique list of matrices representing linkage phases
#' @keywords internal
generate_all_link_phase_elim_equivalent <- function(X, d, sh, m, k1, k2) {
    mat <- matrix(0, nrow = m, ncol = choose(m, d))
    i <- combn(m, d)
    j <- cumsum(c(0, rep(m, choose(m, d) - 1)))
    mat[as.numeric(apply(i, 1, function(x, y) x + y, y = j))] <- 1
    ct <- NULL
    Y <- NULL
    for (i in 1:ncol(mat)) {
        Y[[i]] <- rbind(X, mat[, i])
        so <- sum(apply(Y[[i]][c(k1, k2), ], 2, sum) == 2)
        if (any(so == sh))
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
#' @param m the ploidy level
#' @param mrk the marker to be inserted
#' @return a unique list of matrices representing linkage phases
#'
#' @keywords internal
concatenate_new_marker <- function(X = NULL, d, sh = NULL, seq.num = NULL, m, mrk = 1) {
    if (is.null(X) & is.null(sh) & mrk == 1 & is.null(seq.num)) {
        Y <- numeric(m)
        if (d != 0)
            Y[1:d] <- 1
        return(list(matrix(Y, ncol = m)))
    }
    Y.final <- NULL
    for (i in 1:length(X)) {
        for (j in (mrk - 1):1) {
            id <- paste(sort(c(seq.num[j], seq.num[mrk])), collapse = "-")
            Y.temp <- generate_all_link_phase_elim_equivalent(X[[i]], d = d[mrk], sh = sh[[id]], m = m, k1 = j, k2 = mrk)
            if (length(Y.temp) == 1 && j == (mrk - 1)) {
                Y <- Y.temp
                best.conf <- 1
                (break)()
            }
            if (j == (mrk - 1)) {
                Y <- Y.temp
                best.conf <- 1:length(Y)
            } else {
                best.conf <- na.omit(match(Y, Y.temp))
                if (length(best.conf) != 0)
                  Y <- Y.temp[best.conf]
                if (length(Y) == 1)
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
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
#' @param thres threshold from which the linkage phases can be
#'     discarded (if abs(ph_LOD) > thres)
#' @return a unique list of matrices representing linkage phases
#' @keywords internal
elim_conf_using_two_pts <- function(input.seq, twopt, thres) {
    if (!inherits(input.seq, "mappoly.sequence"))
        stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
    dp <- input.seq$seq.dose.p
    dq <- input.seq$seq.dose.q
    check <- check_pairwise(input.seq, twopt)
    if (any(check != 0)){
        ## cat("There is no information for pairs: \n")
        ## print(check)
        stop("There is no information for pairs: \n", paste(capture.output(print(check)), collapse = "\n"))
    }
    index <- apply(apply(combn(input.seq$seq.num, 2), 2, sort), 2, paste, collapse = "-")
    w <- twopt$pairwise[index]
    sh <- lapply(w, get_indices_from_selected_phases, thres = thres)
    sp <- lapply(sh, function(x) x$P)
    sq <- lapply(sh, function(x) x$Q)
    XP <- concatenate_new_marker(d = dp[1], m = input.seq$m)
    for (i in 2:length(dp)) {
        if (is.null(XP)) {
            warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
            return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
        }
        XP <- concatenate_new_marker(X = XP, d = dp, sh = sp, seq.num = input.seq$seq.num, m = input.seq$m, mrk = i)
    }
    if (is.null(XP)) {
        warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
        return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
    } else for (i in 1:length(XP)) dimnames(XP[[i]]) <- list(input.seq$seq.num, paste("h", 1:input.seq$m, sep = ""))
    XQ <- concatenate_new_marker(d = dq[1], m = input.seq$m)
    for (i in 2:length(dq)) {
        if (is.null(XP)) {
            warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
            return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
        }
        XQ <- concatenate_new_marker(X = XQ, d = dq, sh = sq, seq.num = input.seq$seq.num, m = input.seq$m, mrk = i)
    }
    if (is.null(XQ)) {
        warning("There are too many conflicts for this threshold level.\nTry anotrher one.")
        return(-1)  #stop('There are too many conflicts for this threshold level.\nTry anotrher one.')
    } else for (i in 1:length(XQ)) dimnames(XQ[[i]]) <- list(input.seq$seq.num, paste("h", 1:input.seq$m, sep = ""))
    list(P = XP, Q = XQ)
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
    x <- apply(apply(y, 2, function(x) apply(M[x, ], 2, sum)), 2, function(x) sum(x == 2))
    names(x) <- apply(y, 2, paste, collapse = "-")
    x
}

#' Check if all pairwise combinations of elements of \code{input.seq}
#' are contained in \code{twopt}
#'
#' @param input.seq An object of class \code{mappoly.sequence}
#' @param twopt An object of class \code{poly.est.two.pts.pairwise}
#' @return If all pairwise combinations of elements of
#'     \code{input.seq} are contained in \code{twopt}, the function
#'     returns 0. Otherwise, returns the missing pairs.
#' @keywords internal
check_pairwise <- function(input.seq, twopt) {
    if(!(inherits(input.seq, "mappoly.sequence") || is.integer(input.seq) || is.numeric(input.seq) || is.character(input.seq)))
        stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence', 'numeric' or 'integer'")
    if(!inherits(twopt, "poly.est.two.pts.pairwise"))
        stop(deparse(substitute(twopt)), " is not an object of class 'poly.est.two.pts.pairwise' or 'poly.haplo.est.two.pts.pairwise'")
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
#' object of class \code{poly.est.two.pts.pairwise} and a list
#' containing the linkage phase configuration. This list can be found
#' in any object of class \code{two.pts.linkage.phases}, in
#' x$config.to.test$'Conf-i', where x is the object of class
#' \code{two.pts.linkage.phases} and i is one of the possible
#' configurations.
#'
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
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
        if (all(ph.list$P[[i]] == 0) || all(ph.list$P[[i + 1]] == 0))
            id.ph.P <- "0" else id.ph.P <- sum(!is.na(match(ph.list$P[[i]], ph.list$P[[i + 1]])))
        if (all(ph.list$Q[[i]] == 0) || all(ph.list$Q[[i + 1]] == 0))
            id.ph.Q <- "0" else id.ph.Q <- sum(!is.na(match(ph.list$Q[[i]], ph.list$Q[[i + 1]])))
        rf[i] <- twopt$pairwise[[id]][paste(id.ph.P, id.ph.Q, sep = "-"), 2]
    }
    return(rf)
}


#' List of linkage phases
#'
#' Returns a list of possible linkage phase configurations using
#' the two-point information contained in the object \code{poly.est.two.pts.pairwise}
#' as elimination criteria
#'
#' @param input.seq an object of class \code{mappoly.sequence}
#' 
#' @param thres the LOD threshold used to determine whether linkage phases
#'     compared via two-point analysis should be considered
#'     
#' @param twopt an object of class \code{poly.est.two.pts.pairwise}
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
#'     \item{m}{the ploidy level}
#'     \item{seq.num}{the sequence of markers}
#'     \item{thres}{the LOD threshold}
#'     \item{data.name}{the dataset name}
#'
#' @examples
#' seq.all.mrk <- make_seq_mappoly(hexafake, 'all')
#' id <- get_genomic_order(seq.all.mrk)
#' seq10 <- make_seq_mappoly(hexafake, rownames(id)[1:10])
#' twopt<-est_pairwise_rf(seq10)
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
#' seq5 <- make_seq_mappoly(hexafake, rownames(id)[1:5])
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
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @importFrom utils tail
#' @export ls_linkage_phases

ls_linkage_phases <- function(input.seq, thres, twopt, mrk.to.add = NULL, prev.info = NULL) {
    input_classes <- c("mappoly.sequence")
    if (!inherits(input.seq, input_classes)) {
        stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
    }
    if (!is.null(mrk.to.add) && length(intersect(mrk.to.add, twopt$seq.num)) == 0)
        stop("Marker ", deparse(substitute(mrk.to.add)), " not found in the twopt object")
    if (is.null(mrk.to.add) || is.null(prev.info)) {
        X <- elim_conf_using_two_pts(input.seq, twopt, thres)
        mrk.seq.temp <- input.seq$seq.num
        if (!is.list(X))
            return(structure(list(config.to.test = NULL, rec.frac = NULL, m = input.seq$m, seq.num = input.seq$seq.num, thres = thres, data.name = input.seq$data.name),
                class = "two.pts.linkage.phases"))

    } else {
        LP <- lapply(prev.info$config.to.test, function(x) x$P)
        MP <- lapply(LP, function(x) ph_list_to_matrix(x, input.seq$m))
        LQ <- lapply(prev.info$config.to.test, function(x) x$Q)
        MQ <- lapply(LQ, function(x) ph_list_to_matrix(x, input.seq$m))
        seq.temp <- c(input.seq$seq.num, mrk.to.add)
        dp <- get(input.seq$data.name)$dosage.p[seq.temp]
        dq <- get(input.seq$data.name)$dosage.q[seq.temp]
        check <- check_pairwise(input.seq, twopt)
        if (any(check != 0)){
            ## cat("There is no information for pairs: \n")
            ## print(check)
            stop("There is no information for pairs: \n", paste(capture.output(print(check)), collapse = "\n"))
        }
        index <- apply(apply(combn(seq.temp, 2), 2, sort), 2, paste, collapse = "-")
        w <- twopt$pairwise[index]
        sh <- lapply(w, get_indices_from_selected_phases, thres = thres)
        sp <- lapply(sh, function(x) x$P)
        sq <- lapply(sh, function(x) x$Q)
        XP <- concatenate_new_marker(X = MP, d = dp, sh = sp, seq.num = seq.temp, m = input.seq$m, mrk = length(seq.temp))
        XQ <- concatenate_new_marker(X = MQ, d = dq, sh = sq, seq.num = seq.temp, m = input.seq$m, mrk = length(seq.temp))
        for (i in 1:length(XP)) dimnames(XP[[i]]) <- list(seq.temp, paste("h", 1:input.seq$m, sep = ""))
        for (i in 1:length(XQ)) dimnames(XQ[[i]]) <- list(seq.temp, paste("h", 1:input.seq$m, sep = ""))
        X <- list(P = XP, Q = XQ)
        input.seq <- make_seq_mappoly(get(input.seq$data.name), seq.temp, data.name = input.seq$data.name)
        mrk.seq.temp <- input.seq$seq.num
    }
    ## Insert marker to be added at the end of the linkage group.
    cT <- lapply(X$P, get_ph_conf_ret_sh)
    cP <- NULL
    for (i in 1:length(cT)) cP <- cbind(cP, cT[[i]])
    nC <- rownames(cP)
    rownames(cP) <- apply(apply(sapply(nC, get_ij), 2, sort), 2, paste, collapse = "-")
    cT <- lapply(X$Q, get_ph_conf_ret_sh)
    cQ <- NULL
    for (i in 1:length(cT)) cQ <- cbind(cQ, cT[[i]])
    nC <- rownames(cQ)
    rownames(cQ) <- apply(apply(sapply(nC, get_ij), 2, sort), 2, paste, collapse = "-")
    check <- check_pairwise(input.seq, twopt)
    if (any(check != 0)){
        ## cat("There is no information for pairs: \n")
        ## print(check)
        stop("There is no information for pairs: \n", paste(capture.output(print(check)), collapse = "\n"))
    }
    y <- combn(input.seq$seq.num, 2)
    y <- apply(y, 2, sort)
    index <- apply(apply(y, 2, sort), 2, paste, collapse = "-")
    sel.rf <- twopt$pairwise[index]
    ph.two.pts <- lapply(sel.rf, function(x) rownames(x))
    config.to.test <- NULL
    rf.vec <- NULL
    count <- 0
    for (i in 1:length(X$P)) {
        for (j in 1:length(X$Q)) {
            check.conf <- paste(cP[, i], cQ[, j], sep = "-")
            names(check.conf) <- rownames(cP)
            pos.conf <- ph.two.pts[names(check.conf)]
            res <- NULL
            if (is.null(mrk.to.add) || is.null(prev.info)) {
                for (k in 1:length(pos.conf)) if (any(check.conf[k] == pos.conf[[k]])) {
                  res <- c(res, TRUE)
                } else res <- c(res, FALSE)
                if (sum(res) == length(pos.conf)) {
                  count <- count + 1
                  config.to.test[[count]] <- list(P = ph_matrix_to_list(X$P[[i]]), Q = ph_matrix_to_list(X$Q[[j]]))
                  names(config.to.test[[count]]$P) <- names(config.to.test[[count]]$Q) <- mrk.seq.temp
                  rf.vec <- rbind(rf.vec, get_rf_from_list(twopt = twopt, ph.list = config.to.test[[count]]))
                }
            } else {
                if (any(check.conf[length(pos.conf)] == pos.conf[[length(pos.conf)]])) {
                  count <- count + 1
                  config.to.test[[count]] <- list(P = ph_matrix_to_list(X$P[[i]]), Q = ph_matrix_to_list(X$Q[[j]]))
                  names(config.to.test[[count]]$P) <- names(config.to.test[[count]]$Q) <- mrk.seq.temp
                  rf.vec <- rbind(rf.vec, get_rf_from_list(twopt = twopt, ph.list = config.to.test[[count]]))
                } else stop("Impossible insert marker at this threshold level.")
            }
        }
    }
    # Filtering for identical configurations
    id <- which(!duplicated(config.to.test))
    config.to.test <- config.to.test[id]
    rf.vec <- matrix(rf.vec[id, ], nrow = length(id))
    rownames(rf.vec) <- names(config.to.test) <- paste("Conf", 1:length(config.to.test), sep = "-")
    structure(list(config.to.test = config.to.test, rec.frac = rf.vec, m = input.seq$m, seq.num = input.seq$seq.num, thres = thres, data.name = input.seq$data.name),
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
    if (length(x$config.to.test) == 1) {
        draw_phases(m = get(x$data.name)$m, hom.allele.p = x$config.to.test[1][[1]]$P, hom.allele.q = x$config.to.test[1][[1]]$Q)
    } else {
        n.col <- ceiling(sqrt(length(x$config.to.test)))
        n.row <- ceiling(length(x$config.to.test)/n.col)
        oldpar <- par(xaxt = "n", bty = "n", mar = c(2, 2, 2, 2), mfrow = c(max(n.row, n.col), max(n.row, n.col)))
        on.exit(par(oldpar))
        for (k in names(x$config.to.test)) {
            draw_phases(m = get(x$data.name)$m, hom.allele.p = x$config.to.test[k][[1]]$P, hom.allele.q = x$config.to.test[k][[1]]$Q)
        }
    }

}

#' Plot the linkage phase configuration given a list of homologous chromosomes
#'
#' @param m ploidy level
#' @param hom.allele.p a \code{list} of vectors containing linkage
#'     phase configuration for parent P. Each vector contains the
#'     numbers of the homologous chromosomes in which the alleles are
#'     located.
#' @param hom.allele.q same for parent Q
#' @keywords internal
#' @importFrom graphics lines par plot points text
#'
draw_phases <- function(m, hom.allele.p, hom.allele.q) {
    col1 <- "#e41a1c"
    col2 <- "#377eb8"
    n.mrk <- length(hom.allele.p)
    plot(c(0, 22), c(0, -(m + 15)), type = "n", axes = FALSE, xlab = "", main = "", ylab = "")
    for (i in -(1:m)) {
        lines(c(0, 10), c(i, i), lwd = 1, col = "darkgray", lty = 2)
        lines(c(12, 22), c(i, i), lwd = 1, col = "darkgray", lty = 2)
    }
    pos.p <- cumsum(c(0, rep(1, n.mrk - 1)/sum(rep(1, n.mrk - 1))) * 10)
    for (i in 1:n.mrk) {
        points(x = rep(pos.p[i], m), y = -c(1:m), pch = 15, col = col2, cex = 2)
        if (any(hom.allele.p[[i]] != 0))
            points(x = rep(pos.p[i], length(hom.allele.p[[i]])), y = -hom.allele.p[[i]], col = col1, pch = 15, cex = 2)
    }
    pos.q <- pos.p + 12
    for (i in 1:n.mrk) {
        points(x = rep(pos.q[i], m), y = -c(1:m), col = col2, pch = 15, cex = 2)
        if (any(hom.allele.q[[i]] != 0))
            points(x = rep(pos.q[i], length(hom.allele.q[[i]])), y = -hom.allele.q[[i]], col = col1, pch = 15, cex = 2)
    }
    text(x = 11, y = -(m + 1)/2, labels = "X", cex = 1)
}
