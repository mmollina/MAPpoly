#' Counts for recombinant classes in a polyploid parent.
#'
#' The conditional probability of a genotype at locus \eqn{k+1} given the genotype at locus
#' \eqn{k} is ...
#'
#' @param ploidy Ploidy level
#' @param gen.par.mk1 Genotype of marker mk1 (vector  \eqn{x \in 0, \cdots, m})
#' @param gen.par.mk2 Genotype of marker mk2 (vector  \eqn{x \in 0, \cdots, m})
#' @param gen.prog.mk1 Dose of marker mk1 on progeny
#' @param gen.prog.mk2 Dose of marker mk2 on progeny
#' @return S3 object; a list consisting of
#' \item{counts}{counts for each one of the classes}
#' @keywords internal
#'
#' @useDynLib mappoly
get_counts_one_parent <- function(ploidy, gen.par.mk1, gen.par.mk2, gen.prog.mk1, gen.prog.mk2) {
    res <- .Call("get_counts_one_parent_cpp", as.numeric(ploidy), as.numeric(gen.par.mk1), as.numeric(gen.par.mk2), as.numeric(gen.prog.mk1), as.numeric(gen.prog.mk2), 
        as.numeric(rep(0, ploidy + 1)), PACKAGE = "mappoly")
    return(res[[6]])
}

#' Counts for recombinant classes
#'
#' @param void internal function to be documented
#' @keywords internal
get_counts_two_parents <- function(x = c(2, 2), ploidy, p.k, p.k1, q.k, q.k1, verbose = FALSE, joint.prob = FALSE) {
    gen.prog.mk1 <- x[1]
    gen.prog.mk2 <- x[2]
    if (verbose) {
        cat("Ploidy: ", ploidy, "\n")
        M <- matrix(rep(letters[1:2], 2), ploidy, 4, byrow = TRUE)
        M[1 + p.k, 1] <- "A"
        M[1 + p.k1, 2] <- "B"
        M[1 + q.k, 3] <- "A"
        M[1 + q.k1, 4] <- "B"
        format(apply(M, 1, function(x) cat(c("\t", x[1], "--------", x[2], "               ", x[3], "--------", x[4], "\n"), collapse = "")))
        cat("\n---------------------------------------------------\n\n")
    }
    ## Possible dosages in gametes from P
    dpk <- 0:length(p.k)
    dpk1 <- 0:length(p.k1)
    ## Possible dosages in gametes from Q
    dqk <- 0:length(q.k)
    dqk1 <- 0:length(q.k1)
    ## Combining dosages from P and Q (locus k)
    comb.all.gam.k <- expand.grid(dpk, dqk)
    ## Combining dosages from P and Q (locus k+1)
    comb.all.gam.k1 <- expand.grid(dpk1, dqk1)
    ## Combination of gametes that have x[1] doses in k
    pos.k <- comb.all.gam.k[apply(comb.all.gam.k, 1, sum) == x[1], ]
    ## Combintation of gametes that have x[2] doses in k+1
    pos.k1 <- comb.all.gam.k1[apply(comb.all.gam.k1, 1, sum) == x[2], ]
    r <- NULL
    den <- 0
    for (i in 1:nrow(pos.k)) {
        b <- NULL
        for (j in 1:nrow(pos.k1)) {
            a1 <- get_counts_one_parent(ploidy, p.k, p.k1, pos.k[i, 1], pos.k1[j, 1])
            a2 <- get_counts_one_parent(ploidy, q.k, q.k1, pos.k[i, 2], pos.k1[j, 2])
            r <- rbind(r, kronecker(a1[-(2 + ploidy/2)], a2[-(2 + ploidy/2)]))
            b <- c(b, a1[2 + ploidy/2] * a2[2 + ploidy/2])
        }
        den <- den + mean(b)
    }
    r <- apply(r, 2, sum)
    y <- apply(expand.grid(0:(ploidy/2), 0:(ploidy/2)), 1, function(x) paste(sort(x), collapse = ""))
    names(r) <- y
    res <- NULL
    for (i in sort(unique(y))) res <- c(res, sum(r[names(r) == i]))
    if (!joint.prob) 
        res <- res/den
    names(res) <- sort(unique(y))
    res
}

#' Counts for recombinant classes
#'
#' @param void internal function to be documented
#' @keywords internal
get_counts <- function(m, P.k = NULL, P.k1 = NULL, Q.k = NULL, Q.k1 = NULL, verbose = FALSE, make.names = FALSE, joint.prob = FALSE) {
    if (verbose) {
        cat("Ploidy: ", m, "\n")
        M <- matrix(rep(letters[1:2], 2), m, 4, byrow = TRUE)
        M[1 + P.k, 1] <- "A"
        M[1 + P.k1, 2] <- "B"
        M[1 + Q.k, 3] <- "A"
        M[1 + Q.k1, 4] <- "B"
        format(apply(M, 1, function(x) cat(c("\t", x[1], "--------", x[2], "               ", x[3], "--------", x[4], "\n"), collapse = "")))
        cat("\n---------------------------------------------------\n\n")
    }
    
    if (all(is.null(P.k))) 
        dP.k <- 0 else if (length(P.k) > m/2) 
        dP.k <- (m/2):(m/2 + length(P.k) - m) else dP.k <- 0:length(P.k)
    if (all(is.null(P.k1))) 
        dP.k1 <- 0 else if (length(P.k1) > m/2) 
        dP.k1 <- (m/2):(m/2 + length(P.k1) - m) else dP.k1 <- 0:length(P.k1)
    
    if (all(is.null(Q.k))) 
        dQ.k <- 0 else if (length(Q.k) > m/2) 
        dQ.k <- (m/2):(m/2 + length(Q.k) - m) else dQ.k <- 0:length(Q.k)
    
    if (all(is.null(Q.k1))) 
        dQ.k1 <- 0 else if (length(Q.k1) > m/2) 
        dQ.k1 <- (m/2):(m/2 + length(Q.k1) - m) else dQ.k1 <- 0:length(Q.k1)
    counts <- NULL
    bla <- sort(unique(kronecker(dP.k, dQ.k, "+")))
    ble <- sort(unique(kronecker(dP.k1, dQ.k1, "+")))
    bli <- expand.grid(ble, bla)[, 2:1]
    blo <- bli[1:ceiling(nrow(bli)/2), ]
    if (make.names == TRUE) 
        counts <- matrix(NA, nrow = nrow(bli)) else {
        counts <- t(apply(blo, 1, get_counts_two_parents, ploidy = m, p.k = P.k, p.k1 = P.k1, q.k = Q.k, q.k1 = Q.k1, joint.prob = joint.prob))
        if (nrow(bli) == 1) {
            rownames(counts) <- apply(bli, 1, paste, collapse = " ")
            return(counts)
        }
        if (nrow(bli)%%2 == 1) {
            counts <- rbind(counts, counts[(nrow(counts) - 1):1, ])
        } else {
            counts <- rbind(counts, counts[nrow(counts):1, ])
        }
    }
    rownames(counts) <- apply(bli, 1, paste, collapse = " ")
    return(counts)
}

#' Counts for recombinant classes
#'
#' return the counts of each recombinant class (for two loci) in
#' polyploid cross. The results of this function contains several
#' matrices each one corresponding to one possible linkage phase. The
#' associated names in the matrices indicates the number of shared
#' homologous chromosomes. The row names indicates the dosage in loci
#' k and k+1 respectively
#'
#' @param void internal function to be documented
#' @keywords internal
get_counts_all_phases <- function(x, m, verbose = FALSE, make.names = FALSE, joint.prob = FALSE) {
    pk <- x[1]
    pk1 <- x[2]
    qk <- x[3]
    qk1 <- x[4]
    if (any(is.na(c(m, pk, pk1, qk, qk1)))) 
        return(NULL)
    if (any(c(pk, pk1) == 0)) 
        sh.p <- 0 else {
        sh.p <- min(pk, pk1):0
        if (length(sh.p) > m - max(pk, pk1)) 
            sh.p <- sh.p[1:(m - max(pk, pk1) + 1)]
    }
    if (any(c(qk, qk1) == 0)) 
        sh.q <- 0 else {
        sh.q <- min(qk, qk1):0
        if (length(sh.q) > m - max(qk, qk1)) 
            sh.q <- sh.q[1:(m - max(qk, qk1) + 1)]
    }
    if (pk == 0) 
        pk <- NULL else pk <- 0:(pk - 1)
    if (pk1 == 0) 
        pk1 <- NULL else pk1 <- 0:(pk1 - 1)
    if (qk == 0) 
        qk <- NULL else qk <- 0:(qk - 1)
    if (qk1 == 0) 
        qk1 <- NULL else qk1 <- 0:(qk1 - 1)
    pk.ph <- NULL
    pk1.ph <- NULL
    if (length(pk) < length(pk1)) {
        for (i in 0:(length(sh.p) - 1)) {
            pk.ph <- rbind(pk.ph, pk)
            pk1.ph <- rbind(pk1.ph, pk1 + i)
        }
    } else {
        for (i in 0:(length(sh.p) - 1)) {
            if (!is.null(pk)) 
                pk.ph <- rbind(pk.ph, pk + i)
            pk1.ph <- rbind(pk1.ph, pk1)
        }
    }
    qk.ph <- NULL
    qk1.ph <- NULL
    if (length(qk) < length(qk1)) {
        for (i in 0:(length(sh.q) - 1)) {
            qk.ph <- rbind(qk.ph, qk)
            qk1.ph <- rbind(qk1.ph, qk1 + i)
        }
    } else {
        for (i in 0:(length(sh.q) - 1)) {
            if (!is.null(qk)) 
                qk.ph <- rbind(qk.ph, qk + i)
            qk1.ph <- rbind(qk1.ph, qk1)
        }
    }
    pk.num <- NULL
    if (any(is.null(pk.ph), is.null(pk1.ph))) 
        pk.num <- 0 else {
        for (i in 1:nrow(pk.ph)) pk.num <- c(pk.num, sum(!is.na(match(pk.ph[i, ], pk1.ph[i, ]))))
    }
    qk.num <- NULL
    if (any(is.null(qk.ph), is.null(qk1.ph))) 
        qk.num <- 0 else {
        for (i in 1:nrow(qk.ph)) qk.num <- c(qk.num, sum(!is.na(match(qk.ph[i, ], qk1.ph[i, ]))))
    }
    a.names <- expand.grid(qk.num, pk.num)
    a <- vector("list", length(pk.num) * length(qk.num))
    names(a) <- apply(a.names, 1, function(x) paste(rev(x), collapse = "-"))
    for (i in 1:length(pk.num)) {
        for (j in 1:length(qk.num)) {
            if (verbose) 
                print(names(a)[(i - 1) * length(qk.num) + j])
            a[[(i - 1) * length(qk.num) + j]] <- get_counts(m, pk.ph[i, ], pk1.ph[i, ], qk.ph[j, ], qk1.ph[j, ], verbose = verbose, make.names = make.names, 
                joint.prob = joint.prob)
        }
    }
    a
}
