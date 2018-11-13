#' Extract the LOD Scores in a \code{'mappoly.compare'} or
#' \code{'mappoly.try'} object
#' @param x an object of class \code{'mappoly.compare'} or
#'     \code{'mappoly.try'}
#' @param sorted logical. if \code{TRUE}, the LOD Scores are displayed
#'     in a decreasing order
#' @return a numeric vector containing the LOD Scores
#' @keywords internal
#' @export
get_LOD <- function(x, sorted = TRUE) {
  w <- sapply(x$maps, function(x) x$loglike)
  LOD <- w - max(w)
  if (sorted)
    LOD <- sort(LOD, decreasing = TRUE)
  abs(LOD)
}

#' Get recombination fraction from a matrix
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_rf_from_mat <- function(M){
  r<-numeric(nrow(M)-1)
  for(i in 1:(nrow(M)-1)){
    r[i]<-M[i, i+1]
  }
  r
}

#' Reverse map
#'
#' Provides a reverse version of a given map.
#'
#' @param input.map An object of class \code{mappoly.map}
#' @export
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2017) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_
#'
rev_map<-function(input.map)
{
  output.map<-input.map
  for(i in 1:length(output.map$maps))
  {
    output.map$maps[[i]]$seq.num <- rev(input.map$maps[[i]]$seq.num)
    output.map$maps[[i]]$seq.rf <- rev(input.map$maps[[i]]$seq.rf)
    output.map$maps[[i]]$seq.ph$P <- rev(input.map$maps[[i]]$seq.ph$P)
    output.map$maps[[i]]$seq.ph$Q <- rev(input.map$maps[[i]]$seq.ph$Q)
  }
  return(output.map)
}

#' Selects the class with high probability and return it instead of
#' the distribution of probabilities for all classes
#'
#' @param geno the genotypes contained in the object
#'     \code{'mappoly.data'}
#' @param prob.thres probability threshold to accept the genotyupe as
#'     the correct one. Values below this genotype are assumed as
#'     missing data
#' @return a matrix containing the doses of each genotype and
#'     marker. Missing data are represented by NAs
#' @keywords internal
#' @export
dist_prob_to_class <- function(geno, prob.thres = 0.95) {
  d <- apply(geno[, -c(1:2)], 1, function(x) {
    a <- which(x > prob.thres)
    if (length(a) > 0)
      return(a - 1)
    return(NA)
  })
  d <- cbind(geno[, c(1:2)], d)
  mrk <- unique(d$mrk)
  ind <- unique(d$ind)
  z <- matrix(NA, length(mrk), length(ind))
  dimnames(z) <- list(mrk, ind)
  for (i in mrk){
    s.temp<-subset(d, mrk == i)
    z[i, ] <- as.matrix(s.temp[match(ind, s.temp$ind),]$d)
  }
  z
}


#' Map functions
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
mf_k <- function(d) 0.5 * tanh(d/50)
#'
#' Map functions
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
mf_h <- function(d) 0.5 * (1 - exp(-d/50))
#' Map functions
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
mf_m <- function(d) sapply(d, function(a) min(a/100, 0.5))
#' Map functions
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
imf_k <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  50 * atanh(2 * r)
}
#' Map functions
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}
#' Map functions
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
imf_m <- function(r) sapply(r, function(a) min(a * 100, 50))


#' Compare two polyploid haplotypes stored in list format
#'
#' @param m ploidy level
#' @param h1 homology group 1
#' @param h2 homology group 2
#' @return a numeric vector of size \code{m} indicating which
#'     homologous in h2 represents the homologous in h1. If there is
#'     no correspondence, i.e. different homologous, it returns NA for
#'     that homologous.
#' @keywords internal
#' @export compare_haplotypes
compare_haplotypes <- function(m, h1, h2) {
  I1 <- matrix(0, m, length(h1))
  I2 <- matrix(0, m, length(h2))
  for (i in 1:length(h1)) {
    I1[h1[[i]], i] <- 1
    I2[h2[[i]], i] <- 1
  }
  a <- apply(I1, 1, paste, collapse = "")
  b <- apply(I2, 1, paste, collapse = "")
  haplo.ord <- match(a, b)
  list(is.same.haplo = !any(is.na(haplo.ord)), haplo.ord = haplo.ord)
}

#' Plot two overlapped haplotypes
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export plot_compare_haplotypes
plot_compare_haplotypes <- function(m, hom.allele.p1, hom.allele.q1, hom.allele.p2 = NULL, hom.allele.q2 = NULL) {
  nmmrk<-names(hom.allele.p1)
  op <- par(mar = c(5.1, 4.1, 4.1, 2.1))
  o1 <- order(apply(ph_list_to_matrix(hom.allele.p1, m), 2, paste, collapse = ""), decreasing = TRUE)
  hom.allele.p1 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.p1, m)[, o1])
  o2 <- order(apply(ph_list_to_matrix(hom.allele.q1, m), 2, paste, collapse = ""), decreasing = TRUE)
  hom.allele.q1 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.q1, m)[, o2])

  if (!is.null(hom.allele.p2)) {
    o3 <- order(apply(ph_list_to_matrix(hom.allele.p2, m), 2, paste, collapse = ""), decreasing = TRUE)
    hom.allele.p2 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.p2, m)[, o3])
  }
  if (!is.null(hom.allele.q2)) {
    o4 <- order(apply(ph_list_to_matrix(hom.allele.q2, m), 2, paste, collapse = ""), decreasing = TRUE)
    hom.allele.q2 <- ph_matrix_to_list(ph_list_to_matrix(hom.allele.q2, m)[, o4])
  }
  col2 <- 0  #rgb(1,0,.5,0.5)
  col1 <- rgb(1, 0, 0, 0.5)
  col3 <- rgb(0, 0, 1, 0.5)
  n.mrk <- length(hom.allele.p1)
  plot(c(0, 22), c(0, -(m + 15)), type = "n", axes = FALSE, xlab = "", main = "", ylab = "")
  for (i in -(1:m)) {
    lines(c(0, 10), c(i, i), lwd = 1, col = "darkgray", lty = 2)
    lines(c(12, 22), c(i, i), lwd = 1, col = "darkgray", lty = 2)
  }
  pos.p <- cumsum(c(0, rep(1, n.mrk - 1)/sum(rep(1, n.mrk - 1))) * 10)
  for (i in 1:n.mrk) {
    points(x = rep(pos.p[i], m), y = -c(1:m), pch = 15, col = col2, cex = 2)
    if (any(hom.allele.p1[[i]] != 0))
      points(x = rep(pos.p[i], length(hom.allele.p1[[i]])), y = -hom.allele.p1[[i]], col = col1, pch = 15, cex = 2)
    if (any(hom.allele.p2[[i]] != 0))
      if (!is.null(hom.allele.p2))
        points(x = rep(pos.p[i], length(hom.allele.p2[[i]])), y = -hom.allele.p2[[i]], col = col3, pch = 15, cex = 2)
  }
  text(x = pos.p, y = rep(-(m+1), length(pos.p)), labels = nmmrk, cex = .5)
  pos.q <- pos.p + 12
  for (i in 1:n.mrk) {
    points(x = rep(pos.q[i], m), y = -c(1:m), col = col2, pch = 15, cex = 2)
    if (any(hom.allele.q1[[i]] != 0))
      points(x = rep(pos.q[i], length(hom.allele.q1[[i]])), y = -hom.allele.q1[[i]], col = col1, pch = 15, cex = 2)
    if (any(hom.allele.q2[[i]] != 0))
      if (!is.null(hom.allele.q2))
        points(x = rep(pos.q[i], length(hom.allele.q2[[i]])), y = -hom.allele.q2[[i]], col = col3, pch = 15, cex = 2)
  }
  text(x = pos.q, y = rep(-(m+1) , length(pos.q)), labels = nmmrk, cex = .5)
  text(x = 11, y = -(m + 1)/2, labels = "X", cex = 1)
  par(op)
}


#' Summary of a set of markers
#'
#' @param input.data an object \code{'mappoly.data'}
#' @param mrks marker sequence index
#' @export
print_mrk<-function(input.data, mrks)
{
  for(i in 1:length(mrks))
  {
    x<-input.data$geno.dose[mrks[i], ]
    x[x==input.data$m+1]<-NA
    cat("\n", input.data$mrk.names[mrks[i]])
    cat("\n----------------------------------")
    cat("\n dosage P1: ", input.data$dosage.p[mrks[i]])
    cat("\n dosage P2: ", input.data$dosage.q[mrks[i]])
    cat("\n----")
    cat("\n dosage distribution\n")
    z<-table(as.numeric(x), useNA = "always")
    names(z)[is.na(names(z))]<-"mis"
    print(z)
    cat("----")
    cat("\n expected polysomic segregation\n")
    y<-segreg_poly(input.data$m, dP = input.data$dosage.p[mrks[i]], input.data$dosage.q[mrks[i]])
    names(y)<-0:input.data$m
    print(y)
    cat("----------------------------------\n")
  }
}

#' Check if it is possible to estimate the recombination
#' fraction between neighbour markers using two-point
#' estimation
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
pos_twopt_est<-function(input.seq)
{
  dp<-abs(abs(input.seq$seq.dose.p-(input.seq$m/2))-(input.seq$m/2))
  dq<-abs(abs(input.seq$seq.dose.q-(input.seq$m/2))-(input.seq$m/2))
  y<-numeric(length(input.seq$seq.num)-1)
  for(i in 1:length(y))
    y[i]<-(dp[i] == 0 && dq[i+1] == 0) ||
    (dp[i+1] == 0 && dq[i] == 0) ||
    (dp[i] == 0 && dq[i] == 0) ||
    (dp[i+1] == 0 && dq[i+1] == 0)
  y<-as.logical(y)
  y
}


#' N! combination
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
perm_tot <- function(v) {
    n <- length(v)
    result <- v
    if (n > 1) {
        M <- perm_tot(v[2:n])
        result <- cbind(v[1], M)
        if (n > 2) {
            for (i in 2:(n - 1)) {
                N <- cbind(M[, 1:(i - 1)], v[1], M[, i:(n - 1)])
                result <- rbind(result, N)
            }
        }
        N <- cbind(M, v[1])
        result <- rbind(result, N)
    }
    return(result)
}

#' N!/2 combination
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
perm_pars <- function(v) {
    n <- length(v)
    result <- v
    if (n > 2) {
        Mt <- perm_tot(v[2:n])
        result <- cbind(v[1], Mt)
        f <- floor(n/2)
        c <- ceiling(n/2)
        if (n > 3) {
            for (i in 2:f) {
                N <- cbind(Mt[, 1:(i - 1)], v[1], Mt[, i:(n - 1)])
                result <- rbind(result, N)
            }
        }
        if (c > f) {
            Ms <- perm_pars(v[2:n])
            if (n > 3) {
                N <- cbind(Ms[, 1:f], v[1], Ms[, c:(n - 1)])
            } else {
                N <- cbind(Ms[1:f], v[1], Ms[c:(n - 1)])
            }
            result <- rbind(result, N)
        }
    }
    return(result)
}
