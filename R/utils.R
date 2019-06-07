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
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
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
#' @importFrom magrittr "%>%"
#' @importFrom reshape melt cast
#' @importFrom dplyr group_by filter arrange
dist_prob_to_class <- function(geno, prob.thres = 0.95) {
  a<-reshape::melt(geno, id.vars = c("mrk", "ind"))
  mrk <- ind <- value <- variable <- NULL # Setting the variables to NULL first
  a$variable<-as.numeric(levels(a$variable))[a$variable]
  b<-a %>%
    dplyr::group_by(mrk, ind) %>%
    dplyr::filter(value > prob.thres) %>%
    dplyr::arrange(mrk, ind, variable)
  z<-reshape::cast(data = b[,1:3], formula = mrk ~ ind, value = "variable")
  rownames(z)<-z[,"mrk"]
  z<-data.matrix(frame = z[,-1])
  n<-setdiff(unique(geno$mrk), rownames(z))
  if(length(n) > 0)
  {
    m<-matrix(NA, nrow = length(n), ncol = ncol(z), dimnames = list(n, colnames(z)))
    z<-rbind(z,m)
  }
  return(z[unique(geno$mrk),])
}

#' Msg function
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
#' @importFrom cli rule
msg <- function(text, line = 1)
    cli::rule(line = line, right = text) %>%
  text_col() %>%
  message()


#' @importFrom rstudioapi isAvailable hasFun getThemeInfo
#' @importFrom crayon white black

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
  
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

#' Color pallete ggplot-like
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
#' @importFrom grDevices hcl
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Update missing information
#'
#' Updates the missing data in the dosage matrix of an object of class 
#' \code{mappoly.data} given a new probability threshold
#' 
#' @param input.data an object of class \code{mappoly.data}
#' 
#' @param prob.thres probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than 'prob.thres' 
#'     are considered as missing data for the dosage calling purposes
#'     
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'     
#' @export
update_missing<-function(input.data, 
                         prob.thres = 0.95)
{
  geno.dose <- dist_prob_to_class(geno = input.data$geno, 
                                  prob.thres = prob.thres)
  geno.dose[is.na(geno.dose)] <- input.data$m + 1
  input.data$geno.dose<-geno.dose
  input.data$prob.thres<-prob.thres
  return(input.data)
}

#' Filter non-conforming classes in F1 segregation
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
filter_non_conforming_classes<-function(input.data, 
                                        prob.thres = NULL)
  {
  m<-input.data$m
  dp<-input.data$dosage.p
  dq<-input.data$dosage.q
  Ds<-array(NA, dim = c(m+1, m+1, m+1))
  for(i in 0:m)
    for(j in 0:m)
      Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
  Dpop<-cbind(dp,dq)
  M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
  M[M!=0]<-1
  ##if no prior probabilities
  if(nrow(input.data$geno)==input.data$n.mrk){
    for(i in 1:nrow(M)){
      id0<-!as.numeric(input.data$geno.dose[i,])%in%(which(M[i,]==1)-1)
      if(any(id0))
        input.data$geno.dose[i,id0]<-(m+1)     
    }
    input.data$prob.thres<-NULL
  return(input.data)
  }
  ## 1 represents conforming classes/ 0 represents non-conforming classes
  dp<-rep(dp, input.data$n.ind)
  dq<-rep(dq, input.data$n.ind)
  M<-do.call("rbind", replicate(input.data$n.ind, M, simplify = FALSE))
  R<-input.data$geno[,-c(1:2)] - input.data$geno[,-c(1:2)]*M
  id1<-apply(R, 1, sum) > 0.3 # if the sum of the excluded classes is greater than 0.3, use segreg_poly
  N<-NULL
  for(i in which(id1))
    N<-rbind(N, Ds[dp[i]+1, dq[i]+1, ])
  input.data$geno[id1,-c(1:2)]<-N
  # if the sum of the excluded classes is greater than zero
  # and smaller than 0.3, assign zero to those classes and normalize the vector
  input.data$geno[,-c(1:2)][R > 0]<-0
  input.data$geno[,-c(1:2)]<-sweep(input.data$geno[,-c(1:2)], 1, rowSums(input.data$geno[,-c(1:2)]), FUN="/")
  if(is.null(prob.thres))
    prob.thres<-input.data$prob.thres
  geno.dose <- dist_prob_to_class(input.data$geno, prob.thres)
  geno.dose[is.na(geno.dose)] <- m + 1
  input.data$geno.dose<-geno.dose
  input.data
}

#' Filter missing genotypes
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing<-function(input.data, filter.thres = 0.2, inter = TRUE)
{
  filter.thres<-1-filter.thres
  op<-par(bg = "gray", xpd = TRUE)
  ANSWER <- "flag"
  mrk <- NULL
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
    {
      na.num<-apply(input.data$geno.dose, 1, function(x,m) sum(x!=m+1), m = input.data$m)
      perc.na<-na.num/input.data$n.ind
      plot(sort(perc.na), xlab = "markers", ylab = "frequency of genotyped markers", axes = FALSE);
      axis(1);axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(filter.thres,2), col = 4, lty = 2)
      text(x = input.data$n.mrk/2, y = filter.thres + 0.05, labels = paste0("Included mrks: ", sum(perc.na >= filter.thres)), adj = 0, col = "darkgreen")
      text(x = input.data$n.mrk/2, y = filter.thres - 0.05, labels = paste0("Excluded mrks: ", sum(perc.na < filter.thres)), adj = 0, col = "darkred")
      ANSWER <- readline("Enter 'y' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
      filter.thres  <- as.numeric(ANSWER)
    }
    rm.mrks.id<-which(perc.na < filter.thres)
      if(length(rm.mrks.id)==0) return(input.data)
    rm.mrks<-names(rm.mrks.id)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      filter(!mrk%in%rm.mrks)
    input.data$geno.dose<-input.data$geno.dose[-rm.mrks.id,]
    input.data$n.mrk <- nrow(input.data$geno.dose)
    input.data$mrk.names <- input.data$mrk.names[-rm.mrks.id]
    input.data$dosage.p <- input.data$dosage.p[-rm.mrks.id]
    input.data$dosage.q <- input.data$dosage.q[-rm.mrks.id]
    input.data$sequence <- input.data$sequence[-rm.mrks.id]
    if(!is.null(input.data$chisq.pval)) 
      input.data$chisq.pval <- input.data$chisq.pval[-rm.mrks.id]
    input.data$sequence.pos <- input.data$sequence.pos[-rm.mrks.id]
    par(op)
    return(input.data)
  } else {
    na.num<-apply(input.data$geno.dose, 1, function(x,m) sum(x!=m+1), m = input.data$m)
    perc.na<-na.num/input.data$n.ind
    rm.mrks.id<-which(perc.na < filter.thres)
    if(length(rm.mrks.id)==0) return(input.data)
    rm.mrks<-names(rm.mrks.id)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      filter(!mrk%in%rm.mrks)
    input.data$geno.dose<-input.data$geno.dose[-rm.mrks.id,]
    input.data$n.mrk <- nrow(input.data$geno.dose)
    input.data$mrk.names <- input.data$mrk.names[-rm.mrks.id]
    input.data$dosage.p <- input.data$dosage.p[-rm.mrks.id]
    input.data$dosage.q <- input.data$dosage.q[-rm.mrks.id]
    input.data$sequence <- input.data$sequence[-rm.mrks.id]
    if(!is.null(input.data$chisq.pval)) 
      input.data$chisq.pval <- input.data$chisq.pval[-rm.mrks.id]
    input.data$sequence.pos <- input.data$sequence.pos[-rm.mrks.id]
    par(op)
    return(input.data)
  }
}


#' Chi-square test
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
mrk_chisq_test<-function(x, m){
  y<-x[-c(1:(m+1))]
  y[y==m+1]<-NA
  y<-table(y, useNA = "always")
  names(y)<-c(names(y)[-length(y)], "NA") 
  seg.exp <- x[0:(m+1)]
  seg.exp <- seg.exp[seg.exp!=0]
  seg.obs <- seg.exp
  seg.obs[names(y)[-length(y)]]<-y[-length(y)]
  pval <- suppressWarnings(stats::chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value)
  pval
}


#' marker filter based on chi-square test
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @importFrom graphics axis
#' @export
filter_segregation<-function(input.data, chisq.pval.thres = 10e-5, inter = TRUE){
  ANSWER <- "flag"
  op<-par(bg = "gray", xpd = TRUE)
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
    {
      plot(log10(sort(input.data$chisq.pval, decreasing = TRUE)), xlab = "markers", ylab = "log10(p.val)", axes=F)
      axis(1); axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(log10(chisq.pval.thres),2), col = 4, lty = 2)
      text(x = input.data$n.mrk/2, y = 5, labels = paste0("Included mrks: ", sum(input.data$chisq.pval >= chisq.pval.thres)), adj = .5, col = "darkgreen")
      text(x = input.data$n.mrk/2, y = log10(chisq.pval.thres) - 5, labels = paste0("Excluded mrks: ", sum(input.data$chisq.pval < chisq.pval.thres)), adj = .5, col = "darkred")
      ANSWER <- readline("Enter 'y' to proceed or update the p value threshold: ")
      if(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
        chisq.pval.thres  <- as.numeric(ANSWER)
    }
  }
  keep<-names(which(input.data$chisq.pval >= chisq.pval.thres))
  exclude<-names(which(input.data$chisq.pval < chisq.pval.thres))
  par(op)
  structure(list(keep = keep, exclude = exclude, chisq.pval.thres = chisq.pval.thres, data.name = as.character(sys.call())[2]), class = "mappoly.chitest.seq")
}

#' get genomic order of the markers
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_genomic_order<-function(input.seq){
  if (!class(input.seq) == "mappoly.sequence")
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  if(all(is.na(input.seq$sequence.pos))){
    if(all(is.na(input.seq$sequence))) 
      stop("No sequence or sequence position information found.")
    else{
      message("Ordering markers based on sequence information")
      M<-cbind(input.seq$sequence)
      rownames(M)<-input.seq$seq.mrk.names
      return(M[order(M[,1]),, drop = FALSE])
    }
  } else if(all(is.na(input.seq$sequence))){
    if(all(is.na(input.seq$sequence.pos))) 
      stop("No sequence or sequence position information found.")
    else{
      message("Ordering markers based on sequence position information")
      M<-cbind(input.seq$sequence.pos)
      rownames(M)<-input.seq$seq.mrk.names
      return(M[order(M[,1]),, drop = FALSE])
    }
  } else{
    M<-cbind(input.seq$sequence, input.seq$sequence.pos)
    rownames(M)<-input.seq$seq.mrk.names
    return(M[order(M[,1],M[,2]),])
  }
}
