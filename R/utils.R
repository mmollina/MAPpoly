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

#' Get the number of bivalent configurations
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_w_m <- function(m){
   if(m%%2 != 0) stop("ploidy level should be an even number") 
   if(m <= 0) stop("ploidy level should be greater than zero")
   1/factorial((m/2)) * prod(choose(seq(2, m, 2),2))
}

#' Reverse map
#'
#' Provides a reverse version of a given map.
#'
#' @param input.map An object of class \code{mappoly.map}
#' @export
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
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
filter_non_conforming_classes<-function(input.data, prob.thres = NULL)
  {
  m<-input.data$m
  dp<-input.data$dosage.p
  dq<-input.data$dosage.q
  Ds<-array(NA, dim = c(m+1, m+1, m+1))
  for(i in 0:m)
    for(j in 0:m)
      Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
  Dpop<-cbind(dp,dq)
  #Gathering segregation probabilities given parental dosages
  M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
  M[M!=0]<-1
  dimnames(M)<-list(input.data$mrk.names, 0:m)
  ##if no prior probabilities
  if(nrow(input.data$geno)==input.data$n.mrk){
    for(i in 1:nrow(M)){
      id0<-!as.numeric(input.data$geno.dose[i,])%in%(which(M[i,]==1)-1)
      if(any(id0))
        input.data$geno.dose[i,id0]<-(m+1)     
    }
  return(input.data)
  }
  ## 1 represents conforming classes/ 0 represents non-conforming classes
  dp<-rep(dp, input.data$n.ind)
  dq<-rep(dq, input.data$n.ind)
  M<-M[rep(seq_len(nrow(M)), input.data$n.ind),]
  R<-input.data$geno[,-c(1:2)] - input.data$geno[,-c(1:2)]*M
  id1<-apply(R, 1, function(x) abs(sum(x))) > 0.3 # if the sum of the excluded classes is greater than 0.3, use segreg_poly
  N<-matrix(NA, sum(id1), input.data$m+1)
  ct<-1
  for(i in which(id1)){
    N[ct,] <- Ds[dp[i]+1, dq[i]+1, ]
    ct<-ct+1
  }
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
#' Exclude markers or individuals based on their missing data
#'
#' @param input.data an object of class \code{"mappoly.data"} 
#' @param type one of the following options: \code{marker} (default), 
#'             filter out markers based on their percentage of missing 
#'             data; \code{individual}, 
#'             filter out individuals based on their percentage of missing 
#'             data;
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, it plots markers or individuals vs. frequency of missing data
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing<-function(input.data, 
                             type = c("marker", "individual"), 
                             filter.thres = 0.2, 
                             inter = TRUE)
{
  type <- match.arg(type)
  switch(type,
         marker = filter_missing_mrk(input.data, 
                            filter.thres = filter.thres, 
                            inter = inter),
         individual = filter_missing_ind(input.data, 
                            filter.thres = filter.thres, 
                            inter = inter)
  )
}


#' Filter markers based on missing genotypes
#'
#' @param input.data an object of class \code{"mappoly.data"} 
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, plots markers vs. frequency of genotyped individuals
#' @keywords internal
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing_mrk<-function(input.data, filter.thres = 0.2, inter = TRUE)
{
  ANSWER <- "flag"
  mrk <- NULL
  if(interactive() && inter)
  {
    op<-par(bg = "gray", xpd = TRUE)
    while(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
    {
      na.num<-apply(input.data$geno.dose, 1, function(x,m) sum(x==m+1), m = input.data$m)
      perc.na<-na.num/input.data$n.ind
      plot(sort(perc.na), xlab = "markers", ylab = "frequency of missing data", axes = FALSE);
      axis(1);axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(filter.thres,2), col = 4, lty = 2)
      text(x = input.data$n.mrk/2, y = filter.thres + 0.05, labels = paste0("Excluded mrks: ", sum(perc.na >= filter.thres)), adj = 0, col = "darkred")
      text(x = input.data$n.mrk/2, y = filter.thres - 0.05, labels = paste0("Included mrks: ", sum(perc.na < filter.thres)), adj = 0, col = "darkgreen")
      ANSWER <- readline("Enter 'y' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
      filter.thres  <- as.numeric(ANSWER)
    }
    rm.mrks.id<-which(perc.na > filter.thres)
      if(length(rm.mrks.id)==0) return(input.data)
    rm.mrks<-names(rm.mrks.id)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(!mrk%in%rm.mrks)
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
    na.num<-apply(input.data$geno.dose, 1, function(x,m) sum(x==m+1), m = input.data$m)
    perc.na<-na.num/input.data$n.ind
    rm.mrks.id<-which(perc.na > filter.thres)
    if(length(rm.mrks.id)==0) return(input.data)
    rm.mrks<-names(rm.mrks.id)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(!mrk%in%rm.mrks)
    input.data$geno.dose<-input.data$geno.dose[-rm.mrks.id,]
    input.data$n.mrk <- nrow(input.data$geno.dose)
    input.data$mrk.names <- input.data$mrk.names[-rm.mrks.id]
    input.data$dosage.p <- input.data$dosage.p[-rm.mrks.id]
    input.data$dosage.q <- input.data$dosage.q[-rm.mrks.id]
    input.data$sequence <- input.data$sequence[-rm.mrks.id]
    if(!is.null(input.data$chisq.pval)) 
      input.data$chisq.pval <- input.data$chisq.pval[-rm.mrks.id]
    input.data$sequence.pos <- input.data$sequence.pos[-rm.mrks.id]
    return(input.data)
  }
}

#' Filter individuals based on missing genotypes 
#'
#' @param input.data an object of class \code{"mappoly.data"} 
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, plots markers vs. frequency of genotyped individuals
#' @keywords internal
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing_ind<-function(input.data, filter.thres = 0.2, inter = TRUE)
{
  ANSWER <- "flag"
  ind <- NULL
  if(interactive() && inter)
  {
    op<-par(bg = "gray", xpd = TRUE)
    while(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
    {
      na.num<-apply(input.data$geno.dose, 2, function(x,m) sum(x==m+1), m = input.data$m)
      perc.na<-na.num/input.data$n.mrk
      plot(sort(perc.na), xlab = "individuals", ylab = "frequency of missing data", axes = FALSE);
      axis(1);axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(filter.thres,2), col = 4, lty = 2)
      text(x = input.data$n.ind/2, y = filter.thres + 0.05, labels = paste0("Excluded individuals: ", sum(perc.na >= filter.thres)), adj = 0, col = "darkred")
      text(x = input.data$n.ind/2, y = filter.thres - 0.05, labels = paste0("Included individuals: ", sum(perc.na < filter.thres)), adj = 0, col = "darkgreen")
      ANSWER <- readline("Enter 'y' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1) != "y" && ANSWER !="")
        filter.thres  <- as.numeric(ANSWER)
    }
    rm.ind.id<-which(perc.na > filter.thres)
    if(length(rm.ind.id)==0) return(input.data)
    rm.ind<-names(rm.ind.id)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(!ind%in%rm.ind)
    input.data$geno.dose<-input.data$geno.dose[,-rm.ind.id]
    input.data$ind.names<-input.data$ind.names[-rm.ind.id]
    input.data$n.ind <- ncol(input.data$geno.dose)
    ##Computing chi-square p.values
    if(!is.null(input.data$chisq.pval)){
      m<-input.data$m
      Ds <- array(NA, dim = c(m+1, m+1, m+1))
      for(i in 0:m)
        for(j in 0:m)
          Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
      Dpop<-cbind(input.data$dosage.p, input.data$dosage.q)
      M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
      dimnames(M)<-list(input.data$mrk.names, c(0:m))
      M<-cbind(M, input.data$geno.dose)
      input.data$chisq.pval<-apply(M, 1, mrk_chisq_test, m = m)
    }
    par(op)
    return(input.data)
  } else {
    na.num<-apply(input.data$geno.dose, 2, function(x,m) sum(x==m+1), m = input.data$m)
    perc.na<-na.num/input.data$n.mrk
    rm.ind.id<-which(perc.na > filter.thres)
    if(length(rm.ind.id)==0) return(input.data)
    rm.ind<-names(rm.ind.id)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(!ind%in%rm.ind)
    input.data$geno.dose<-input.data$geno.dose[,-rm.ind.id]
    input.data$ind.names<-input.data$ind.names[-rm.ind.id]
    input.data$n.ind <- ncol(input.data$geno.dose)
    ##Computing chi-square p.values
    if(!is.null(input.data$chisq.pval)){
      m<-input.data$m
      Ds <- array(NA, dim = c(m+1, m+1, m+1))
      for(i in 0:m)
        for(j in 0:m)
          Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
        Dpop<-cbind(input.data$dosage.p, input.data$dosage.q)
        M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
        dimnames(M)<-list(input.data$mrk.names, c(0:m))
        M<-cbind(M, input.data$geno.dose)
        input.data$chisq.pval<-apply(M, 1, mrk_chisq_test, m = m)
    }
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


#' Filter markers based on chi-square test
#'
#' This function filter markers based on chi-square test p-values. The chi-square tests assume that markers follow the expected segregation patterns under Mendelian inheritance, only random chromosome bivalent pairing occurs and there is no double reduction.
#'
#' @param input.data name of input object (class \code{mappoly.data})
#' @param chisq.pval.thres p-value threshold used for chi-square tests
#' @param inter if TRUE (default), plots distorted vs. non-distorted markers 
#'
#' @return An object of class \code{mappoly.chitest.seq} which contains a list with the following components:
#' \item{keep}{markers that follow Mendelian segregation pattern}
#' \item{exclude}{markers with distorted segregation}
#' \item{chisq.pval.thres}{threshold p-value used for chi-square tests}
#' \item{data.name}{input dataset used to perform the chi-square tests}
#'@examples
#' \dontrun{
#'     mydata <- filter_segregation(mydata, chisq.pval.thres = 0.05/mydata$n.mrk, inter = TRUE)
#'     seq.init<-make_seq_mappoly(mrks.chi.filt)
#'}
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
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
      M<-data.frame(seq = input.seq$sequence, row.names = input.seq$seq.mrk.names)
      return(M[order(M[,1]),])
    }
  } else if(all(is.na(input.seq$sequence))){
  if(all(is.na(input.seq$sequence.pos))) 
      stop("No sequence or sequence position information found.")
    else{
      message("Ordering markers based on sequence position information")
      M<-data.frame(seq.pos = input.seq$sequence.pos, row.names = input.seq$seq.mrk.names)
      return(M[order(M[,1]),])
    }
  } else{
    M<-data.frame(seq = input.seq$sequence, 
                  seq.pos = input.seq$sequence.pos, 
                  row.names = input.seq$seq.mrk.names)
    return(M[order(M[,1], as.numeric(M[,2])),])
  }
}


#' Remove markers from a map
#' 
#' This function creates a new map by removing markers from an existing one.
#'
#' @param input.map an object of class \code{"mappoly.map"} 
#' @param mrk a vector containing markers to be removed from the input map, identified by their names or positions in the map.
#'
##' @return An object of class 'mappoly.map'
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#' \dontrun{
#' input.map <- maps.hexafake[[1]]
#' mrk <- c(4, 9, 23)
#' mrk <- c("M_85", "M_152")
#' output.map <- drop_marker(input.map, mrk)
#' plot(output.map, mrk.names = TRUE, cex = .5)
#'}
#' 
#' @export
drop_marker<-function(input.map, mrk)
{
  ## Checking class of arguments
  if (!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  ## Checking markers to be removed
  if(is.numeric(mrk) & all(mrk >= 0)){
    if(any(mrk > input.map$info$n.mrk))
      stop("At least one marker position is not contained in the map.")
  } else if(is.character(mrk)){
    if(any(!mrk%in%input.map$info$mrk.names)){
      stop("At least one marker position is not contained in the map.")
    } else {
      mrk <- which(input.map$info$mrk.names%in%mrk)
    }
  }
  ## Getting new map
  suppressMessages(output.map<-get_submap(input.map = input.map,
                                          mrk.pos = c(1:input.map$info$n.mrk)[-mrk], 
                                          phase.config = 1, 
                                          reestimate.rf = FALSE))
  if(length(input.map$maps) > 1){
    for(i in 2:length(input.map$maps)){
      suppressMessages(temp.map<-get_submap(input.map = input.map,
                                            mrk.pos = c(1:input.map$info$n.mrk)[-mrk], 
                                            phase.config = 1, 
                                            reestimate.rf = FALSE))
      output.map$maps[[i]] <- temp.map$maps[[1]]
    }
  }
  message("
    INFO:
    -----------------------------------------
    The recombination fractions provided were
    obtained using the marker positions in the 
    input map; For accurate values, plese 
    reestimate the map using functions 'reest_rf', 
    'est_full_hmm_with_global_error' or 
    'est_full_hmm_with_prior_dist'")
  return(output.map)
}

#' Add a single marker to a map
#' 
#' Creates a new map by adding a marker in a given position in a pre-built map. 
#'
#' \code{add_marker} splits the input map in two submaps to the left and to the right 
#' of the given position, and using the genotype probabilities, computes the log-likelihhod
#' of all possible linkage phases under a two-point threshold. 
#'
#' @param input.map an object of class \code{"mappoly.map"} 
#' @param mrk the name of the marker to be inserted
#' @param pos the name of the marker after which the new marker should be added.
#'            One also can inform the numeric position (between markers) were the 
#'            new marker should be added. To insert a marker at the beginning of a 
#'            map, one should use \code{pos = 0}.
#' @param rf.matrix an object of class \code{"mappoly.rf.matrix"}
#' @param genoprob an object of class \code{"mappoly.genoprob"}
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration.
#' @param tol the desired accuracy.
#'                     
##' @return An object of class 'mappoly.map'
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#' \dontrun{
#' input.map <- maps.hexafake[[1]]
#' mrk <- c(4, 9, 23)
#' mrk <- c("M_85", "M_152")
#' output.map <- drop_marker(input.map, mrk)
#' plot(output.map, mrk.names = TRUE, cex = .5)
#'}
#' 
#' @export
add_marker <- function(input.map, 
                       mrk,
                       pos,
                       rf.matrix, 
                       genoprob = NULL,
                       phase.config = "best",
                       tol = 10e-4){
  ## Checking class of arguments
  if(!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  if(!inherits(mrk, "character")) {
    stop("Please, provide the marker name")
  }
  if(!inherits(rf.matrix, "mappoly.rf.matrix")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.rf.matrix'")
  }
  if(is.null(genoprob)){
    message("Calculating genoprob.")
    genoprob <- calc_genoprob(input.map)
  }
  if(!inherits(genoprob, "mappoly.genoprob")) {
    stop(deparse(substitute(genoprob)), " is not an object of class 'mappoly.genoprob'")
  }
  if(any(names(genoprob$map)!=input.map$info$mrk.names)){
    warning(deparse(substitute(genoprob)), "is inconsistent with 'input.map'.\nRecalculating genoprob.")
    genoprob <- calc_genoprob(input.map)
  }
  ## ploidy
  m <- input.map$info$m
  ## Number of genotypes in the offspring
  ngen <- dim(genoprob$probs)[1]
  ## number of markers
  nmrk <- dim(genoprob$probs)[2]
  ## number of individuals
  nind <- dim(genoprob$probs)[3]
  ## number of genotipic states
  ngam <- choose(m, m/2)
  ## the threshold for visiting states: 1/ngen
  thresh.cut.path <- 1/ngen

  ## Indexing position were the marker should be inserted
  if(is.numeric(pos)){
    if(!(pos >=0 & pos <= (nmrk+1))){
      stop(deparse(substitute(pos)), " out of bounds.")
    } 
  } else if(is.character(pos)){
    if(!pos%in%input.map$info$mrk.names){
      stop("The reference marker is not contained in the map.")
    } else {
      pos <- which(input.map$info$mrk.names%in%pos)
    }
  }
  
  ## Hash table: homolog combination --> states to visit in both parents
  A<-as.matrix(expand.grid(c(1:ngam) - 1, c(1:ngam) - 1)[,2:1])
  rownames(A) <- dimnames(genoprob$probs)[[1]]
  
  ## Indexing marker to be inserted
  if(!mrk%in%colnames(rf.matrix$rec.mat)){
    stop(deparse(substitute(mrk)), " is not in 'rf.matrix'")
  }
  mrk.id <- as.numeric(colnames(rf.matrix$ShP)[match(mrk, colnames(rf.matrix$rec.mat))])
  
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  }  else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  
  ## Adding marker: beginning of sequence
  if(pos == 0){
    ## h: states to visit in both parents
    ## e: probability distribution (ignored in this version) 
    e.right <- h.right <- vector("list", nind)
    for(i in 1:nind){
      a.right <-  genoprob$probs[,pos+1,i]
      e.right[[i]] <- a.right[a.right > thresh.cut.path]
      h.right[[i]] <- A[names(e.right[[i]]), , drop = FALSE]
    }
    r.test<-generate_all_link_phases_elim_equivalent_haplo(block1 = input.map$maps[[i.lpc]], 
                                                           block2 = mrk.id, 
                                                           rf.matrix = rf.matrix, 
                                                           m = m, max.inc = 0)
  } else if(pos > 0 & pos < nmrk){   ## Adding marker: middle positions
    ## h: states to visit in both parents
    ## e: probability distribution (ignored in this version) 
    e.left <- h.left <- e.right <- h.right <- vector("list", nind)
    for(i in 1:nind){
      a.left <- genoprob$probs[,pos,i]  
      a.right <-  genoprob$probs[,pos+1,i]
      e.left[[i]] <- a.left[a.left > thresh.cut.path]
      h.left[[i]] <- A[names(e.left[[i]]), , drop = FALSE]
      e.right[[i]] <- a.right[a.right > thresh.cut.path]
      h.right[[i]] <- A[names(e.right[[i]]), , drop = FALSE]
    }
    ## Info to the left
    suppressMessages(map.left <- get_submap(input.map = input.map,
                                            mrk.pos = 1:pos,
                                            phase.config = i.lpc, 
                                            reestimate.rf = FALSE))
    r.left <- generate_all_link_phases_elim_equivalent_haplo(block1 = map.left$maps[[i.lpc]], 
                                                             block2 = mrk.id, 
                                                             rf.matrix = rf.matrix, 
                                                             m = m, max.inc = 0)
    ## Info to the right
    suppressMessages(map.right <- get_submap(input.map = input.map,
                                             mrk.pos = (pos+1):input.map$info$n.mrk,
                                             phase.config = i.lpc, 
                                             reestimate.rf = FALSE))
    r.right <- generate_all_link_phases_elim_equivalent_haplo(block1 = map.right$maps[[i.lpc]], 
                                                              block2 = mrk.id, 
                                                              rf.matrix = rf.matrix, 
                                                              m = m, max.inc = 0)
    ## Using both information sources, left and right 
    r.test <- unique(r.left, r.right)
  } else if(pos == nmrk){
    ## h: states to visit in both parents
    ## e: probability distribution (ignored in this version) 
    e.left <- h.left <- vector("list", nind)
    for(i in 1:nind){
      a.left <- genoprob$probs[,pos,i]  
      e.left[[i]] <- a.left[a.left > thresh.cut.path]
      h.left[[i]] <- A[names(e.left[[i]]), , drop = FALSE]
    }
    r.test<-generate_all_link_phases_elim_equivalent_haplo(block1 = input.map$maps[[i.lpc]], 
                                                           block2 = mrk.id, 
                                                           rf.matrix = rf.matrix, 
                                                           m = m, max.inc = 0)
  } else stop("should not get here!")
  ## gathering maps to test and conditional probabilities
  test.maps <- mrk.genoprobs <- vector("list", length(r.test))
  for(i in 1:length(test.maps))
  {
    suppressMessages(hap.temp <- get_submap(input.map, 
                                            c(1,1), 
                                            reestimate.rf = FALSE))
    hap.temp$maps[[1]]$seq.num<-rep(mrk.id, 2)
    hap.temp$maps[[1]]$seq.ph <- list(P = c(r.test[[i]]$P, r.test[[i]]$P),
                                      Q = c(r.test[[i]]$Q, r.test[[i]]$Q))
    hap.temp$maps[[1]]$seq.rf <- 10e-6
    test.maps[[i]] <- hap.temp
    ## States to visit for inserted biallelic SNP 
    mrk.genoprobs[[i]] <- calc_genoprob(test.maps[[i]], verbose = FALSE)
  }
  ## Inserted marker  
  ## h: states to visit in both parents
  ## e: probability distribution 
  h <- e <- vector("list", length(r.test))
  for(j in 1:length(r.test)){
    etemp<-htemp<-vector("list", nind)
    for(i in 1:nind){
      a <- mrk.genoprobs[[j]]$probs[,1,i]  
      etemp[[i]] <- a[a > thresh.cut.path]
      htemp[[i]] <- A[names(etemp[[i]]), , drop = FALSE]
    }
    h[[j]] <- htemp
    e[[j]] <- etemp
  }
  configs<-vector("list", length(test.maps))
  names(configs) <- paste0("Phase_config.", 1:length(test.maps))
  res<-matrix(NA, nrow = length(h), ncol = 3, dimnames = list(names(configs), c("log_like", "rf1", "rf2")))
  ## Computing three-point multiallelic map using HMM
  for(i in 1:length(h)){
    cat(".")
    if(pos == 0){
      h.test<-c(h[i], list(h.right))
      names(h.test) <- c("hap1", "hap2")
      e.test<-c(e[i], list(e.right))
      restemp<-est_haplo_hmm(m = m, 
                             n.mrk = length(h.test), 
                             n.ind = nind, 
                             haplo = h.test, 
                             #emit = e.test, 
                             rf_vec = rep(0.01, length(h.test)-1), 
                             verbose = FALSE, 
                             use_H0 = FALSE, 
                             tol = tol) 
      temp<-unlist(restemp)
      res[i,1:length(temp)] <- temp
      P<-c(test.maps[[i]]$maps[[1]]$seq.ph$P[1],
           input.map$maps[[1]]$seq.ph$P)
      Q<-c(test.maps[[i]]$maps[[1]]$seq.ph$Q[1], 
           input.map$maps[[1]]$seq.ph$Q)
      names(P)<-names(Q)<-c(test.maps[[i]]$maps[[1]]$seq.num[1], 
                            input.map$maps[[1]]$seq.num)
      configs[[i]]<-list(P = P, Q = Q)
    } else if(pos > 0 & pos < nmrk){
      h.test<-c(list(h.left), h[i], list(h.right))
      names(h.test) <- c("hap1", "hap2", "hap3")
      e.test<-c(list(e.left), e[i], list(e.right))
      restemp<-est_haplo_hmm(m = m, 
                             n.mrk = length(h.test), 
                             n.ind = nind, 
                             haplo = h.test, 
                             #emit = e.test, 
                             rf_vec = rep(0.01, length(h.test)-1), 
                             verbose = FALSE, 
                             use_H0 = FALSE, 
                             tol = tol) 
      temp<-unlist(restemp)
      res[i,1:length(temp)] <- temp
      P<-c(map.left$maps[[1]]$seq.ph$P, 
           test.maps[[i]]$maps[[1]]$seq.ph$P[1],
           map.right$maps[[1]]$seq.ph$P)
      Q<-c(map.left$maps[[1]]$seq.ph$Q, 
           test.maps[[i]]$maps[[1]]$seq.ph$Q[1], 
           map.right$maps[[1]]$seq.ph$Q)
      names(P)<-names(Q)<-c(map.left$maps[[1]]$seq.num, 
                            test.maps[[i]]$maps[[1]]$seq.num[1], 
                            map.right$maps[[1]]$seq.num)
      configs[[i]]<-list(P = P, Q = Q)
    } else if(pos == nmrk){
      h.test<-c(list(h.left), h[i])
      names(h.test) <- c("hap1", "hap2")
      e.test<-c(list(e.left), e[i])
      restemp<-est_haplo_hmm(m = m, 
                             n.mrk = length(h.test), 
                             n.ind = nind, 
                             haplo = h.test, 
                             #emit = e.test, 
                             rf_vec = rep(0.01, length(h.test)-1), 
                             verbose = FALSE, 
                             use_H0 = FALSE, 
                             tol = tol) 
      temp<-unlist(restemp)
      res[i,1:length(temp)] <- temp
      P<-c(input.map$maps[[1]]$seq.ph$P, 
           test.maps[[i]]$maps[[1]]$seq.ph$P[1])
      Q<-c(input.map$maps[[1]]$seq.ph$Q, 
           test.maps[[i]]$maps[[1]]$seq.ph$Q[1])
      names(P)<-names(Q)<-c(input.map$maps[[1]]$seq.num, 
                            test.maps[[i]]$maps[[1]]$seq.num[1])
      configs[[i]]<-list(P = P, Q = Q)
    }
  }
  cat("\n")
  res<-res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
  return(list(stats = res, phase.config  = configs))  
}

#' Data sanity check
#'
#' Checks the consistency of dataset
#'
#' @param x an object of class \code{mappoly.data}
#' 
#' @examples
#' \dontrun{
#' 
#'     solcap.dose.file <- system.file('extdata', 'tetra_solcap_geno', package = 'mappoly')
#'     dat.dose <- read_geno(file.in  = solcap.dose.file)
#'     check_data_sanity(dat.dose)
#' 
#'     solcap.file <- system.file('extdata', 'tetra_solcap_geno_dist.bz2', package = 'mappoly')
#'     dat.dist <- read_geno_dist(file.in  = solcap.file)
#'     check_data_sanity(dat.dist)
#'}
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'     
#' @keywords internal
#' @export check_data_sanity
check_data_sanity<-function(x){
  if(ncol(x$geno) == x$n.ind)
    check_data_dose_sanity(x)
  else if(ncol(x$geno) == x$m + 3)
    check_data_dist_sanity(x)
  else
    stop("Inconsistent genotypic information.")
}

#' Checks the consistency of dataset (dosage)
#'
#' @param void interfunction to be documented
#' @keywords internal
check_data_dose_sanity <- function(x){
  test<-logical(24L)
  names(test) <- 1:24

  # ploidy
  test[1] <- x$m%%2 != 0
  test[2] <- any(sapply(x$dosage.p, function(y) max(y) > x$m | min(y) < 0))
  test[3] <- any(sapply(x$dosage.q, function(y) max(y) > x$m | min(y) < 0))
  test[4] <- max(x$geno.dose) > x$m + 1
  test[5] <- min(x$geno.dose) < 0
  
  # number of individuals
  test[6] <- x$n.ind < 0
  test[7] <- length(x$ind.names) != x$n.ind
  test[8] <- ncol(x$geno.dose) != x$n.ind
  
  # number of markers
  test[9] <- x$n.mrk < 0
  test[10] <- length(x$mrk.names) != x$n.mrk
  test[11] <- length(x$dosage.p) != x$n.mrk
  test[12] <- length(x$dosage.q) != x$n.mrk
  test[13] <- length(x$sequence) != x$n.mrk
  test[14] <- length(x$sequence.pos) != x$n.mrk
  test[15] <- nrow(x$geno.dose) != x$n.mrk
  test[16] <- length(x$chisq.pval) != x$n.mrk
  
  # individual names in the probability dataset
  test[17] <- !is.character(x$ind.names)
  test[18] <- !identical(colnames(x$geno.dose), x$ind.names)
  
  # individual names in the dosage dataset
  test[19] <- !is.character(x$mrk.names)
  test[20] <- !identical(rownames(x$geno.dose), x$mrk.names)
  
  # dosage in both parents
  test[21] <- !is.integer(x$dosage.p)
  test[22] <- !is.integer(x$dosage.q)
  test[23] <- !identical(names(x$dosage.p), x$mrk.names)
  test[24] <- !identical(names(x$dosage.q), x$mrk.names)
  
  if(any(test))
    return(test)
  else
    return(0)
}

#' Checks the consistency of dataset (probability distribution)
#'
#' @param void interfunction to be documented
#' @keywords internal
check_data_dist_sanity <- function(x){
  test<-logical(29L)
  names(test) <- 1:29
  # ploidy
  test[1] <- x$m%%2 != 0
  test[2] <- any(sapply(x$dosage.p, function(y) max(y) > x$m | min(y) < 0))
  test[3] <- any(sapply(x$dosage.q, function(y) max(y) > x$m | min(y) < 0))
  test[4] <- ncol(x$geno) > x$m + 3
  test[5] <- max(x$geno.dose) > x$m + 1
  test[6] <- min(x$geno.dose) < 0
  
  # number of individuals
  test[7] <- x$n.ind < 0
  test[8] <- length(x$ind.names) != x$n.ind
  test[9] <- length(unique(x$geno$ind)) != x$n.ind
  test[10] <- ncol(x$geno.dose) != x$n.ind
  
  # number of markers
  test[11] <- x$n.mrk < 0
  test[12] <- length(x$mrk.names) != x$n.mrk
  test[13] <- length(x$dosage.p) != x$n.mrk
  test[14] <- length(x$dosage.q) != x$n.mrk
  test[15] <- length(x$sequence) != x$n.mrk
  test[16] <- length(x$sequence.pos) != x$n.mrk
  test[17] <- nrow(x$geno)/x$n.ind != x$n.mrk
  test[18] <- nrow(x$geno.dose) != x$n.mrk
  test[19] <- length(x$chisq.pval) != x$n.mrk
  
  # individual names in the probability dataset
  test[20] <- !is.character(x$ind.names)
  test[21] <- !identical(x$geno$ind, rep(x$ind.names, each = x$n.mrk))
  test[22] <- !identical(colnames(x$geno.dose), x$ind.names)
  
  # individual names in the dosage dataset
  test[23] <- !is.character(x$mrk.names)
  test[24] <- !identical(x$geno$mrk, rep(x$mrk.names, x$n.ind))
  test[25] <- !identical(rownames(x$geno.dose), x$mrk.names)
  
  # dosage in both parents
  test[26] <- !is.integer(x$dosage.p)
  test[27] <- !is.integer(x$dosage.q)
  test[28] <- !identical(names(x$dosage.p), x$mrk.names)
  test[29] <- !identical(names(x$dosage.q), x$mrk.names)
  
  if(any(test))
    return(test)
  else
    return(0)
}
