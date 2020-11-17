#' Extract the LOD Scores in a \code{'mappoly.map'} object
#' @param x an object of class \code{mappoly.map}
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
#' @param void internal function to be documented
#' @keywords internal
#' @export
get_rf_from_mat <- function(M){
  r<-numeric(nrow(M)-1)
  for(i in 1:(nrow(M)-1)){
    r[i]<-M[i, i+1]
  }
  r
}


#' Is it a probability dataset?
#'
#' @param void internal function to be documented
#' @keywords internal
is.prob.data <- function(x){
  exists('geno', where = x)
}

#' Get the number of bivalent configurations
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
get_w_m <- function(m){
  if(m%%2 != 0) stop("ploidy level should be an even number") 
  if(m <= 0) stop("ploidy level should be greater than zero")
  1/factorial((m/2)) * prod(choose(seq(2, m, 2),2))
}

#' Reverse map
#'
#' Provides the reverse of a given map.
#'
#' @param input.map an object of class \code{mappoly.map}
#' 
#' @examples
#' plot_genome_vs_map(solcap.mds.map[[1]])
#' plot_genome_vs_map(rev_map(solcap.mds.map[[1]]))
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#'@export
#'
rev_map<-function(input.map)
{
  output.map<-input.map
  output.map$info <- lapply(output.map$info, rev)
  for(i in 1:length(output.map$maps))
  {
    output.map$maps[[i]]$seq.num <- rev(input.map$maps[[i]]$seq.num)
    output.map$maps[[i]]$seq.rf <- rev(input.map$maps[[i]]$seq.rf)
    output.map$maps[[i]]$seq.ph$P <- rev(input.map$maps[[i]]$seq.ph$P)
    output.map$maps[[i]]$seq.ph$Q <- rev(input.map$maps[[i]]$seq.ph$Q)
  }
  return(output.map)
}

#' Returns the class with the highest probability in 
#' a genotype probability distribution
#'
#' @param geno the probabilistic genotypes contained in the object
#'     \code{'mappoly.data'}
#' @param prob.thres probability threshold to select the genotype. 
#'     Values below this genotype are assumed as missing data
#' @return a matrix containing the doses of each genotype and
#'     marker. Markers are disposed in rows and individuals are 
#'     disposed in columns. Missing data are represented by NAs
#' @keywords internal
#' @examples
#' \donttest{
#' geno.dose <- dist_prob_to_class(hexafake.geno.dist$geno)
#' geno.dose$geno.dose[1:10,1:10]
#'}   
#' @importFrom magrittr "%>%"
#' @importFrom reshape2 melt dcast
#' @importFrom dplyr group_by filter arrange
#' @export
dist_prob_to_class <- function(geno, prob.thres = 0.9) {
  a<-reshape2::melt(geno, id.vars = c("mrk", "ind"))
  mrk <- ind <- value <- variable <- NULL # Setting the variables to NULL first
  a$variable<-as.numeric(levels(a$variable))[a$variable]
  b<-a %>%
    dplyr::group_by(mrk, ind) %>%
    dplyr::filter(value > prob.thres) %>%
    dplyr::arrange(mrk, ind, variable)
  z<-reshape2::dcast(data = b[,1:3], formula = mrk ~ ind, value.var = "variable")
  rownames(z)<-z[,"mrk"]
  z<-data.matrix(frame = z[,-1])
  n<-setdiff(unique(geno$mrk), rownames(z))
  if(length(n) > 0)
  {
    m<-matrix(NA, nrow = length(n), ncol = ncol(z), dimnames = list(n, colnames(z)))
    z<-rbind(z,m)
  }
  rm.ind<-setdiff(unique(geno$ind), colnames(z))
  flag <- FALSE
  if(length(rm.ind) > 0){
    flag <- TRUE
    warning("Inividual(s) ", paste(rm.ind, collapse = " "), 
            "\n  did not meet the 'prob.thres' criteria for any of\n  the markers and was (were) removed.")
    geno <- geno %>% dplyr::filter(ind %in% colnames(z))
  }
  z <- z[as.character(unique(geno$mrk)), as.character(unique(geno$ind))]
  list(geno.dose = z, geno = geno, flag = flag)
}

#' Export data to \code{polymapR}
#' @param data.in an object of class \code{mappoly.data}
#' @return a dosage \code{matrix} 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @examples
#' \donttest{
#' require(polymapR)
#' dat<-export_data_to_polymapR(hexafake)
#' F1checked <- checkF1(dosage_matrix = dat, 
#'                      parent1 = "P1",
#'                      parent2 = "P2",
#'                      F1 = colnames(dat)[-c(1:2)],
#'                      polysomic = TRUE, 
#'                      disomic = FALSE, 
#'                      mixed = FALSE, 
#'                      ploidy = 6)
#'  head(F1checked$checked_F1)
#'}  
#' @export export_data_to_polymapR
export_data_to_polymapR <- function(data.in)
{
  data.out<-as.matrix(data.frame(P1 = data.in$dosage.p,     
                                 P2 = data.in$dosage.q,
                                 data.in$geno.dose))
  data.out[data.out == (data.in$m + 1)] <- NA
  return(data.out)
}

#' Msg function
#'
#' @param void internal function to be documented
#' @keywords internal
#' @importFrom cli rule
msg <- function(text, line = 1){
  cli::rule(line = line, right = text) %>%
    text_col() %>%
    message()
}

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
#' @param void internal function to be documented
#' @keywords internal
#' @export
mf_k <- function(d) 0.5 * tanh(d/50)
#'
#' Map functions
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
mf_h <- function(d) 0.5 * (1 - exp(-d/50))
#' Map functions
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
mf_m <- function(d) sapply(d, function(a) min(a/100, 0.5))
#' Map functions
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
imf_k <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  50 * atanh(2 * r)
}
#' Map functions
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
imf_h <- function(r) {
  r[r >= 0.5] <- 0.5 - 1e-14
  -50 * log(1 - 2 * r)
}
#' Map functions
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
imf_m <- function(r) sapply(r, function(a) min(a * 100, 50))


#' Compare two polyploid haplotypes stored in list format
#'
#' @param m ploidy level
#' @param h1 homology group 1
#' @param h2 homology group 2
#' @return a numeric vector of size \code{m} indicating which
#'     homolog in h2 represents the homolog in h1. If there is
#'     no correspondence, i.e. different homolog, it returns NA for
#'     that homolog.
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
#' @param void internal function to be documented
#' @keywords internal
#' @export plot_compare_haplotypes
plot_compare_haplotypes <- function(m, hom.allele.p1, hom.allele.q1, hom.allele.p2 = NULL, hom.allele.q2 = NULL) {
  nmmrk<-names(hom.allele.p1)
  oldpar <- par(mar = c(5.1, 4.1, 4.1, 2.1))
  on.exit(par(oldpar))
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
}

#' Check if it is possible to estimate the recombination
#' fraction between neighbor markers using two-point
#' estimation
#'
#' @param void internal function to be documented
#' @keywords internal
pos_twopt_est<-function(input.seq){
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
#' @param void internal function to be documented
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
#' @param void internal function to be documented
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

#' Color pallet ggplot-like
#'
#' @param void internal function to be documented
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
#' @param input.data an object of class \code{mappoly.data}
#' @param prob.thres probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than 'prob.thres' 
#'     are considered as missing data for the dosage calling purposes
#' @examples
#' \donttest{
#' data.updated = update_missing(hexafake.geno.dist, prob.thres = 0.5)
#' print(hexafake.geno.dist)
#' print(data.updated)
#' }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
update_missing<-function(input.data, prob.thres = 0.95){
  geno.dose <- dist_prob_to_class(geno = input.data$geno, prob.thres = prob.thres)
  if(geno.dose$flag)
  {
    geno <- geno.dose$geno
    geno.dose <- geno.dose$geno.dose
  } else {
    geno.dose <- geno.dose$geno.dose
  }
  geno.dose[is.na(geno.dose)] <- input.data$m + 1
  input.data$geno.dose<-geno.dose
  input.data$prob.thres<-prob.thres
  return(input.data)
}

#' Chi-square test
#'
#' @param void internal function to be documented
#' @keywords internal
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

#' Get the genomic position of markers in a sequence
#'
#' This functions gets the genomic position of markers in a sequence and
#' return an ordered data frame with the name and position of each marker
#' @param input.seq a sequence object of class \code{mappoly.sequence}
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @examples
#'  s1<-make_seq_mappoly(tetra.solcap, "all")
#'  o1<-get_genomic_order(s1)
#'  head(o1)
#' @export
get_genomic_order<-function(input.seq, verbose = TRUE){
  if (!inherits(input.seq, "mappoly.sequence")) {
    stop(deparse(substitute(input.seq)), 
         " is not an object of class 'mappoly.sequence'")
  }
  if(all(is.na(input.seq$sequence.pos))){
    if(all(is.na(input.seq$sequence))) 
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on sequence information")
      M<-data.frame(seq = input.seq$sequence, row.names = input.seq$seq.mrk.names)
      return(M[order(M[,1]),])
    }
  } else if(all(is.na(input.seq$sequence))){
    if(all(is.na(input.seq$sequence.pos))) 
      stop("No sequence or sequence position information found.")
    else{
      if (verbose) message("Ordering markers based on sequence position information")
      M<-data.frame(seq.pos = input.seq$sequence.pos, row.names = input.seq$seq.mrk.names)
      return(M[order(as.numeric(M[,1])),])
    }
  } else{
    M<-data.frame(seq = input.seq$sequence, 
                  seq.pos = input.seq$sequence.pos, 
                  row.names = input.seq$seq.mrk.names)
    return(M[order(as.numeric(M[,1]), as.numeric(M[,2])),])
  }
}

#' Remove markers from a map
#' 
#' This function creates a new map by removing markers from an existing one.
#'
#' @param input.map an object of class \code{mappoly.map}
#' @param mrk a vector containing markers to be removed from the input map, 
#'            identified by their names or positions
#' @return an object of class \code{mappoly.map}
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#'  sub.map<-get_submap(maps.hexafake[[1]], 1:50, reestimate.rf = FALSE)
#'  plot(sub.map, mrk.names = TRUE)
#'  mrk.to.remove <- c("M_1", "M_23", "M_34")
#'  red.map <- drop_marker(sub.map, mrk.to.remove)
#'  plot(red.map, mrk.names = TRUE)
#' 
#' @export
drop_marker<-function(input.map, mrk, verbose = TRUE)
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
  if (verbose) message("
    INFO:
    -----------------------------------------
    The recombination fractions provided were
    obtained using the marker positions in the 
    input map; For accurate values, plese 
    reestimate the map using functions 'reest_rf', 
    'est_full_hmm_with_global_error' or 
    'est_full_hmm_with_prior_prob'")
  return(output.map)
}

#' Add a single marker to a map
#' 
#' Creates a new map by adding a marker in a given position in a pre-built map. 
#'
#' \code{add_marker} splits the input map into two sub-maps to the left and the 
#' right of the given position. Using the genotype probabilities, it computes 
#' the log-likelihood of all possible linkage phases under a two-point threshold
#' inherited from function \code{\link[mappoly]{rf_list_to_matrix}}. 
#'
#' @param input.map an object of class \code{mappoly.map}
#'  
#' @param mrk the name of the marker to be inserted
#' 
#' @param pos the name of the marker after which the new marker should be added.
#'            One also can inform the numeric position (between markers) were the 
#'            new marker should be added. To insert a marker at the beginning of a 
#'            map, use \code{pos = 0}
#'            
#' @param rf.matrix an object of class \code{mappoly.rf.matrix} containing the recombination
#'                  fractions and the number of homologues sharing alleles between pairwise 
#'                  markers on \code{input.map}. It is important that \code{shared.alleles = TRUE}
#'                  in function \code{\link[mappoly]{rf_list_to_matrix}} when computing \code{rf.matrix}.
#' 
#' @param genoprob an object of class \code{mappoly.genoprob} containing the genotype probabilities
#' for all marker positions on \code{input.map}
#' 
#' @param phase.config which phase configuration should be used. "best" (default) 
#'                     will choose the maximum likelihood configuration
#'                     
#' @param tol the desired accuracy (default = 10e-04)
#' 
#' @param r.test for internal use only
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return A list of class \code{mappoly.map} with two elements: 
#' 
#' i) info:  a list containing information about the map, regardless of the linkage phase configuration:
#' \item{m}{the ploidy level}
#' \item{n.mrk}{number of markers}
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{mrk.names}{the names of markers in the map}
#' \item{seq.dose.p}{a vector containing the dosage in parent 1 for all markers in the map}
#' \item{seq.dose.q}{a vector containing the dosage in parent 2 for all markers in the map}
#' \item{sequence}{a vector indicating the sequence (usually chromosome) each marker belongs 
#'                 as informed in the input file. If not available, 
#'                 \code{sequence = NULL}}
#' \item{sequence.pos}{physical position (usually in megabase) of the markers into the sequence}
#' \item{seq.ref}{reference base used for each marker (i.e. A, T, C, G). If not available, 
#'                 \code{seq.ref = NULL}}                 
#' \item{seq.alt}{alternative base used for each marker (i.e. A, T, C, G). If not available, 
#'                 \code{seq.ref = NULL}}
#' \item{chisq.pval}{a vector containing p-values of the chi-squared test of Mendelian 
#'                   segregation for all markers in the map}                 
#' \item{data.name}{name of the dataset of class \code{mappoly.data}}
#' \item{ph.thres}{the LOD threshold used to define the linkage phase configurations to test}
#' 
#' ii) a list of maps with possible linkage phase configuration. Each map in the list is also a 
#'    list containing
#' \item{seq.num}{a vector containing the (ordered) indices of markers in the map, 
#'                according to the input file}
#' \item{seq.rf}{a vector of size (\code{n.mrk - 1}) containing a sequence of recombination 
#'               fraction between the adjacent markers in the map}
#' \item{seq.ph}{linkage phase configuration for all markers in both parents}
#' \item{loglike}{the hmm-based multipoint likelihood}
#' 
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @examples
#' \donttest{
#' sub.map<-get_submap(maps.hexafake[[1]], 1:20, reestimate.rf = FALSE)
#' plot(sub.map, mrk.names = TRUE)
#' s<-make_seq_mappoly(hexafake, sub.map$info$mrk.names)
#' tpt <- est_pairwise_rf(s)
#' rf.matrix <- rf_list_to_matrix(input.twopt = tpt,
#'                                thresh.LOD.ph = 3, 
#'                                thresh.LOD.rf = 3,
#'                                shared.alleles = TRUE)
#' ###### Removing marker "M_1" (first) #######
#' mrk.to.remove <- "M_1"
#' input.map <- drop_marker(sub.map, mrk.to.remove)
#' plot(input.map, mrk.names = TRUE)
#' ## Computing conditional probabilities using the resulting map
#' genoprob <- calc_genoprob(input.map)
#' res.add.M_1<-add_marker(input.map = input.map,
#'                         mrk = "M_1",
#'                         pos = 0,
#'                         rf.matrix = rf.matrix,
#'                         genoprob = genoprob,
#'                         tol = 10e-4)  
#'  plot(res.add.M_1, mrk.names = TRUE)                       
#'  best.phase <- res.add.M_1$maps[[1]]$seq.ph
#'  names.id<-names(best.phase$P)
#'  plot_compare_haplotypes(m = 6,
#'                          hom.allele.p1 = best.phase$P[names.id],
#'                          hom.allele.q1 = best.phase$Q[names.id],
#'                          hom.allele.p2 = sub.map$maps[[1]]$seq.ph$P[names.id],
#'                          hom.allele.q2 = sub.map$maps[[1]]$seq.ph$Q[names.id])
#'                          
#' ###### Removing marker "M_10" (middle or last) #######
#' mrk.to.remove <- "M_10"
#' input.map <- drop_marker(sub.map, mrk.to.remove)
#' plot(input.map, mrk.names = TRUE)
#' # Computing conditional probabilities using the resulting map
#' genoprob <- calc_genoprob(input.map)
#' res.add.M_10<-add_marker(input.map = input.map,
#'                         mrk = "M_10",
#'                         pos = "M_9",
#'                         rf.matrix = rf.matrix,
#'                         genoprob = genoprob,
#'                         tol = 10e-4)  
#'  plot(res.add.M_10, mrk.names = TRUE)                       
#'  best.phase <- res.add.M_10$maps[[1]]$seq.ph
#'  names.id<-names(best.phase$P)
#'  plot_compare_haplotypes(m = 6,
#'                          hom.allele.p1 = best.phase$P[names.id],
#'                          hom.allele.q1 = best.phase$Q[names.id],
#'                          hom.allele.p2 = sub.map$maps[[1]]$seq.ph$P[names.id],
#'                          hom.allele.q2 = sub.map$maps[[1]]$seq.ph$Q[names.id]) 
#' }
#' @export
add_marker <- function(input.map,  mrk, pos, rf.matrix, genoprob = NULL, 
                       phase.config = "best", tol = 10e-4, r.test = NULL, verbose = TRUE){
  ## Checking class of arguments
  if(!inherits(input.map, "mappoly.map")) {
    stop(deparse(substitute(input.map)), " is not an object of class 'mappoly.map'")
  }
  if(!inherits(mrk, "character")) {
    stop("Please, provide the marker name")
  }
  if(!inherits(rf.matrix, "mappoly.rf.matrix")) {
    stop(deparse(substitute(rf.matrix)), " is not an object of class 'mappoly.rf.matrix'")
  }
  if(is.null(rf.matrix$ShP))
    stop(deparse(substitute(rf.matrix)), " should contain the number of homologues sharing alleles.
    Use 'shared.alleles = TRUE' in function 'rf_list_to_matrix'")
  if(!all(c(input.map$info$mrk.names, mrk)%in%colnames(rf.matrix$rec.mat))){
    stop(deparse(substitute(rf.matrix)), " does not contain all necessary information about 'input.map' and 'mrk'.")
  }
  ## choosing the linkage phase configuration
  LOD.conf <- get_LOD(input.map, sorted = FALSE)
  if(phase.config == "best") {
    i.lpc <- which.min(LOD.conf)
  }  else if (phase.config > length(LOD.conf)) {
    stop("invalid linkage phase configuration")
  } else i.lpc <- phase.config
  
  ## Checking genoprob
  if(is.null(genoprob)){
    if (verbose) message("Calculating genoprob.")
    genoprob <- calc_genoprob(input.map, phase.config = i.lpc, verbose = FALSE)
  }
  if(!inherits(genoprob, "mappoly.genoprob")) {
    stop("'", deparse(substitute(genoprob)), "' is not an object of class 'mappoly.genoprob'")
  }
  if(!identical(names(genoprob$map), input.map$info$mrk.names)){
    warning("'", deparse(substitute(genoprob)), "' is inconsistent with 'input.map'.\n  Recalculating genoprob.")
    genoprob <- calc_genoprob(input.map, phase.config = i.lpc, verbose = FALSE)
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
    if(is.null(r.test)){
      r.test<-generate_all_link_phases_elim_equivalent_haplo(block1 = input.map$maps[[i.lpc]], 
                                                             block2 = mrk.id, 
                                                             rf.matrix = rf.matrix, 
                                                             m = m, max.inc = 0)
      
    }
  } 
  else if(pos > 0 & pos < nmrk){   ## Adding marker: middle positions
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
  } 
  else if(pos == nmrk){
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
  } 
  else stop("should not get here!")
  ## gathering maps to test and conditional probabilities
  test.maps <- mrk.genoprobs <- vector("list", length(r.test))
  for(i in 1:length(test.maps))
  {
    ## This sub-map is just to create a framework
    suppressMessages(hap.temp <- get_submap(input.map, 
                                            c(1,1), 
                                            reestimate.rf = FALSE))
    hap.temp <- filter_map_at_hmm_thres(hap.temp, thres.hmm = 10e-10)
    hap.temp$maps[[1]]$seq.num<-rep(mrk.id, 2)
    hap.temp$maps[[1]]$seq.ph <- list(P = c(r.test[[i]]$P, r.test[[i]]$P),
                                      Q = c(r.test[[i]]$Q, r.test[[i]]$Q))
    hap.temp$maps[[1]]$seq.rf <- 10e-6
    hap.temp$info$mrk.names <- rep(mrk,2)
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
    if (verbose) cat(".")
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
    } 
    else if(pos > 0 & pos < nmrk){
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
    } 
    else if(pos == nmrk){
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
  if (verbose) cat("\n")
  res<-res[order(res[,"log_like"], decreasing = TRUE),,drop = FALSE]
  
  
  ## Updating map
  output.map <- input.map
  seq.num<-as.numeric(names(configs[[1]]$P))
  output.map$info$seq.num <- seq.num
  output.map$info$mrk.names <- colnames(rf.matrix$rec.mat)[match(seq.num, colnames(rf.matrix$ShP))]
  output.map$info$n.mrk <- length(output.map$info$mrk.names)
  output.map$maps <- vector("list", nrow(res))
  for(i in 1:nrow(res))
  {
    ## Updating recombination fractions (approximated)
    if(pos == 0){
      seq.rf <- as.numeric(c(res[i, "rf1"], input.map$maps[[i.lpc]]$seq.rf))
    } else if(pos > 0 & pos < nmrk){
      seq.rf <- as.numeric(c(head(input.map$maps[[i.lpc]]$seq.rf, n = pos - 1),
                  res[i, c("rf1", "rf2")], 
                  tail(input.map$maps[[i.lpc]]$seq.rf, n = input.map$info$n.mrk - pos - 1)))
    } else if(pos == nmrk){
      seq.rf <- as.numeric(c(input.map$maps[[i.lpc]]$seq.rf, res[i, "rf1"]))
    }
    output.map$maps[[i]] <- list(seq.num = seq.num, 
                                 seq.rf = seq.rf, 
                                 seq.ph = configs[[rownames(res)[i]]],
                                 loglike = res[i, "log_like"])
  }
  dat<-get(input.map$info$data.name, pos = 1)
  output.map$info$seq.dose.p <- dat$dosage.p[output.map$info$mrk.names]
  output.map$info$seq.dose.q <- dat$dosage.q[output.map$info$mrk.names]
  output.map$info$sequence <- dat$sequence[output.map$info$mrk.names]
  output.map$info$sequence.pos <- dat$sequence.pos[output.map$info$mrk.names]
  output.map$info$chisq.pval <- dat$chisq.pval[output.map$info$mrk.names]
  return(output.map)
}

#' Data sanity check
#'
#' Checks the consistency of a dataset
#'
#' @param x an object of class \code{mappoly.data}
#' 
#' @return if consistent, returns 0. If not consistent, returns a 
#'         vector with a number of tests, where \code{TRUE} indicates
#'         a failed test.
#' 
#' @examples
#' check_data_sanity(tetra.solcap)
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378}
#'     
#' @export check_data_sanity
check_data_sanity<-function(x){
  if(exists('geno', where = x)){
    check_data_dist_sanity(x)
  } else if (exists('geno.dose', where = x)){
    check_data_dose_sanity(x)
  } else
    stop("Inconsistent genotypic information.")
}

#' Checks the consistency of dataset (dosage)
#'
#' @param void internal function to be documented
#' @keywords internal
check_data_dose_sanity <- function(x){
  test<-logical(24L)
  names(test) <- 1:24
  
  # ploidy
  test[1] <- x$m%%2 != 0 #is ploidy even?
  test[2] <- any(sapply(x$dosage.p, function(y) max(y) > x$m | min(y) < 0)) #are dosages in P higher than ploidy?
  test[3] <- any(sapply(x$dosage.q, function(y) max(y) > x$m | min(y) < 0)) #are dosages in Q higher than ploidy?
  test[4] <- max(x$geno.dose) > x$m + 1 #is there any dose in offspring higher than ploidy?
  test[5] <- min(x$geno.dose) < 0 #is there any negative dose in offspring?
  
  # number of individuals
  test[6] <- x$n.ind < 0 #is the number of individuals greater than zero?
  test[7] <- length(x$ind.names) != x$n.ind #is the number of individual names equal to the number of individuals?
  test[8] <- ncol(x$geno.dose) != x$n.ind #is the number of columns in the dosage matrix equal to the number of individuals?
  
  # number of markers
  test[9] <- x$n.mrk < 0 #is the number of markers greater than zero?
  test[10] <- length(x$mrk.names) != x$n.mrk #is the number of marker names equal to the number of markers?
  test[11] <- length(x$dosage.p) != x$n.mrk #is the number of marker dosages in P equal to the number of markers?
  test[12] <- length(x$dosage.q) != x$n.mrk #is the number of marker dosages in Q equal to the number of markers?
  if(length(x$sequence) > 0)
    test[13] <- length(x$sequence) != x$n.mrk #is the number of sequences equal to the number of markers?
  if(length(x$sequence.pos) > 0)
    test[14] <- length(x$sequence.pos) != x$n.mrk #is the number of sequence positions equal to the number of markers?
  test[15] <- nrow(x$geno.dose) != x$n.mrk #is the number of rows in the dosage matrix equal to the number of markers?
  if(length(x$chisq.pval) > 0)
    test[16] <- length(x$chisq.pval) != x$n.mrk #is the number of chi-square tests equal to the number of markers?
  
  # individual names in the dosage dataset
  test[17] <- !is.character(x$ind.names) # are individual's names characters
  test[18] <- !identical(colnames(x$geno.dose), x$ind.names) # are column names in dosage matrix identical to individual names?
  
  # markers names in the dosage dataset
  test[19] <- !is.character(x$mrk.names)# are marker's names characters
  test[20] <- !identical(rownames(x$geno.dose), x$mrk.names)# are row names in dosage matrix identical to marker names?
  
  # dosage in both parents
  test[21] <- !is.integer(x$dosage.p) # are dosages in P numeric
  test[22] <- !is.integer(x$dosage.q) # are dosages in Q numeric 
  test[23] <- !identical(names(x$dosage.p), x$mrk.names) # are names in P's dosage vector identical to marker names?
  test[24] <- !identical(names(x$dosage.q), x$mrk.names) # are names in Q's dosage vector identical to marker names?
  
  if(any(test))
    return(test)
  else
    return(0)
}

#' Checks the consistency of dataset (probability distribution)
#'
#' @param void internal function to be documented
#' @keywords internal
check_data_dist_sanity <- function(x){
  test<-logical(29L)
  names(test) <- 1:29
  
  # ploidy
  test[1] <- x$m%%2 != 0 #is ploidy even?
  test[2] <- any(sapply(x$dosage.p, function(y) max(y) > x$m | min(y) < 0)) #are dosages in P higher than ploidy?
  test[3] <- any(sapply(x$dosage.q, function(y) max(y) > x$m | min(y) < 0)) #are dosages in Q higher than ploidy?
  test[4] <- ncol(x$geno) > x$m + 3 #is the number of columns in the probability data frame correct? (ploidy + 3)
  test[5] <- max(x$geno.dose) > x$m + 1 #is there any dose in offspring higher than ploidy?
  test[6] <- min(x$geno.dose) < 0 #is there any negative dose in offspring?
  
  # number of individuals
  test[7] <- x$n.ind < 0 #is the number of individuals greater than zero?
  test[8] <- length(x$ind.names) != x$n.ind #is the number of individual names equal to the number of individuals?
  test[9] <- length(unique(x$geno$ind)) != x$n.ind #is the number of individuals in the probability data frame equal to the number of individuals?
  test[10] <- ncol(x$geno.dose) != x$n.ind #is the number of columns in the dosage matrix equal to the number of individuals?
  
  # number of markers
  test[11] <- x$n.mrk < 0 #is the number of markers greater than zero?
  test[12] <- length(x$mrk.names) != x$n.mrk #is the number of marker names equal to the number of markers?
  test[13] <- length(x$dosage.p) != x$n.mrk #is the number of marker dosages in P equal to the number of markers?
  test[14] <- length(x$dosage.q) != x$n.mrk #is the number of marker dosages in Q equal to the number of markers?
  if(length(x$sequence) > 0)
    test[15] <- length(x$sequence) != x$n.mrk #is the number of sequences equal to the number of markers?
  if(length(x$sequence.pos) > 0)
    test[16] <- length(x$sequence.pos) != x$n.mrk #is the number of sequence positions equal to the number of markers?
  test[17] <- nrow(x$geno)/x$n.ind != x$n.mrk#is the number of rows in the probability data frame divided by n.ind equal to the number of markers?
  test[18] <- nrow(x$geno.dose) != x$n.mrk#is the number of rows in the dosage matrix equal to the number of markers?
  if(length(x$chisq.pval) > 0)
    test[19] <- length(x$chisq.pval) != x$n.mrk#is the number of chi-square tests equal to the number of markers?
  
  # individual names in the dosage and probability dataset
  test[20] <- !is.character(x$ind.names) # are individual's names characters
  test[21] <- !identical(x$geno$ind, rep(x$ind.names, each = x$n.mrk)) #are individual's names in the probability data frame properly 
                                                                       #arranged and consistent with the informed individual's names 
  test[22] <- !identical(colnames(x$geno.dose), x$ind.names)# are column names in dosage matrix identical to individual names?
  
  # marker names in the dosage and probability dataset
  test[23] <- !is.character(x$mrk.names)# are marker's names characters
  test[24] <- !identical(x$geno$mrk, rep(x$mrk.names, x$n.ind))#are marker names in the probability data frame properly 
                                                               #arranged and consistent with the informed marker names 
  test[25] <- !identical(rownames(x$geno.dose), x$mrk.names)# are row names in dosage matrix identical to marker names?
  
  # dosage in both parents
  test[26] <- !is.integer(x$dosage.p)# are dosages in P numeric
  test[27] <- !is.integer(x$dosage.q)# are dosages in Q numeric 
  test[28] <- !identical(names(x$dosage.p), x$mrk.names)# are names in P's dosage vector identical to marker names?
  test[29] <- !identical(names(x$dosage.q), x$mrk.names)# are names in Q's dosage vector identical to marker names?
  
  if(any(test))
    return(test)
  else
    return(0)
}

#' Merge datasets
#'
#' This function merges two datasets of class \code{mappoly.data}. This can be useful
#' when individuals of a population were genotyped using two or more techniques
#' and have datasets in different files or formats. Please notice that the datasets 
#' should contain the same number of individuals and they must be represented identically 
#' in both datasets  (e.g. \code{Ind_1} in both datasets, not \code{Ind_1}
#' in one dataset and \code{ind_1} or \code{Ind.1} in the other).
#'
#' @param dat.1 the first dataset of class \code{mappoly.data} to be merged
#'
#' @param dat.2 the second dataset of class \code{mappoly.data} to be merged (default = NULL);
#' if \code{dat.2 = NULL}, the function returns \code{dat.1} only
#'
#' @return An object of class \code{mappoly.data} which contains all markers
#' from both datasets. It will be a list with the following components:
#'     \item{m}{ploidy level}
#'     \item{n.ind}{number individuals}
#'     \item{n.mrk}{total number of markers}
#'     \item{ind.names}{the names of the individuals}
#'     \item{mrk.names}{the names of the markers}
#'     \item{dosage.p}{a vector containing the dosage in
#'       parent P for all \code{n.mrk} markers}
#'     \item{dosage.q}{a vector containing the dosage in
#'       parent Q for all \code{n.mrk} markers}
#'     \item{sequence}{a vector indicating which sequence each marker
#'       belongs. Zero indicates that the marker was not assigned to any
#'       sequence}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{seq.ref}{if one or both datasets originated from read_vcf, it keeps
#' reference alleles from sequencing platform, otherwise is NULL}
#'     \item{seq.alt}{if one or both datasets originated from read_vcf, it keeps
#' alternative alleles from sequencing platform, otherwise is NULL}
#'     \item{all.mrk.depth}{if one or both datasets originated from read_vcf, it keeps
#' marker read depths from sequencing, otherwise is NULL}
#'     \item{prob.thres}{(unused field)}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{geno}{if both datasets contain genotype distribution information,
#' the final object will contain 'geno'. This is set to NULL otherwise}
#'     \item{nphen}{(0)}
#'     \item{phen}{(NULL)}
#'     \item{chisq.pval}{a vector containing p-values related to the chi-squared 
#'     test of Mendelian segregation performed for all markers in both datasets}
#'     \item{kept}{if elim.redundant=TRUE when reading any dataset, holds all non-redundant markers}
#'     \item{elim.correspondence}{if elim.redundant=TRUE when reading any dataset,
#' holds all non-redundant markers and its equivalence to the redundant ones}
#' 
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @examples
#' \donttest{
#' ## Loading a subset of SNPs from chromosomes 3 and 12 of sweetpotato dataset 
#' ## (SNPs anchored to Ipomoea trifida genome)
#' dat <- NULL
#' for(i in c(3, 12)){
#'   cat("Loading chromosome", i, "...\n")
#'     tempfl <- tempfile(pattern = paste0("ch", i), fileext = ".vcf.gz")
#'     x <- "https://github.com/mmollina/MAPpoly_vignettes/raw/master/data/sweet_sample_ch"
#'     address <- paste0(x, i, ".vcf.gz")
#'     download.file(url = address, destfile = tempfl)
#'     dattemp <- read_vcf(file = tempfl, parent.1 = "PARENT1", parent.2 = "PARENT2",
#'                         ploidy = 6, verbose = FALSE)
#'     dat <- merge_datasets(dat, dattemp)
#'   cat("\n")
#' }
#' dat
#' plot(dat)
#'}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @export merge_datasets
#' @importFrom dplyr bind_rows arrange
merge_datasets = function(dat.1 = NULL, dat.2 = NULL){
  ## Check objects class
  if (is.null(dat.1)){
    if (is.null(dat.2)) return(dat.1)
    if (!inherits(dat.2, "mappoly.data")) {
      stop(deparse(substitute(dat.2)), " is not an object of class 'mappoly.data'")
    } else {
      return(dat.2)
    }
  } 
  if (!inherits(dat.1, "mappoly.data")) {
    stop(deparse(substitute(dat.1)), " is not an object of class 'mappoly.data'")
  }
  if (is.null(dat.2)) return(dat.1)
  if (!inherits(dat.2, "mappoly.data")) {
    stop(deparse(substitute(dat.2)), " is not an object of class 'mappoly.data'")
  }
  ## Check ploidy
  if (dat.1$m != dat.2$m){
    stop("The ploidy levels in the datasets do not match. Please check your files and try again.")
  }
  ## Check individuals
  if (dat.1$n.ind != dat.2$n.ind){
    stop("The datasets have different number of individuals. Please check your files and try again.")
  }
  ## Check individual names
  if (!all(dat.1$ind.names %in% dat.2$ind.names)){
    nmi.1 = dat.1$ind.names[!(dat.1$ind.names %in% dat.2$ind.names)]
    nmi.2 = dat.2$ind.names[!(dat.2$ind.names %in% dat.1$ind.names)]
    cat("Some individuals' names don't match:\n")
    cat("Dataset 1\tDataset 2\n")
    for (i in 1:length(nmi.1)) cat(nmi.1[i], "\t\t", nmi.2[i], "\n")
    stop("Your datasets have the same number of individuals, but some of them are not the same. Please check the list above, fix your files and try again.")
  }
  ## Checking (and fixing) individuals order in geno.dose
  if (!all(colnames(dat.1$geno.dose) == colnames(dat.2$geno.dose))){
    dat.2$geno.dose = dat.2$geno.dose[,colnames(dat.1$geno.dose)]
  }
  ## Merging all items
  dat.1$geno.dose = rbind(dat.1$geno.dose, dat.2$geno.dose)
  dat.1$dosage.p = c(dat.1$dosage.p, dat.2$dosage.p)
  dat.1$dosage.q = c(dat.1$dosage.q, dat.2$dosage.q)
  ##dat.1$mrk.names = c(dat.1$mrk.names, dat.2$mrk.names) ## Fixing possible name incompatibilities  
  dat.1$mrk.names = rownames(dat.1$geno.dose) ## Fixing possible name incompatibilities  
  dat.1$n.mrk = dat.1$n.mrk + dat.2$n.mrk
  dat.1$sequence = c(dat.1$sequence, dat.2$sequence)
  dat.1$sequence.pos = c(dat.1$sequence.pos, dat.2$sequence.pos)
  # VCF info
  ## seq.ref and seq.alt
  if (!is.null(dat.1$seq.ref) && !is.null(dat.2$seq.ref)){
    dat.1$seq.ref = c(dat.1$seq.ref,dat.2$seq.ref)
    dat.1$seq.alt = c(dat.1$seq.alt,dat.2$seq.alt)  
  }  else if (is.null(dat.1$seq.ref) && is.null(dat.2$seq.ref)){
    dat.1$seq.ref = dat.1$seq.alt = NULL
  }  else if (!is.null(dat.1$seq.ref) && is.null(dat.2$seq.ref)){
    dat.1$seq.ref = c(dat.1$seq.ref,rep(NA,length(dat.2$n.mrk)))
    dat.1$seq.alt = c(dat.1$seq.alt,rep(NA,length(dat.2$n.mrk)))
  } else if (is.null(dat.1$seq.ref) && !is.null(dat.2$seq.ref)){
    dat.1$seq.ref = c(rep(NA,length(dat.1$n.mrk)),dat.2$seq.ref)
    dat.1$seq.alt = c(rep(NA,length(dat.1$n.mrk)),dat.2$seq.alt)
  }
  ## all.mrk.depth
  if (!is.null(dat.1$all.mrk.depth) && !is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = c(dat.1$all.mrk.depth,dat.2$all.mrk.depth)
  }  else if (is.null(dat.1$all.mrk.depth) && is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = NULL
  }  else if (!is.null(dat.1$all.mrk.depth) && is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = c(dat.1$all.mrk.depth,rep(NA,length(dat.2$n.mrk)))
  } else if (is.null(dat.1$all.mrk.depth) && !is.null(dat.2$all.mrk.depth)){
    dat.1$all.mrk.depth = c(rep(NA,length(dat.1$n.mrk)),dat.2$all.mrk.depth)
  }
  # Chisq info
  if (!is.null(dat.1$chisq.pval) && !is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = c(dat.1$chisq.pval,dat.2$chisq.pval)
  }  else if (is.null(dat.1$chisq.pval) && is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = NULL
  }  else if (!is.null(dat.1$chisq.pval) && is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = c(dat.1$chisq.pval,rep(NA,length(dat.2$n.mrk)))
  } else if (is.null(dat.1$chisq.pval) && !is.null(dat.2$chisq.pval)){
    dat.1$chisq.pval = c(rep(NA,length(dat.1$n.mrk)),dat.2$chisq.pval)
  }
  # Redundant info
  ## Kept info
  if (!is.null(dat.1$kept) && !is.null(dat.2$kept)){
    dat.1$kept = c(dat.1$kept,dat.2$kept)
  }  else if (is.null(dat.1$kept) && is.null(dat.2$kept)){
    dat.1$kept = NULL
  }  else if (!is.null(dat.1$kept) && is.null(dat.2$kept)){
    dat.1$kept = c(dat.1$kept, dat.2$mrk.names)
  } else if (is.null(dat.1$kept) && !is.null(dat.2$kept)){
    dat.1$kept = c(dat.1$mrk.names,dat.2$kept)
  }  
  ## elim.correspondence info
  if (!is.null(dat.1$elim.correspondence) && !is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = rbind(dat.1$elim.correspondence,dat.2$elim.correspondence)
  }  else if (is.null(dat.1$elim.correspondence) && is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = NULL
  }  else if (!is.null(dat.1$elim.correspondence) && is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = dat.1$elim.correspondence
  } else if (is.null(dat.1$elim.correspondence) && !is.null(dat.2$elim.correspondence)){
    dat.1$elim.correspondence = dat.2$elim.correspondence
  }  
  ## geno dist info (just keep if both datasets contain this information)
  if (exists('geno', where = dat.1) && exists('geno', where = dat.2)){
    #dat.1$geno = rbind(dat.1$geno,dat.2$geno)
    ind <- NULL
    dat.1$geno <- dplyr::bind_rows(dat.1$geno, dat.2$geno) %>% arrange(ind)
  } else dat.1$geno = NULL
  if(!is.null(dat.1$chisq.pval) | !any(is.na(dat.1$chisq.pval)))
    dat.1$chisq.pval <- dat.1$chisq.pval[dat.1$mrk.names]
  ## Fixing possible issues with marker names
  names(dat.1$dosage.p) = names(dat.1$dosage.q) = names(dat.1$sequence.pos) = names(dat.1$sequence) = names(dat.1$chisq.pval) = dat.1$mrk.names
  return(dat.1)
}

#' Summary maps
#'
#' This function generates a brief summary table of a list of \code{mappoly.map} objects
#' @param map.list a list of objects of class \code{mappoly.map}
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' @return a data frame containing a brief summary of all maps contained in \code{map.list}
#' @examples
#' tetra.sum <- summary_maps(solcap.err.map)
#' tetra.sum
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @export summary_maps
summary_maps = function(map.list, verbose = TRUE){
  ## Check data
  if (any(!sapply(map.list, inherits, "mappoly.map"))) 
    stop(deparse(substitute(map.list)), 
         " is not a list containing 'mappoly.map' objects.")
  results = data.frame("LG" = as.character(seq(1,length(map.list),1)),
                       "Genomic sequence" = as.character(unlist(lapply(map.list, function(x) paste(unique(x$info$sequence), collapse = "-")))),
                       "Map size (cM)" = unlist(lapply(map.list, function(x) round(sum(c(0, imf_h(x$maps[[1]]$seq.rf))), 2))),
                       "Markers/cM" = round(unlist(lapply(map.list, function(x) x$info$n.mrk/(round(sum(c(0, imf_h(x$maps[[1]]$seq.rf))), 2)))),2),
                       "Simplex" = unlist(lapply(map.list, function(x) sum(get_tab_mrks(x)[rbind(c(1,2),c(2,1),c(x$info$m,(x$info$m+1)),c((x$info$m+1),x$info$m))]))),
                       "Double-simplex" = unlist(lapply(map.list, function(x) sum(get_tab_mrks(x)[rbind(c(2,2),c(x$info$m,x$info$m))]))),
                       "Multiplex" = unlist(lapply(map.list, function(x) sum(get_tab_mrks(x)) - sum(get_tab_mrks(x)[rbind(c(2,2),c(x$info$m,x$info$m))]) - sum(get_tab_mrks(x)[rbind(c(1,2),c(2,1),c(x$info$m,(x$info$m+1)),c((x$info$m+1),x$info$m))]))),
                       "Total" = unlist(lapply(map.list, function(x) x$info$n.mrk)),
                       "Max gap" = unlist(lapply(map.list, function(x) round(imf_h(max(x$maps[[1]]$seq.rf)),2))),
                       check.names = FALSE, stringsAsFactors = F)
  results = rbind(results, c('Total', NA, sum(as.numeric(results$`Map size (cM)`)), round(mean(as.numeric(results$`Markers/cM`)),2), sum(as.numeric(results$Simplex)), sum(as.numeric(results$`Double-simplex`)), sum(as.numeric(results$Multiplex)), sum(as.numeric(results$Total)), round(mean(as.numeric(results$`Max gap`)),2)))
  if (verbose){
    all.mrks = unlist(lapply(map.list, function(x) return(x$info$mrk.names)))
    if (!any(get(map.list[[1]]$info$data.name, pos = 1)$elim.correspondence$elim %in% all.mrks))
      message("\nYour dataset contains removed (redundant) markers. Once finished the maps, remember to add them back with the function 'update_map'.\n")
  }
  return(results)
}

#' Get table of dosage combinations
#' 
#' Internal function
#' @param x an object of class \code{mappoly.map}
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
get_tab_mrks = function(x){
  tab = table(get(x$info$data.name, pos = 1)$dosage.p[which(get(x$info$data.name, pos = 1)$mrk.names %in% x$info$mrk.names)], get(x$info$data.name, pos = 1)$dosage.q[which(get(x$info$data.name, pos = 1)$mrk.names %in% x$info$mrk.names)])
  doses = as.character(seq(0,x$info$m,1))
  # checking dimensions
  if (!all(doses %in% colnames(tab))){
    newmat = matrix(NA, nrow(tab), length(doses))
    newmat[,which(doses %in% colnames(tab))] = tab[,which(colnames(tab) %in% doses)]
    newmat[,which(!doses %in% colnames(tab))] = 0
    rownames(newmat) = rownames(tab)
    colnames(newmat) = doses
    tab = newmat
  }
  if (!all(doses %in% rownames(tab))){
    newmat = matrix(NA, length(doses), ncol(tab))
    newmat[which(doses %in% rownames(tab)),] = tab[which(rownames(tab) %in% doses),]
    newmat[which(!doses %in% rownames(tab)),] = 0
    colnames(newmat) = colnames(tab)
    rownames(newmat) = doses
    tab = newmat
  }
  return(tab)
}

#' Update map
#' 
#' This function takes an object of class \code{mappoly.map} and checks for
#' removed redundant markers in the original dataset. Once redundant markers
#' are found, they are re-added to the map in their respective equivalent positions
#' and another HMM round is performed.
#' 
#' @param input.maps a single map or a list of maps of class \code{mappoly.map}
#' @param verbose if TRUE (default), shows information about each update process
#' @return an updated map (or list of maps) of class \code{mappoly.map}, containing the original map(s) plus redundant markers
#' @author Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#' @examples
#' orig.map <- solcap.err.map
#' up.map <- lapply(solcap.err.map, update_map)
#' summary_maps(orig.map)
#' summary_maps(up.map)
#' @export update_map
#' 
update_map = function(input.maps, verbose = TRUE){
  ## Checking object
  if (inherits(input.maps, "mappoly.map"))
    input.maps = list(input.maps)
  if(!inherits(input.maps, "list"))
    stop(deparse(substitute(input.maps)), 
         " is not an object of class 'mappoly.map' neither a list containing 'mappoly.map' objects.")
  if (any(!sapply(input.maps, inherits, "mappoly.map"))) 
    stop(deparse(substitute(input.maps)), 
         " is not an object of class 'mappoly.map' neither a list containing 'mappoly.map' objects.")
  ## Checking the existence of redundant markers
  if (is.null(get(input.maps[[1]]$info$data.name, pos = 1)$elim.correspondence))
    stop('Your dataset does not contain redundant markers. Please check it and try again.')
  ## Creating list to handle results
  results = list()
  for (i in 1:length(input.maps)){
    if (verbose) cat("Updating map", i , "\n")    
    input.map = input.maps[[i]]
    ## Checking if redundant markers belong to the informed map
    map.kept.mrks = which(as.character(get(input.map$info$data.name, pos = 1)$elim.correspondence$kept) %in% input.map$info$mrk.names)
    if (is.null(map.kept.mrks)){
      if (verbose) cat("There is no redundant marker on map ",i,". Skipping it.\n")
      next
    }
    corresp = get(input.map$info$data.name, pos = 1)$elim.correspondence[which(as.character(get(input.map$info$data.name, pos = 1)$elim.correspondence$kept) %in% input.map$info$mrk.names),]
    ## Check if redundant markers were not already added to the map
    if (any(corresp$elim %in% input.map$info$mrk.names)) {
      if (verbose) cat("Some redundant markers were already added to the map ", i,". These markers will be skipped.\n")
      corresp = corresp[!(corresp$elim %in% input.map$info$mrk.names),]
    }
    ## Updating number of markers
    input.map$info$n.mrk = input.map$info$n.mrk + nrow(corresp)
    ## Adding markers to the sequence
    while (nrow(corresp) > 0){
      pos.kep = match(as.character(corresp$kept), get(input.map$info$data.name, pos = 1)$mrk.names)
      input.map$info$seq.num = append(input.map$info$seq.num, NA, after = which(input.map$info$seq.num == pos.kep[1]))
      input.map$info$chisq.pval = append(input.map$info$chisq.pval, input.map$info$chisq.pval[which(input.map$info$seq.num == pos.kep[1])], after = which(input.map$info$seq.num == pos.kep[1]))
      input.map$info$seq.dose.p = append(input.map$info$seq.dose.p, input.map$info$seq.dose.p[which(input.map$info$seq.num == pos.kep[1])], after = which(input.map$info$seq.num == pos.kep[1]))
      input.map$info$seq.dose.q = append(input.map$info$seq.dose.q, input.map$info$seq.dose.q[which(input.map$info$seq.num == pos.kep[1])], after = which(input.map$info$seq.num == pos.kep[1]))
      input.map$info$sequence = as.numeric(append(input.map$info$sequence, as.character(corresp$sequence[1]), after = which(input.map$info$seq.num == pos.kep[1])))
      input.map$info$sequence.pos = as.numeric(append(input.map$info$sequence.pos, as.character(corresp$sequence.pos[1]), after = which(input.map$info$seq.num == pos.kep[1])))
      input.map$info$mrk.names = append(input.map$info$mrk.names, as.character(corresp$elim[1]), after = which(input.map$info$mrk.names == as.character(corresp$kept[1])))
      if (!is.null(input.map$info$seq.ref) && !is.null(input.map$info$seq.alt)){
        input.map$info$seq.ref = append(input.map$info$seq.ref, as.character(corresp$seq.ref[1]), after = which(input.map$info$seq.num == pos.kep[1]))
        input.map$info$seq.alt = append(input.map$info$seq.alt, as.character(corresp$seq.alt[1]), after = which(input.map$info$seq.num == pos.kep[1]))
      }
      for (j in 1:length(input.map$maps)){
        input.map$maps[[j]]$seq.rf = append(input.map$maps[[j]]$seq.rf, 0.000001, after = which(input.map$info$seq.num == pos.kep[1]))
        input.map$maps[[j]]$seq.ph$P = append(input.map$maps[[j]]$seq.ph$P, input.map$maps[[j]]$seq.ph$P[paste0(pos.kep[1])], after = which(input.map$info$seq.num == pos.kep[1]))
        input.map$maps[[j]]$seq.ph$Q = append(input.map$maps[[j]]$seq.ph$Q, input.map$maps[[j]]$seq.ph$Q[paste0(pos.kep[1])], after = which(input.map$info$seq.num == pos.kep[1]))
      }
      corresp = corresp[-1,]
    }
    ## Renaming objects and updating map information
    names(input.map$info$sequence) = names(input.map$info$sequence.pos) = names(input.map$info$chisq.pval) = names(input.map$info$seq.dose.p) = names(input.map$info$seq.dose.q) = names(input.map$info$seq.num) = input.map$info$mrk.names
    for (j in 1:length(input.map$maps)){
      names(input.map$maps[[j]]$seq.ph$P) = names(input.map$maps[[j]]$seq.ph$Q) = as.character(input.map$info$seq.num)
      input.map$maps[[j]]$seq.num = input.map$info$seq.num
    }
    if (!is.null(input.map$info$seq.ref) && !is.null(input.map$info$seq.alt))
      names(input.map$info$seq.ref) = names(input.map$info$seq.alt) = input.map$info$mrk.names
    results[[i]] = input.map
  }
  if (length(results) == 1) return(results[[1]])
  else return(results)
}

#' Random sampling of dataset
#' @param x an object  of class \code{mappoly.data}
#' @param n number of individuals or markers to be sampled
#' @param percentage if \code{n==NULL}, the percentage of individuals or markers to be sampled
#' @param type should sample individuals or markers?
#' @param selected.ind a vector containing the name of the individuals to select. Only has effect 
#' if \code{type = "individual"}, \code{n = NULL} and \code{percentage = NULL}
#' @param selected.mrk a vector containing the name of the markers to select. Only has effect 
#' if \code{type = "marker"}, \code{n = NULL} and \code{percentage = NULL}
#' @return an object  of class \code{mappoly.data}
#' @keywords internal
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
sample_data <- function(input.data, n = NULL, 
                        percentage = NULL, 
                        type = c("individual", "marker"),
                        selected.ind = NULL,
                        selected.mrk = NULL){
  type <- match.arg(type)
  if(type == "individual"){
    if(!is.null(n)){
      selected.ind.id <- sort(sample(input.data$n.ind, n))
    } else if(!is.null(percentage)){
      selected.ind.id <- sample(input.data$n.ind, ceiling(input.data$n.ind * percentage/100))
    } else if(!is.null(selected.ind)){
      selected.ind.id <- match(selected.ind, input.data$ind.names)
    } else {
      stop("Inform 'n', 'percentage' or selected.ind.")
    }
    ind <- mrk <- NULL
    selected.ind <- input.data$ind.names[selected.ind.id]
    if(length(selected.ind.id) >= input.data$n.ind) return(input.data)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(ind%in%selected.ind)
    input.data$geno.dose<-input.data$geno.dose[,selected.ind.id]
    input.data$ind.names<-input.data$ind.names[selected.ind.id]
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
  else if(type == "marker"){
    if(!is.null(n)){
      selected.mrk.id <- sort(sample(input.data$n.mrk, n))
    } else if(!is.null(percentage)){
      selected.mrk.id <- sort(sample(input.data$n.mrk, ceiling(input.data$n.mrk * percentage/100)))
    } else if(!is.null(selected.mrk)){
      selected.mrk.id <- match(selected.mrk, input.data$mrk.names)
    } else{
      stop("Inform 'n', 'percentage' or selected.mrk.")
    }
    selected.mrk <- input.data$mrk.names[selected.mrk.id]
    if(length(selected.mrk.id) >= input.data$n.mrk) return(input.data)
    if(nrow(input.data$geno)!=input.data$n.mrk)
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(mrk%in%selected.mrk)
    input.data$geno.dose<-input.data$geno.dose[selected.mrk.id,]
    input.data$n.mrk <- nrow(input.data$geno.dose)
    input.data$mrk.names <- input.data$mrk.names[selected.mrk.id]
    input.data$dosage.p <- input.data$dosage.p[selected.mrk.id]
    input.data$dosage.q <- input.data$dosage.q[selected.mrk.id]
    input.data$sequence <- input.data$sequence[selected.mrk.id]
    input.data$sequence.pos <- input.data$sequence.pos[selected.mrk.id]
    input.data$seq.ref <- input.data$seq.ref[selected.mrk.id]
    input.data$seq.alt <- input.data$seq.alt[selected.mrk.id]
    input.data$all.mrk.depth <- input.data$all.mrk.depth[selected.mrk.id]
    input.data$kept <- intersect(input.data$mrk.names, input.data$kept)
    input.data$elim.correspondence <- input.data$elim.correspondence[input.data$elim.correspondence$kept%in%input.data$mrk.names,]
    input.data$chisq.pval <- input.data$chisq.pval[names(input.data$chisq.pval)%in%input.data$mrk.names]
    if(!is.null(input.data$chisq.pval)) 
      input.data$chisq.pval <- input.data$chisq.pval[names(input.data$chisq.pval)%in%input.data$mrk.names]
    return(input.data)  
  }
  else stop("Inform type")
}


#' Get weighted ordinary least squared map give a sequence and rf matrix
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
#' @importFrom zoo na.approx
get_ols_map <- function(input.seq, input.map, weight = TRUE){
  id <- input.seq$seq.mrk.names
  y <- as.numeric((imf_h(as.dist(input.map$rec.mat[id,id]))))
  w <- as.numeric((imf_h(as.dist(input.map$lod.mat[id,id]))))
  v <- t(combn(id,2))
  rf <- get_rf_from_mat(input.map$rec.mat[id,id])
  rf <- zoo::na.approx(rf)
  z<-cumsum(imf_h(c(0,rf)))
  names(z)<-id
  x<-numeric(nrow(v))
  names(x)<-names(y)<-apply(v, 1, paste0, collapse="-")
  for(i in 1:nrow(v))
    x[i]<-z[v[i,2]]-z[v[i,1]]
  if(weight)
    model <- lm(y ~ x-1, weights=w)
  else
    model <- lm(y ~ x-1)
  new <- data.frame(x = z)
  u<-predict(model, new)
  d <- cumsum(imf_h(c(0, mf_h(diff(u)))))
  names(d) <- id
  d
}

#' Get dosage type in a sequence
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
get_dosage_type <- function(input.seq){
  p <- abs(abs(input.seq$seq.dose.p - input.seq$m/2) - input.seq$m/2)
  q <- abs(abs(input.seq$seq.dose.q - input.seq$m/2) - input.seq$m/2)
  s.p <- p == 1 & q == 0
  s.q <- p == 0 & q == 1
  ds <- p == 1 & q == 1
  list(simplex.p = input.seq$seq.mrk.names[s.p],
       simplex.q = input.seq$seq.mrk.names[s.q], 
       double.simplex = input.seq$seq.mrk.names[ds],
       multiplex = input.seq$seq.mrk.names[!(s.p | s.q | ds)])
}

#' Aggregate matrix cells (lower the resolution by a factor)
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
aggregate_matrix <- function(M, fact){
  id <- seq(1,ncol(M), by = fact)
  id <- cbind(id, c(id[-1]-1, ncol(M)))
  R<-matrix(NA, nrow(id), nrow(id))
  for(i in 1:(nrow(id)-1)){
    for(j in (i+1):nrow(id)){
      R[j,i] <-  R[i,j] <- mean(M[id[i,1]:id[i,2], id[j,1]:id[j,2]], na.rm = TRUE)
    }
  }
  R
}








