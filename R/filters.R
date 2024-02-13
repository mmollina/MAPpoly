#' Filter non-conforming classes in F1, non double reduced population.
#' @keywords internal
filter_non_conforming_classes <- function(input.data, prob.thres = NULL)
{
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  ploidy <- input.data$ploidy
  dp <- input.data$dosage.p1
  dq <- input.data$dosage.p2
  Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
  for(i in 0:ploidy)
    for(j in 0:ploidy)
      Ds[i+1,j+1,] <- segreg_poly(ploidy = ploidy, dP = i, dQ = j)
  Dpop <- cbind(dp,dq)
  #Gathering segregation probabilities given parental dosages
  M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
  M[M != 0] <- 1
  dimnames(M) <- list(input.data$mrk.names, 0:ploidy)
  ##if no prior probabilities
  if(!is.prob.data(input.data)){
    for(i in 1:nrow(M)){
      id0 <- !as.numeric(input.data$geno.dose[i,])%in%(which(M[i,] == 1)-1)
      if(any(id0))
        input.data$geno.dose[i,id0] <- (ploidy+1)
    }
    return(input.data)
  }
  ## 1 represents conforming classes/ 0 represents non-conforming classes
  dp <- rep(dp, input.data$n.ind)
  dq <- rep(dq, input.data$n.ind)
  M <- M[rep(seq_len(nrow(M)), input.data$n.ind),]
  R <- input.data$geno[,-c(1:2)] - input.data$geno[,-c(1:2)]*M
  id1 <- apply(R, 1, function(x) abs(sum(x))) > 0.3 # if the sum of the excluded classes is greater than 0.3, use segreg_poly
  N <- matrix(NA, sum(id1), input.data$ploidy+1)
  ct <- 1
  for(i in which(id1)){
    N[ct,] <- Ds[dp[i]+1, dq[i]+1, ]
    ct <- ct+1
  }
  input.data$geno[id1,-c(1:2)] <- N
  # if the sum of the excluded classes is greater than zero
  # and smaller than 0.3, assign zero to those classes and normalize the vector
  input.data$geno[,-c(1:2)][R > 0] <- 0
  input.data$geno[,-c(1:2)] <- sweep(input.data$geno[,-c(1:2)], 1, rowSums(input.data$geno[,-c(1:2)]), FUN = "/")
  if(is.null(prob.thres))
    prob.thres <- input.data$prob.thres
  geno.dose <- dist_prob_to_class(geno = input.data$geno, prob.thres = prob.thres)
  if(geno.dose$flag)
  {
    input.data$geno <- geno.dose$geno
    input.data$geno.dose <- geno.dose$geno.dose
  } else {
    input.data$geno.dose <- geno.dose$geno.dose
  }
  input.data$geno.dose[is.na(input.data$geno.dose)] <- ploidy + 1
  input.data$n.ind <- ncol(input.data$geno.dose)
  input.data$ind.names <- colnames(input.data$geno.dose)
  return(input.data)
}

#' Filter missing genotypes
#'
#' Excludes markers or individuals based on their proportion of missing data.
#'
#' @param input.data an object of class \code{mappoly.data}.
#'
#' @param type one of the following options:
#' \enumerate{
#'   \item \code{"marker"}: filter out markers based on their percentage of missing data (default).
#'   \item \code{"individual"}: filter out individuals based on their percentage of missing data.
#' }
#' Please notice that removing individuals with certain amount of data can change some marker parameters
#' (such as depth), and can also change the estimated genotypes for other individuals.
#' So, be careful when removing individuals.
#'
#' @param filter.thres maximum percentage of missing data (default = 0.2).
#'
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering.
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}.
#'
#' @examples
#' plot(tetra.solcap)
#' dat.filt.mrk <- filter_missing(input.data = tetra.solcap,
#'                                type = "marker",
#'                                filter.thres = 0.1,
#'                                inter = TRUE)
#' plot(dat.filt.mrk)
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing <- function(input.data,
                           type = c("marker", "individual"),
                           filter.thres = 0.2,
                           inter = TRUE)
{
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
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
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering
#' @keywords internal
filter_missing_mrk <- function(input.data, filter.thres = 0.2, inter = TRUE)
{
  op <- par(pty="s")
  on.exit(par(op))
  ANSWER <- "flag"
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
    {
      na.num <- apply(input.data$geno.dose, 1, function(x,ploidy) sum(x == ploidy+1), ploidy = input.data$ploidy)
      perc.na <- na.num/input.data$n.ind
      x <- sort(perc.na)
      plot(x,
           xlab = "markers",
           ylab = "frequency of missing data",
           col = ifelse(x <= filter.thres, 4, 2),
           pch =ifelse(x <= filter.thres, 1, 4))
      abline(h = filter.thres, lty = 2)
      f<-paste0("Filtered out: ", sum(perc.na > filter.thres))
      i<-paste0("Included: ", sum(perc.na <= filter.thres))
      legend("topleft",  c(f, i) , col = c(2,4), pch = c(4,1))
      ANSWER <- readline("Enter 'Y/n' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N")
        stop("Stop function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
        filter.thres  <- as.numeric(ANSWER)
    }
    mrks.id <- which(perc.na <= filter.thres)
    if(length(mrks.id) == input.data$n.mrk){
      return(input.data)
    }
    out.dat <- sample_data(input.data, type = "markers", selected.mrk = names(mrks.id))
    return(out.dat)
  }
  else{
    na.num <- apply(input.data$geno.dose, 1, function(x,ploidy) sum(x == ploidy+1), ploidy = input.data$ploidy)
    perc.na <- na.num/input.data$n.ind
    mrks.id <- which(perc.na <= filter.thres)
    if(length(mrks.id) == input.data$n.mrk){
      return(input.data)
    }
    out.dat <- sample_data(input.data, type = "markers", selected.mrk = names(mrks.id))
    return(out.dat)
  }
}

#' Filter individuals based on missing genotypes
#'
#' @param input.data an object of class \code{"mappoly.data"}
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering
#' @keywords internal
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing_ind <- function(input.data, filter.thres = 0.2, inter = TRUE)
{
  op <- par(pty="s")
  on.exit(par(op))
  ANSWER <- "flag"
  ind <- NULL
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
    {
      na.num <- apply(input.data$geno.dose, 2, function(x,ploidy) sum(x == ploidy+1), ploidy = input.data$ploidy)
      perc.na <- na.num/input.data$n.mrk
      x <- sort(perc.na)
      plot(x,
           xlab = "offspring",
           ylab = "frequency of missing data",
           col = ifelse(x <= filter.thres, 4, 2),
           pch =ifelse(x <= filter.thres, 1, 4));
      abline(h = filter.thres, lty = 2)
      f<-paste0("Filtered out: ", sum(perc.na > filter.thres))
      i<-paste0("Included: ", sum(perc.na <= filter.thres))
      legend("topleft",  c(f, i) , col = c(2,4), pch = c(4,1))
      ANSWER <- readline("Enter 'Y/n' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N")
        stop("You decided to stop the function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
        filter.thres  <- as.numeric(ANSWER)
    }
    ind.id <- which(perc.na <= filter.thres)
    if(length(ind.id) == input.data$n.ind){
      return(input.data)
    }
    ind <- names(ind.id)
    out.dat <- sample_data(input.data, type = "individual", selected.ind = names(ind.id))
    return(out.dat)
  }
  else{
    na.num <- apply(input.data$geno.dose, 2, function(x,ploidy) sum(x == ploidy+1), ploidy = input.data$ploidy)
    perc.na <- na.num/input.data$n.mrk
    ind.id <- which(perc.na <= filter.thres)
    if(length(ind.id) == input.data$n.ind){
      return(input.data)
    }
    ind <- names(ind.id)
    out.dat <- sample_data(input.data, type = "individual", selected.ind = names(ind.id))
    return(out.dat)
  }
}

#' Filter markers based on chi-square test
#'
#' This function filter markers based on p-values of a chi-square test.
#' The chi-square test assumes that markers follow the expected segregation
#'  patterns under Mendelian inheritance, random chromosome bivalent
#'  pairing and no double reduction.
#'
#' @param input.obj name of input object (class \code{mappoly.data})
#'
#' @param chisq.pval.thres p-value threshold used for chi-square tests
#'  (default = Bonferroni aproximation with global alpha of 0.05, i.e.,
#'  0.05/n.mrk)
#'
#' @param inter if TRUE (default), plots distorted vs. non-distorted markers
#'
#' @return An object of class \code{mappoly.chitest.seq} which contains a list with the following components:
#' \item{keep}{markers that follow Mendelian segregation pattern}
#' \item{exclude}{markers with distorted segregation}
#' \item{chisq.pval.thres}{threshold p-value used for chi-square tests}
#' \item{data.name}{input dataset used to perform the chi-square tests}
#'
#'@examples
#' mrks.chi.filt <- filter_segregation(input.obj = tetra.solcap,
#'                                     chisq.pval.thres = 0.05/tetra.solcap$n.mrk,
#'                                     inter = TRUE)
#' seq.init <- make_seq_mappoly(mrks.chi.filt)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @importFrom graphics axis
#' @export
filter_segregation <- function(input.obj, chisq.pval.thres = NULL, inter = TRUE){
  op <- par(pty="s")
  on.exit(par(op))
  if(inherits(input.obj, c("mappoly.data"))){
    chisq.val <- input.obj$chisq.pval
    n.mrk <- input.obj$n.mrk
    data.name <- as.character(sys.call())[2]
  } else if (inherits(input.obj, c("mappoly.sequence"))){
    chisq.val <- input.obj$chisq.pval[input.obj$seq.mrk.names]
    n.mrk <- length(input.obj$seq.num)
    data.name <- input.obj$data.name
  } else {
    stop(deparse(substitute(input.obj)),
         " is not an object of class 'mappoly.data' or 'mappoly.sequence'")
  }
  ##Bonferroni approx
  if(is.null(chisq.pval.thres))
    chisq.pval.thres <- 0.05/n.mrk
  ANSWER <- "flag"
  
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
    {
      x <- log10(sort(chisq.val, decreasing = TRUE))
      th <- log10(chisq.pval.thres)
      plot(x,
           xlab = "markers",
           ylab = bquote(log[10](P)),
           col = ifelse(x <= th, 2, 4),
           pch =ifelse(x <= th, 4, 1))
      abline(h = th, lty = 2)
      f<-paste0("Filtered out: ", sum(x < th))
      i<-paste0("Included: ", sum(x >= th))
      legend("bottomleft",  c(f, i) , col = c(2,4), pch = c(4,1))
      ANSWER <- readline("Enter 'Y/n' to proceed or update the p value threshold: ")
      if(substr(ANSWER, 1, 1)  ==  "n" | substr(ANSWER, 1, 1)  ==  "no" | substr(ANSWER, 1, 1)  ==  "N")
        stop("You decided to stop the function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER  != "")
        chisq.pval.thres  <- as.numeric(ANSWER)
    }
  }
  keep <- names(which(chisq.val >= chisq.pval.thres))
  exclude <- names(which(chisq.val < chisq.pval.thres))
  structure(list(keep = keep, exclude = exclude, chisq.pval.thres = chisq.pval.thres,
                 data.name = data.name),
            class = "mappoly.chitest.seq")
}

#' Filter out individuals
#'
#' This function removes individuals from the data set. Individuals can be
#' user-defined or can be accessed via interactive kinship analysis.
#'
#' @param input.data name of input object (class \code{mappoly.data})
#'
#' @param ind.to.remove individuals to be removed. If \code{NULL} it opens
#'                      an interactive graphic to proceed with the individual
#'                      selection
#' @param inter if \code{TRUE}, expects user-input to proceed with filtering
#'
#' @param type A character string specifying the procedure to be used for 
#'             detecting outlier offspring. Options include "Gmat", which 
#'             utilizes the genomic kinship matrix, and "PCA", which employs 
#'             principal component analysis on the dosage matrix.
#' 
#'   coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#'
#' @param verbose if \code{TRUE} (default), shows the filtered out individuals
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @importFrom stats prcomp
#' @export
filter_individuals <- function(input.data, 
                               ind.to.remove = NULL, 
                               inter = TRUE,
                               type = c("Gmat", "PCA"),
                               verbose = TRUE){
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  type <- match.arg(type)
  op <- par(pty="s")
  on.exit(par(op))
  D <- t(input.data$geno.dose)
  D[D == input.data$ploidy+1] <- NA
  D <- rbind(input.data$dosage.p1, input.data$dosage.p2, D)
  rownames(D)[1:2] <- c("P1", "P2")
  if(type == "Gmat"){
    G  <- AGHmatrix::Gmatrix(D, method = "VanRaden",ploidy = input.data$ploidy)
    x <- G[1,]
    y <- G[2,]
    df <- data.frame(x = x, y = y, type = c(2, 2, rep(4, length(x)-2)))
    plot(df[,1:2], col = df$type, pch = 19,
         xlab = "relationships between the offspring and P1",
         ylab = "relationships between the offspring and P2")
    abline(c(0,1), lty = 2)
    abline(c(-0.4,1), lty = 2, col = "gray")
    abline(c(0.4,1), lty = 2, col = "gray")
    legend("topright",  c("Parents", "Offspring") , col = c(2,4), pch = 19)
    if(!is.null(ind.to.remove)){
      out.data <- sample_data(input.data, selected.ind = setdiff(input.data$ind.names, ind.to.remove))
      return(out.data)
    }
    if(interactive() && inter)
    {
      #   if (!require(gatepoints))
      #     stop("Please install package 'gatepoints' to proceed")
      ANSWER <- readline("Enter 'Y/n' to proceed with interactive filtering or quit: ")
      if(substr(ANSWER, 1, 1)  ==  "y" | substr(ANSWER, 1, 1)  ==  "yes" | substr(ANSWER, 1, 1)  ==  "Y" | ANSWER  == "")
      {
        ind.to.remove <- gatepoints::fhs(df, mark = TRUE)
        ind.to.remove <- setdiff(ind.to.remove, c("P1", "P2"))
        ind.to.include <- setdiff(rownames(df)[-c(1:2)], ind.to.remove)
        if(verbose){
          cat("Removing individual(s): \n")
          print(ind.to.remove)
          cat("...\n")
        }
        if(length(ind.to.remove) == 0){
          warning("No individuals removed. Returning original data set.")
          return(input.data)
        }
        out.data <- sample_data(input.data, selected.ind = ind.to.include)
        return(out.data)
      } else{
        warning("No individuals removed. Returning original data set.")
        return(input.data)
      }
    }    
  }
  else{
    row_means <- rowMeans(D, na.rm = TRUE)
    for (i in 1:nrow(D))
      D[i, is.na(D[i, ])] <- row_means[i]
    pc <- prcomp(D)
    x <- pc$x[,"PC1"]
    y <- pc$x[,"PC2"]
    a <- diff(range(x))*0.05
    b <- diff(range(y))*0.05
    df <- data.frame(x = x, y = y, type = c(2, 2, rep(4, length(x)-2)))
    plot(df[,1:2], col = df$type, pch = 19,
         xlab = "PC1",
         ylab = "PC2", 
         xlim = c(min(x)-a, max(x)+a),
         ylim = c(min(y)-b, max(y)+b))
    points(df[1:2,1:2], col = 2, pch = 19)
    legend("bottomleft",  c("Parents", "Offspring") , col = c(2,4), pch = 19)
    if(!is.null(ind.to.remove)){
      out.data <- sample_data(input.data, selected.ind = setdiff(input.data$ind.names, ind.to.remove))
      return(out.data)
    }
    if(interactive() && inter)
    {
      ANSWER <- readline("Enter 'Y/n' to proceed with interactive filtering or quit: ")
      if(substr(ANSWER, 1, 1)  ==  "y" | substr(ANSWER, 1, 1)  ==  "yes" | substr(ANSWER, 1, 1)  ==  "Y" | ANSWER  == "")
      {
        ind.to.remove <- gatepoints::fhs(df, mark = TRUE)
        ind.to.remove <- setdiff(ind.to.remove, c("P1", "P2"))
        ind.to.include <- setdiff(rownames(df)[-c(1:2)], ind.to.remove)
        if(verbose){
          cat("Removing individual(s): \n")
          print(ind.to.remove)
          cat("...\n")
        }
        if(length(ind.to.remove) == 0){
          warning("No individuals removed. Returning original data set.")
          return(input.data)
        }
        out.data <- sample_data(input.data, selected.ind = ind.to.include)
        return(out.data)
      } else{
        warning("No individuals removed. Returning original data set.")
        return(input.data)
      }
    }    
  }
}
#'  Remove markers that do not meet a LOD criteria
#'
#'  Remove markers that do not meet a LOD and recombination fraction
#'  criteria for at least a percentage of the pairwise marker
#'  combinations. It also removes markers with strong evidence of
#'  linkage across the whole linkage group (false positive).
#'
#' \code{thresh.LOD.ph} should be set in order to only select
#'     recombination fractions that have LOD scores associated to the
#'     linkage phase configuration higher than \code{thresh_LOD_ph}
#'     when compared to the second most likely linkage phase configuration.
#'     That action usually eliminates markers that are unlinked to the
#'     set of analyzed markers.
#'
#' @param input.twopt an object of class \code{mappoly.twopt}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase configuration
#' (default = 5)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fraction
#' (default = 5)
#'
#' @param thresh.rf threshold for recombination fractions (default = 0.15)
#'
#' @param probs indicates the probability corresponding to the filtering
#' quantiles. (default = c(0.05, 1))
#'
#' @param diag.markers A window where marker pairs should be considered.
#'    If NULL (default), all markers are considered.
#'
#' @param mrk.order marker order. Only has effect if 'diag.markers' is not NULL
#'
#' @param ncpus number of parallel processes (i.e. cores) to spawn
#' (default = 1)
#'
#' @param diagnostic.plot if \code{TRUE} produces a diagnostic plot
#'
#' @param breaks number of cells for the histogram
#'
#' @return A filtered object of class \code{mappoly.sequence}.
#' See \code{\link[mappoly]{make_seq_mappoly}} for details
#'
#' @examples
#'     all.mrk <- make_seq_mappoly(hexafake, 1:20)
#'     red.mrk <- elim_redundant(all.mrk)
#'     unique.mrks <- make_seq_mappoly(red.mrk)
#'     all.pairs <- est_pairwise_rf(input.seq = unique.mrks,
#'                                ncpus = 1,
#'                                verbose = TRUE)
#'
#'     ## Full recombination fraction matrix
#'     mat.full <- rf_list_to_matrix(input.twopt = all.pairs)
#'     plot(mat.full)
#'
#'     ## Removing disruptive SNPs
#'     tpt.filt <- rf_snp_filter(all.pairs, 2, 2, 0.07, probs = c(0.15, 1))
#'     p1.filt <- make_pairs_mappoly(input.seq = tpt.filt, input.twopt = all.pairs)
#'     m1.filt <- rf_list_to_matrix(input.twopt = p1.filt)
#'     plot(mat.full, main.text = "LG1")
#'     plot(m1.filt, main.text = "LG1.filt")
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} with updates by Gabriel Gesteira, \email{gdesiqu@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_.
#'     \doi{10.1534/g3.119.400378}
#'
#' @export rf_snp_filter
#' @importFrom ggplot2 ggplot geom_histogram aes scale_fill_manual xlab ggtitle
#' @importFrom graphics hist
rf_snp_filter <- function(input.twopt,
                          thresh.LOD.ph = 5,
                          thresh.LOD.rf = 5,
                          thresh.rf = 0.15,
                          probs = c(0.05, 1),
                          diag.markers = NULL,
                          mrk.order = NULL,
                          ncpus = 1L,
                          diagnostic.plot = TRUE,
                          breaks = 100)
{
  
  input_classes <- c("mappoly.twopt", "mappoly.twopt2")
  if (!inherits(input.twopt, input_classes)) {
    stop(deparse(substitute(input.twopt)), paste0(" is not an object of class ", paste0(input_classes, collapse =  " or ")))
  }
  probs <- range(probs)
  ## Getting filtered rf matrix
  rf_mat <-  rf_list_to_matrix(input.twopt = input.twopt, thresh.LOD.ph = thresh.LOD.ph,
                               thresh.LOD.rf = thresh.LOD.rf, thresh.rf = thresh.rf,
                               ncpus = ncpus, verbose = FALSE)
  M <- rf_mat$rec.mat
  if(!is.null(mrk.order))
    M <- M[mrk.order, mrk.order]
  if(!is.null(diag.markers))
    M[abs(col(M) - row(M)) > diag.markers] <- NA
  x <- apply(M, 1, function(x) sum(!is.na(x)))
  w <- hist(x, breaks = breaks, plot = FALSE)
  th <- quantile(x, probs = probs)
  rem <- c(which(x < th[1]), which(x > th[2]))
  ids <- names(which(x >= th[1] & x <= th[2]))
  value <- type <- NULL
  if(diagnostic.plot){
    d <- rbind(data.frame(type = "original", value = x),
               data.frame(type = "filtered", value = x[ids]))
    p <- ggplot2::ggplot(d, ggplot2::aes(value)) +
      ggplot2::geom_histogram(ggplot2::aes(fill = type),
                              alpha = 0.4, position = "identity", binwidth = diff(w$mids)[1]) +
      ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      ggplot2::ggtitle( paste0("Filtering probs: [", probs[1], " : ", probs[2], "] - Non NA values by row in rf matrix - b width: ", diff(w$mids)[1])) +
      ggplot2::xlab(paste0("Non 'NA' values at LOD.ph = ", thresh.LOD.ph, ", LOD.rf = ", thresh.LOD.rf, ", and thresh.rf = ", thresh.rf))
    print(p)
  }
  ## Returning sequence object
  ch_filt <- make_seq_mappoly(input.obj = get(input.twopt$data.name, pos = 1), arg = ids, data.name = input.twopt$data.name)
  return(ch_filt)
}

#' Edit sequence ordered by reference genome positions
#' comparing to another set order
#'
#' @param input.seq object of class mappoly.sequence with alternative order (not genomic order)
#' @param invert vector of marker names to be inverted
#' @param remove vector of marker names to be removed
#'
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#'
#' @examples
#'  \donttest{
#'   dat <- filter_segregation(tetra.solcap, inter = FALSE)
#'   seq_dat <- make_seq_mappoly(dat)
#'   seq_chr <- make_seq_mappoly(seq_dat, arg = seq_dat$seq.mrk.names[which(seq_dat$chrom=="1")])
#' 
#'   tpt <- est_pairwise_rf(seq_chr)
#'   seq.filt <- rf_snp_filter(tpt, probs = c(0.05, 0.95))
#'   mat <- rf_list_to_matrix(tpt)
#'   mat2 <- make_mat_mappoly(mat, seq.filt)
#' 
#'   seq_test_mds <- mds_mappoly(mat2)
#'   seq_mds <- make_seq_mappoly(seq_test_mds)
#'   edit_seq <- edit_order(input.seq = seq_mds)
#'  }
#'  
#' @return object of class \code{mappoly.edit.order}: a list containing
#'         vector of marker names ordered according to editions (`edited_order`);
#'         vector of removed markers names (`removed`);
#'         vector of inverted markers names (`inverted`).
#'  
#' @export
edit_order <- function(input.seq, invert= NULL, remove = NULL){
  
  if (!inherits(input.seq, "mappoly.sequence")) {
    stop(deparse(substitute(input.seq)), " is not an object of class 'mappoly.sequence'")
  }
  
  get_weird <- data.frame(x = 1:length(input.seq$genome.pos),
                          y = input.seq$genome.pos)
  
  rownames(get_weird) <- input.seq$seq.mrk.names
  get_weird <- get_weird[order(get_weird$y),]
  plot(get_weird$x, get_weird$y, xlab="input sequence order", ylab = "genomic position (bp)")
  
  inverted <- removed <- vector()
  if(!is.null(invert) | !is.null(remove)){
    if(!is.null(invert)){
      inverted <- c(inverted, as.vector(invert))
      repl <- get_weird[rev(match(as.vector(invert),rownames(get_weird))),]
      get_weird[match(as.vector(invert),rownames(get_weird)),2] <- repl[,2]
      rownames(get_weird)[match(as.vector(invert), rownames(get_weird))] <- rownames(repl) 
    }
    if(!is.null(remove)){
      removed <- c(removed, as.vector(remove))
      get_weird <- get_weird[-match(remove, rownames(get_weird)),]
    }
    plot(get_weird$x, get_weird$y, xlab="input sequence order", ylab = "genomic position (bp)")
  } else {
    cat("Mark at least three points on the plot and press `Esc` to continue.")
    if(interactive()){
      ANSWER <- "Y"
      while(substr(ANSWER, 1, 1)  ==  "y" | substr(ANSWER, 1, 1)  ==  "yes" | substr(ANSWER, 1, 1)  ==  "Y" | ANSWER  == ""){
        plot(get_weird$x, get_weird$y, xlab="input sequence order", ylab = "genomic position (bp)")
        mks.to.remove <- gatepoints::fhs(get_weird, mark = TRUE)
        if(length(which(rownames(get_weird) %in% mks.to.remove)) > 0){
          ANSWER2 <- readline("Enter 'invert/remove' to proceed with the edition: ")
          if(ANSWER2 == "invert"){
            inverted <- c(inverted, as.vector(mks.to.remove))
            repl <- get_weird[rev(match(as.vector(mks.to.remove),rownames(get_weird))),]
            get_weird[match(as.vector(mks.to.remove),rownames(get_weird)),2] <- repl[,2]
            rownames(get_weird)[match(as.vector(mks.to.remove), rownames(get_weird))] <- rownames(repl) 
          } else {
            removed <- c(removed, as.vector(mks.to.remove))
            get_weird <- get_weird[-match(mks.to.remove, rownames(get_weird)),]
          }
        }
        ANSWER <- readline("Enter 'Y/n' to proceed with interactive edition or quit: ")
      }
      plot(get_weird$x, get_weird$y, xlab="input sequence order", ylab = "genomic position (bp)")
    }
  }
  
  return(structure(list(edited_order = rownames(get_weird),
                        removed = removed,
                        inverted = inverted,
                        data.name = input.seq$data.name), class = "mappoly.edit.order"))
}

#' Filter aneuploid chromosomes from progeny individuals
#'
#' @param input.data name of input object (class \code{mappoly.data})
#'
#' @param aneuploid.info data.frame with ploidy information by chromosome (columns) for each individual
#'  in progeny (rows). The chromosome and individuals names must match the ones in the file used as input
#'  in mappoly.
#'
#' @param ploidy main ploidy
#' 
#' @param rm_missing remove also genotype information from chromosomes with missing data (NA) in the aneuploid.info file
#'
#' @examples
#'      aneuploid.info <- matrix(4, nrow=tetra.solcap$n.ind, ncol = 12)
#'      set.seed(8080)
#'      aneuploid.info[sample(1:length(aneuploid.info), round((4*length(aneuploid.info))/100),0)] <- 3
#'      aneuploid.info[sample(1:length(aneuploid.info), round((4*length(aneuploid.info))/100),0)] <- 5
#'
#'      colnames(aneuploid.info) <- paste0(1:12)
#'      aneuploid.info <- cbind(inds = tetra.solcap$ind.names, aneuploid.info)
#'
#'      filt.dat <- filter_aneuploid(input.data = tetra.solcap, 
#'      aneuploid.info = aneuploid.info, ploidy = 4)
#'
#' @author Cristiane Taniguti, \email{chtaniguti@tamu.edu}
#' 
#' @return object of class \code{mappoly.data}
#'
#' @export
filter_aneuploid <- function(input.data, aneuploid.info, ploidy, rm_missing = TRUE){
  
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  
  aneu.chroms <- colnames(aneuploid.info)[-1]
  
  keep.ind <- match(input.data$ind.names, aneuploid.info[,1])
  
  if(length(keep.ind) != length(input.data$ind.names) | any(is.na(keep.ind))) 
    stop("The aneuploid information ploidy is missing at least one sample.")
  
  aneuploid.info.filt <- aneuploid.info[keep.ind,]
  
  idx.list <- apply(aneuploid.info.filt[,-1], 1, function(x) which(x != ploidy))
  names(idx.list) <- aneuploid.info.filt[,1]
  
  idx.list <- idx.list[-which(!sapply(idx.list, function(x) length(x) > 0))]
  
  if(rm_missing){
    idx.list.na <- apply(aneuploid.info.filt[,-1], 1, function(x) which(is.na(x)))
    if(length(idx.list.na) > 0){
      names(idx.list.na) <- aneuploid.info.filt[,1]
      idx.list.na <- idx.list.na[-which(!sapply(idx.list.na, function(x) length(x) > 0))]
      idx.list <- c(idx.list, idx.list.na)
    }
    cat(round((length(unlist(idx.list))/(dim(aneuploid.info.filt)[1]*(dim(aneuploid.info.filt)[2]-1)))*100,2), "% of the chromosomes x individuals are aneuploids or has missing ploidy information\n")
  } else 
    cat(round((length(unlist(idx.list))/(dim(aneuploid.info.filt)[1]*(dim(aneuploid.info.filt)[2]-1)))*100,2), "% of the chromosomes x individuals are aneuploids\n")
  
  for(i in 1:length(idx.list)){
    idx <- which(input.data$chrom %in% names(idx.list[[i]]))
    input.data$geno.dose[idx,names(idx.list[i])] <- ploidy + 1
  }
  
  return(input.data)
}

