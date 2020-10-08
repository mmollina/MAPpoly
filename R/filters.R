#' Filter non-conforming classes in F1, non double reduced population.
#'
#' @param void internal function to be documented
#' @keywords internal
#' @export
filter_non_conforming_classes<-function(input.data, prob.thres = NULL)
{
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
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
  if(!is.prob.data(input.data)){
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
  geno.dose <- dist_prob_to_class(geno = input.data$geno, prob.thres = prob.thres)
  if(geno.dose$flag)
  {
    input.data$geno <- geno.dose$geno
    input.data$geno.dose <- geno.dose$geno.dose
  } else {
    input.data$geno.dose <- geno.dose$geno.dose
  }
  input.data$geno.dose[is.na(input.data$geno.dose)] <- m + 1
  input.data$n.ind <- ncol(input.data$geno.dose)
  input.data$ind.names <- colnames(input.data$geno.dose)
  return(input.data)
}

#' Filter missing genotypes
#' 
#' Excludes markers or individuals based on their proportion of missing data
#'
#' @param input.data an object of class \code{mappoly.data} 
#' 
#' @param type one of the following options: 
#' \code{'marker'}{filter out markers based on their percentage of missing data (default)}
#' \code{'individual'}{filter out individuals based on their percentage of missing data}
#' Please notice that removing individuals with certain amount of data can change some marker parameters
#' (such as depth), and can also change the estimated genotypes for other individuals.
#' So be careful when removing individuals.
#' 
#' @param filter.thres maximum percentage of missing data (default = 0.2)
#' 
#' @param inter if \code{TRUE} (default), it plots markers or individuals vs. frequency of missing data
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @examples
#' plot(tetra.solcap)
#' dat.filt.mrk <- filter_missing(input.data = tetra.solcap,
#'                                type = "marker", 
#'                                filter.thres = 0.1,
#'                                inter = TRUE)
#' plot(dat.filt.mrk)
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing<-function(input.data, 
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
#' @param inter if \code{TRUE}, plots markers vs. frequency of genotyped individuals
#' @keywords internal
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing_mrk<-function(input.data, filter.thres = 0.2, inter = TRUE)
{
  ANSWER <- "flag"
  mrk <- NULL
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
    {
      na.num<-apply(input.data$geno.dose, 1, function(x,m) sum(x==m+1), m = input.data$m)
      perc.na<-na.num/input.data$n.ind
      plot(sort(perc.na), xlab = "markers", ylab = "frequency of missing data", axes = FALSE);
      axis(1);axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(filter.thres,2), col = 4, lty = 2)
      text(x = input.data$n.mrk/2, y = filter.thres + 0.05, labels = paste0("Excluded mrks: ", sum(perc.na >= filter.thres)), adj = 0, col = "darkred")
      text(x = input.data$n.mrk/2, y = filter.thres - 0.05, labels = paste0("Included mrks: ", sum(perc.na < filter.thres)), adj = 0, col = "darkgreen")
      ANSWER <- readline("Enter 'Y/n' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "no" | substr(ANSWER, 1, 1) == "N")
          stop("You decided to stop the function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
        filter.thres  <- as.numeric(ANSWER)
    }
    rm.mrks.id<-which(perc.na > filter.thres)
    if(length(rm.mrks.id)==0){
      return(input.data)
    } 
    rm.mrks<-names(rm.mrks.id)
    if(is.prob.data(input.data))
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
  } else {
    na.num<-apply(input.data$geno.dose, 1, function(x,m) sum(x==m+1), m = input.data$m)
    perc.na<-na.num/input.data$n.ind
    rm.mrks.id<-which(perc.na > filter.thres)
    if(length(rm.mrks.id)==0) return(input.data)
    rm.mrks<-names(rm.mrks.id)
    if(is.prob.data(input.data))
      input.data$geno <-  input.data$geno %>%
      dplyr::filter(!mrk%in%rm.mrks)
    input.data$geno.dose<-input.data$geno.dose[-rm.mrks.id,]
    input.data$n.mrk <- nrow(input.data$geno.dose)
    input.data$mrk.names <- input.data$mrk.names[-rm.mrks.id]
    input.data$dosage.p <- input.data$dosage.p[-rm.mrks.id]
    input.data$dosage.q <- input.data$dosage.q[-rm.mrks.id]
    input.data$sequence <- input.data$sequence[-rm.mrks.id]
    input.data$sequence.pos <- input.data$sequence.pos[-rm.mrks.id]
    input.data$seq.ref <- input.data$seq.ref[-rm.mrks.id]
    input.data$seq.alt <- input.data$seq.alt[-rm.mrks.id]
    input.data$all.mrk.depth <- input.data$all.mrk.depth[-rm.mrks.id]
    input.data$kept <- intersect(input.data$mrk.names, input.data$kept)
    input.data$elim.correspondence <- input.data$elim.correspondence[input.data$elim.correspondence$kept%in%input.data$mrk.names,]
    input.data$chisq.pval <- input.data$chisq.pval[names(input.data$chisq.pval)%in%input.data$mrk.names]
    if(!is.null(input.data$chisq.pval)) 
      input.data$chisq.pval <- input.data$chisq.pval[names(input.data$chisq.pval)%in%input.data$mrk.names]
    return(input.data)
  }
}

#' Filter individuals based on missing genotypes 
#'
#' @param input.data an object of class \code{"mappoly.data"} 
#' @param filter.thres maximum percentage of missing data
#' @param inter if \code{TRUE}, plots markers vs. frequency of genotyped individuals
#' @keywords internal
#' @importFrom magrittr "%>%"
#' @importFrom dplyr filter
#' @importFrom graphics axis
filter_missing_ind<-function(input.data, filter.thres = 0.2, inter = TRUE)
{
  ANSWER <- "flag"
  ind <- NULL
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
    {
      na.num<-apply(input.data$geno.dose, 2, function(x,m) sum(x==m+1), m = input.data$m)
      perc.na<-na.num/input.data$n.mrk
      plot(sort(perc.na), xlab = "individuals", ylab = "frequency of missing data", axes = FALSE);
      axis(1);axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(filter.thres,2), col = 4, lty = 2)
      text(x = input.data$n.ind/2, y = filter.thres + 0.05, labels = paste0("Excluded individuals: ", sum(perc.na >= filter.thres)), adj = 0, col = "darkred")
      text(x = input.data$n.ind/2, y = filter.thres - 0.05, labels = paste0("Included individuals: ", sum(perc.na < filter.thres)), adj = 0, col = "darkgreen")
      ANSWER <- readline("Enter 'Y/n' to proceed or update the filter threshold: ")
      if(substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "no" | substr(ANSWER, 1, 1) == "N")
          stop("You decided to stop the function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
        filter.thres  <- as.numeric(ANSWER)
    }
    rm.ind.id<-which(perc.na > filter.thres)
    if(length(rm.ind.id)==0){
      return(input.data)
    } 
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

#' Filter markers based on chi-square test
#'
#' This function filter markers based on p-values of a chi-square test. 
#' The chi-square test assumes that markers follow the expected segregation
#'  patterns under Mendelian inheritance, random chromosome bivalent 
#'  pairing and no double reduction.
#'
#' @param input.data name of input object (class \code{mappoly.data})
#' 
#' @param chisq.pval.thres p-value threshold used for chi-square tests (default = 10e-05)
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
#' mrks.chi.filt <- filter_segregation(input.data = tetra.solcap,
#'                                     chisq.pval.thres = 0.05/tetra.solcap$n.mrk,
#'                                     inter = TRUE)
#' seq.init<-make_seq_mappoly(mrks.chi.filt)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' 
#' @importFrom graphics axis
#' @export
filter_segregation<-function(input.data, chisq.pval.thres = 10e-5, inter = TRUE){
  ANSWER <- "flag"
  if (!inherits(input.data, "mappoly.data")) {
    stop(deparse(substitute(input.data)), " is not an object of class 'mappoly.data'")
  }
  if(interactive() && inter)
  {
    while(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
    {
      plot(log10(sort(input.data$chisq.pval, decreasing = TRUE)), xlab = "markers", ylab = "log10(p.val)", axes=F)
      axis(1); axis(2)
      lines(x = c(0, input.data$n.mrk), y = rep(log10(chisq.pval.thres),2), col = 4, lty = 2)
      text(x = input.data$n.mrk/2, y = 5, labels = paste0("Included mrks: ", sum(input.data$chisq.pval >= chisq.pval.thres)), adj = .5, col = "darkgreen")
      text(x = input.data$n.mrk/2, y = log10(chisq.pval.thres) - 5, labels = paste0("Excluded mrks: ", sum(input.data$chisq.pval < chisq.pval.thres)), adj = .5, col = "darkred")
      ANSWER <- readline("Enter 'Y/n' to proceed or update the p value threshold: ")
      if(substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "no" | substr(ANSWER, 1, 1) == "N")
          stop("You decided to stop the function.")
      if(substr(ANSWER, 1, 1) != "y" && substr(ANSWER, 1, 1) != "yes" && substr(ANSWER, 1, 1) != "Y" && ANSWER !="")
        chisq.pval.thres  <- as.numeric(ANSWER)
    }
  }
  keep<-names(which(input.data$chisq.pval >= chisq.pval.thres))
  exclude<-names(which(input.data$chisq.pval < chisq.pval.thres))
  structure(list(keep = keep, exclude = exclude, chisq.pval.thres = chisq.pval.thres, data.name = as.character(sys.call())[2]), class = "mappoly.chitest.seq")
}
