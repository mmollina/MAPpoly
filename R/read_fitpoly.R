#' Data Input in fitPoly format
#'
#' Reads an external data file generated as output of \code{\link[fitPoly]{saveMarkerModels}}. 
#' This function creates an object of class \code{mappoly.data}.
#' 
#' @param file.in a character string with the name of (or full path to) the input file 
#' 
#' @param ploidy the ploidy level
#' 
#' @param parent1 a character string containing the name (or pattern of genotype IDs) of parent 1
#' 
#' @param parent2 a character string containing the name (or pattern of genotype IDs) of parent 2
#' 
#' @param offspring a character string containing the name (or pattern of genotype IDs) of the offspring 
#'                  individuals. If \code{NULL} (default) it considers all individuals as offsprings, except 
#'                  \code{parent1} and \code{parent2}.
#'
#' @param filter.non.conforming if \code{TRUE} (default) converts data points with unexpected 
#'        genotypes (i.e. no double reduction) to 'NA'. See function \code{\link[mappoly]{segreg_poly}} 
#'        for information on expected classes and their respective frequencies.  
#'
#' @param elim.redundant logical. If \code{TRUE} (default), removes redundant markers
#'        during map construction, keeping them annotated to in order to include them in the final map.
#'        
#' @param parent.geno indicates whether to use the joint probability \code{'joint'} (default) or the 
#'         maximum probability of multiple replicates (if available) to assign dosage to parents. 
#'         If there is one observation per parent, both options will yield the same results.
#' 
#' @param thresh.parent.geno threshold probability to assign a dosage to parents. If the probability 
#'        is smaller than \code{thresh.parent.geno}, the marker is discarded.
#' 
#' @param prob.thres threshold probability to assign a dosage to offspring. If the probability 
#'        is smaller than \code{thresh.parent.geno}, the data point is converted to 'NA'.
#' 
#' @param  file.type indicates whether the characters in the input file are separated by 
#'                  'white spaces' ("table") or by commas ("csv").
#' 
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return An object of class \code{mappoly.data} which contains a
#'     list with the following components:
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
#'     \item{seq.ref}{NULL (unused in this type of data)}
#'     \item{seq.alt}{NULL (unused in this type of data)}
#'     \item{all.mrk.depth}{NULL (unused in this type of data)}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{n.phen}{number of phenotypic traits}
#'     \item{phen}{a matrix containing the phenotypic data. The rows
#'                 correspond to the traits and the columns correspond
#'                 to the individuals}
#'     \item{kept}{if elim.redundant=TRUE, holds all non-redundant markers}
#'     \item{elim.correspondence}{if elim.redundant=TRUE, holds all non-redundant markers and
#' its equivalence to the redundant ones}
#' @examples
#' \donttest{
#' #### Tetraploid Example
#' ft <- "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/fitpoly.dat"
#' tempfl <- tempfile()
#' download.file(ft, destfile = tempfl)
#' fitpoly.dat <- read_fitpoly(file.in = tempfl, ploidy = 4, 
#'                             parent1 = "P1", parent2 = "P2", 
#'                             verbose = TRUE)
#' print(fitpoly.dat, detailed = TRUE)
#' plot(fitpoly.dat)
#' plot_mrk_info(fitpoly.dat, 37)
#'}
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#' 
#'     Voorrips, R.E., Gort, G. & Vosman, B. (2011) Genotype calling 
#'     in tetraploid species from bi-allelic marker data using mixture 
#'     models. _BMC Bioinformatics_.
#'     \url{https://doi.org/10.1186/1471-2105-12-172}
#'
#' @export read_fitpoly
#' @importFrom dplyr filter group_by summarise across
#' @importFrom utils read.delim

read_fitpoly <- function(file.in, ploidy, parent1, parent2, offspring = NULL, 
                         filter.non.conforming = TRUE, elim.redundant = TRUE, 
                         parent.geno = c("joint", "max"), thresh.parent.geno = 0.95,
                         prob.thres = 0.95, file.type = c("table", "csv"), verbose = TRUE) {
  file.type <- match.arg(file.type)
  if(file.type == "table")
    dat <- read.delim(file = file.in, header = TRUE, stringsAsFactors = FALSE)
  else if(file.type == "csv")
    dat <- read.csv(file = file.in, header = TRUE, stringsAsFactors = FALSE)
  p1 <- unique(grep(pattern = parent1, dat[,"SampleName"], value = TRUE))
  p2 <- unique(grep(pattern = parent2, dat[,"SampleName"], value = TRUE))
  if(is.null(offspring)){
    offspring <- setdiff(unique(dat[,"SampleName"]), c(p1, p2))    
  } else {
    offspring <- unique(grep(pattern = offspring, dat[,"SampleName"], value = TRUE))
  }
  parent.geno <- match.arg(parent.geno)
  dat<-dat[c(2,3,5:(5 + ploidy))]
  if(parent.geno == "joint"){
    dat.p1 <- dat %>%
      filter(SampleName %in% p1) %>%
      group_by(MarkerName) %>%
      summarise(across(2:(2+ploidy), prod), .groups = 'drop') 
    dat.p1 <- data.frame(dat.p1[,-1], row.names = dat.p1$MarkerName)
    dat.p1 <- as.data.frame(t(apply(dat.p1, 1, function(x) round(x/sum(x),3))))
    genoP1 <- apply(dat.p1, 1, function(x) { 
      if(any(is.na(x))) return(NA)
      if(max(x) < thresh.parent.geno) 
        return(NA) 
      return(which.max(x)-1)
    })
    dat.p2 <- dat %>%
      filter(SampleName %in% p2) %>%
      group_by(MarkerName) %>%
      summarise(across(2:(2+ploidy), prod), .groups = 'drop') 
    dat.p2 <- data.frame(dat.p2[,-1], row.names = dat.p2$MarkerName)
    dat.p2 <- as.data.frame(t(apply(dat.p2, 1, function(x) round(x/sum(x),3))))
    genoP2 <- apply(dat.p2, 1, function(x) { 
      if(any(is.na(x))) return(NA)
      if(max(x) < thresh.parent.geno) 
        return(NA) 
      return(which.max(x)-1)
    })
  }
  if(parent.geno == "max"){
    dat.p1 <- dat %>%
      filter(SampleName %in% p1) %>%
      group_by(MarkerName) %>%
      summarise(across(2:(2+ploidy), max), .groups = 'drop') 
    dat.p1 <- data.frame(dat.p1[,-1], row.names = dat.p1$MarkerName)
    genoP1 <- apply(dat.p1, 1, function(x) { 
      if(any(is.na(x))) return(NA)
      if(max(x) < thresh.parent.geno) 
        return(NA) 
      return(which.max(x)-1)
    })
    dat.p2 <- dat %>%
      filter(SampleName %in% p2) %>%
      group_by(MarkerName) %>%
      summarise(across(2:(2+ploidy), max), .groups = 'drop') 
    dat.p2 <- data.frame(dat.p2[,-1], row.names = dat.p2$MarkerName)
    genoP2 <- apply(dat.p2, 1, function(x) { 
      if(any(is.na(x))) return(NA)
      if(max(x) < thresh.parent.geno) 
        return(NA) 
      return(which.max(x)-1)
    })
  }
  ## get marker names ----------------------
  mrk.names <- names(which(!is.na(genoP1 + genoP2)))
  ## get number of individuals -------------
  n.ind <- length(offspring)
  ## get number of markers -----------------
  n.mrk <- length(mrk.names)
  ## get individual names ------------------
  ind.names <- offspring
  ## get dosage in parent P ----------------
  dosage.p <- as.integer(genoP1[mrk.names])
  names(dosage.p)<-mrk.names
  ## get dosage in parent Q ----------------
  dosage.q <- as.integer(genoP2[mrk.names])
  names(dosage.q)<-mrk.names
  ## monomorphic markers
  dp<-abs(abs(dosage.p-(ploidy/2))-(ploidy/2))
  dq<-abs(abs(dosage.q-(ploidy/2))-(ploidy/2))
  mrk.names<-names(which(dp+dq!=0))
  dosage.p <- dosage.p[mrk.names]
  dosage.q <- dosage.q[mrk.names]
  nphen <- 0
  phen <- NULL
  if (verbose){
      cat("Reading the following data:")
      cat("\n    Ploidy level:", ploidy)
      cat("\n    No. individuals: ", n.ind)
      cat("\n    No. markers: ", n.mrk) 
      cat("\n    No. informative markers:  ", length(mrk.names), " (", round(100*length(mrk.names)/n.mrk,1), "%)", sep = "")
      cat("\n    ...")
  }

  ## get genotypic info --------------------
  MarkerName <- SampleName <- NULL
  geno <- dat %>%
    filter(SampleName %in% offspring)  %>%
    filter(MarkerName %in% mrk.names) %>%
    arrange(SampleName, MarkerName)

  colnames(geno) <- c("mrk", "ind", as.character(0:ploidy))
  ind.names <- unique(geno$ind)
  mrk.names <- unique(geno$mrk)
  dosage.p <- dosage.p[mrk.names]
  dosage.q <- dosage.q[mrk.names]
  
  ## transforming na's in expected genotypes using Mendelian segregation
  i.na <- which(apply(geno, 1, function(x) any(is.na(x))))
  if (length(i.na) > 0) {
    m.na <- match(geno[i.na, 1], mrk.names)
    dp.na <- dosage.p[m.na]
    dq.na <- dosage.q[m.na]
    for (i in 1:length(m.na)) geno[i.na[i], -c(1, 2)] <- segreg_poly(ploidy, dp.na[i], dq.na[i])
  }
  ## dosage info
  if(filter.non.conforming){
    geno.dose <- matrix(NA,1,1)      
  } else {
    geno.dose <- dist_prob_to_class(geno = geno, prob.thres = prob.thres)
    if(geno.dose$flag)
    {
      geno <- geno.dose$geno
      geno.dose <- geno.dose$geno.dose
      n.ind <- ncol(geno.dose)
      ind.names <- colnames(geno.dose)
    } else {
      geno.dose <- geno.dose$geno.dose
    }
    geno.dose[is.na(geno.dose)] <- ploidy + 1
  }
  ## returning the 'mappoly.data' object
  if (verbose) cat("\n    Done with reading.\n")
  res<-structure(list(m = ploidy,
                      n.ind = n.ind,
                      n.mrk = length(mrk.names),
                      ind.names = ind.names,
                      mrk.names = mrk.names,
                      dosage.p = dosage.p,
                      dosage.q = dosage.q,
                      sequence = rep(NA, length(mrk.names)),
                      sequence.pos = rep(NA, length(mrk.names)),
                      seq.ref = NULL,
                      seq.alt = NULL,
                      all.mrk.depth = NULL,
                      prob.thres = prob.thres,
                      geno = geno,
                      geno.dose = geno.dose,
                      nphen = nphen,
                      phen = phen,
                      chisq.pval = NULL,
                      kept = NULL,
                      elim.correspondence = NULL),
                 class = "mappoly.data")
  
  if(filter.non.conforming){
    if (verbose) cat("    Filtering non-conforming markers.\n    ...")
    res<-filter_non_conforming_classes(res)
    if (verbose) cat("\n    Performing chi-square test.\n    ...")
    ##Computing chi-square p.values
    Ds <- array(NA, dim = c(ploidy+1, ploidy+1, ploidy+1))
    for(i in 0:ploidy)
      for(j in 0:ploidy)
        Ds[i+1,j+1,] <- segreg_poly(m = ploidy, dP = i, dQ = j)
    Dpop<-cbind(res$dosage.p, res$dosage.q)
    M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M)<-list(res$mrk.names, c(0:ploidy))
    M<-cbind(M, res$geno.dose)
    #res$chisq.pval<-apply(M, 1, mrk_chisq_test, m = ploidy)
    res$chisq.pval<-apply(M, 1, mrk_chisq_test, m = ploidy)
    if (verbose) cat("\n    Done.\n")
  }
  if (elim.redundant){
    seqred = make_seq_mappoly(res, arg = 'all', data.name = res)
    redun = elim_redundant(seqred, data = res)
    if (nrow(redun$elim.correspondence) < 1) return(res)
    res$kept = redun$kept
    res$elim.correspondence = redun$elim.correspondence
    mrks.rem = match(res$elim.correspondence$elim, res$mrk.names)
    res$elim.correspondence$sequence = res$sequence[c(mrks.rem)]
    res$elim.correspondence$sequence.pos = res$sequence.pos[c(mrks.rem)]
    res$elim.correspondence$seq.ref = NULL
    res$elim.correspondence$seq.alt = NULL
    res$elim.correspondence$all.mrk.depth = NULL
    res$n.mrk = length(res$kept)
    res$mrk.names = res$mrk.names[-c(mrks.rem)]
    res$geno.dose = res$geno.dose[-c(mrks.rem),]
    res$geno = res$geno[which(res$geno$mrk %in% rownames(res$geno.dose)),]
    res$dosage.p = res$dosage.p[-c(mrks.rem)]
    res$dosage.q = res$dosage.q[-c(mrks.rem)]
    res$sequence = res$sequence[-c(mrks.rem)]
    res$sequence.pos = res$sequence.pos[-c(mrks.rem)]
    res$chisq.pval = res$chisq.pval[-c(mrks.rem)]
  }
  return(res)
}
