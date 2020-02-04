#' Data Input in CSV format
#'
#' Reads an external comma-separated values (CSV) data file. The format of the file is described in the \code{Details}
#' section. This function creates an object of class \code{mappoly.data}.
#' 
#' This is an alternative and a somewhat more straightforward version of the function 
#'  \code{\link[mappoly]{read_geno}}. The input is a standard CSV file where the rows 
#'  represent the markers, except for the first row which is used as a header. 
#'  The first five columns contain the marker names, the dosage in parents 1 and 2, 
#'  the sequence information (i.e. chromosome,  scaffold, contig, etc) and the 
#'  position of the marker within the sequence. The remaining columns contain 
#'  the dosage of the full-sib population. A tetraploid example of such file 
#'  can be found in \code{ins/extdata/tetra_solcap.csv}
#'
#' @param file.in a character string with the name of (or full path to) the input file which contains the data to
#'     be read
#'     
#' @param ploidy the ploidy level
#'     
#' @param filter.non.conforming if \code{TRUE} (default) exclude samples with non 
#'     expected genotypes under random chromosome pairing and no double reduction 
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
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{n.phen}{number of phenotypic traits}
#'     \item{phen}{a matrix containing the phenotypic data. The rows
#'                 corespond to the trais and the columns correspond
#'                 to the individuals}
#' @examples
#' \dontrun{
#'     solcap.file <- system.file('extdata', 'tetra_solcap.csv', package = 'mappoly')
#'     tetra.solcap <- read_geno_csv(file.in  = solcap.file, ploidy = 4)
#'     print(tetra.solcap, detailed = TRUE)
#'
#'     ## Same thing
#'     data("tetra.solcap")
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
#' @export read_geno_csv

read_geno_csv <- function(file.in, ploidy, filter.non.conforming = TRUE) {
  m <- ploidy
  dat<-read.csv(file = file.in, header = TRUE, stringsAsFactors = FALSE)
  ## get number of individuals -------------
  n.ind <- ncol(dat) - 5
  ## get number of markers -----------------
  n.mrk <- nrow(dat)
  ## get marker names ----------------------
  mrk.names <- dat[,1]
  ## get individual names ------------------
  ind.names <- colnames(dat)[-c(1:5)]
  ## get dosage in parent P ----------------
  dosage.p <- as.integer(dat[,2])
  ## get dosage in parent Q ----------------
  dosage.q <- as.integer(dat[,3])
  ## monomorphic markers
  dp<-abs(abs(dosage.p-(m/2))-(m/2))
  dq<-abs(abs(dosage.q-(m/2))-(m/2))
  id<-dp+dq!=0
  ## get sequence info ---------------------
  sequence <- as.character(dat[,4])
  ## get sequence position info ------------
  sequencepos <- as.numeric(dat[,5])
  names(sequencepos) <- names(sequence) <- names(dosage.q) <- names(dosage.p) <-  mrk.names
  nphen <- 0
  phen <- NULL
  cat("Reading the following data:")
  cat("\n    Ploidy level:", m)
  cat("\n    No. individuals: ", n.ind)
  cat("\n    No. markers: ", n.mrk) 
  cat("\n    No. informative markers:  ", sum(id), " (", round(100*sum(id)/n.mrk,1), "%)", sep = "")
  if (all(unique(nphen) != 0))
    cat("\n    This dataset contains phenotypic information.")
  
  if (length(sequence) > 1)
    cat("\n    This dataset contains sequence information.")
  cat("\n    ...")
  ## get genotypic info --------------------
  geno.dose <- dat[,-c(1:5)]
  dimnames(geno.dose)<-list(mrk.names, ind.names)
  geno.dose[is.na(geno.dose)] <- m + 1
  ## returning the 'mappoly.data' object
  cat("\n    Done with reading.\n")
  geno.dose<-geno.dose[id,]
  res <- structure(list(m = m,
                 n.ind = n.ind,
                 n.mrk = sum(id),
                 ind.names = ind.names,
                 mrk.names = mrk.names[id],
                 dosage.p = dosage.p[id],
                 dosage.q = dosage.q[id],
                 sequence = sequence[id],
                 sequence.pos = sequencepos[id],
                 prob.thres = NULL,
                 geno.dose = geno.dose,
                 nphen = nphen,
                 phen = phen),
            class = "mappoly.data")
  
  if(filter.non.conforming){
    cat("    Filtering non-conforming markers.\n    ...")
    res<-filter_non_conforming_classes(res)
    ##Computing chi-square p.values
    Ds <- array(NA, dim = c(m+1, m+1, m+1))
    for(i in 0:m)
      for(j in 0:m)
        Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
    Dpop<-cbind(res$dosage.p, res$dosage.q)
    M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M)<-list(res$mrk.names, c(0:m))
    M<-cbind(M, res$geno.dose)
    res$chisq.pval<-apply(M, 1, mrk_chisq_test, m = m)
    cat("\n    Done with filtering.\n")
    return(res)
  }
  return(res)
}