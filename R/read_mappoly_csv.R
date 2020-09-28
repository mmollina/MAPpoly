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
#'  can be found in the \code{Examples} section.
#'
#' @param file.in a character string with the name of (or full path to) the input file 
#'        containing the data to be read
#'     
#' @param ploidy the ploidy level
#'     
#' @param filter.non.conforming if \code{TRUE} (default) converts data points with unexpected 
#'        genotypes (i.e. no double reduction) to 'NA'. See function \code{\link[mappoly]{segreg_poly}} 
#'        for information on expected classes and their respective frequencies.  
#'
#' @param elim.redundant logical. If \code{TRUE} (default), removes redundant markers
#'        during map construction, keeping them annotated to export to the final map.
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
#' ft="https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/tetra_solcap.csv"
#' tempfl <- tempfile()
#' download.file(ft, destfile = tempfl)
#' SolCAP.dose <- read_geno_csv(file.in  = tempfl, ploidy = 4)
#' print(SolCAP.dose, detailed = TRUE)
#' plot(SolCAP.dose)
#'}
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#' 
#'     Mollinari M., Olukolu B. A.,  Pereira G. da S., 
#'     Khan A., Gemenet D., Yencho G. C., Zeng Z-B. (2020), 
#'     Unraveling the Hexaploid Sweetpotato Inheritance 
#'     Using Ultra-Dense Multilocus Mapping, 
#'     _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400620} 
#'     
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @export read_geno_csv

read_geno_csv <- function(file.in, ploidy, filter.non.conforming = TRUE, elim.redundant = TRUE, verbose = TRUE) {
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
  if (verbose){
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
  }

  ## get genotypic info --------------------
  geno.dose <- dat[,-c(1:5)]
  dimnames(geno.dose)<-list(mrk.names, ind.names)
  geno.dose[is.na(geno.dose)] <- m + 1
  ## returning the 'mappoly.data' object
  if (verbose) cat("\n    Done with reading.\n")
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
                        seq.ref = NULL,
                        seq.alt = NULL,
                        all.mrk.depth = NULL,
                        prob.thres = NULL,
                        geno.dose = geno.dose,
                        nphen = nphen,
                        phen = phen,
                        kept = NULL,
                        elim.correspondence = NULL),
                   class = "mappoly.data")
  
  if(filter.non.conforming){
    if (verbose) cat("    Filtering non-conforming markers.\n    ...")
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
    if (verbose) cat("\n    Done with filtering.\n")
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
    res$elim.correspondence$seq.ref = NA
    res$elim.correspondence$seq.alt = NA
    res$elim.correspondence$all.mrk.depth = NA
    res$n.mrk = length(res$kept)
    res$mrk.names = res$mrk.names[-c(mrks.rem)]
    res$geno.dose = res$geno.dose[-c(mrks.rem),]
    res$dosage.p = res$dosage.p[-c(mrks.rem)]
    res$dosage.q = res$dosage.q[-c(mrks.rem)]
    res$sequence = res$sequence[-c(mrks.rem)]
    res$sequence.pos = res$sequence.pos[-c(mrks.rem)]
    res$chisq.pval = res$chisq.pval[-c(mrks.rem)]
  }
  return(res)
}
