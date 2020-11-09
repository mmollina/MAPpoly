#' Data Input
#'
#' Reads an external data file. The format of the file is described in the \code{Details}
#' section. This function creates an object of class \code{mappoly.data}
#'
#' The first line of the input file contains the string \code{ploidy} followed by the ploidy level of the parents.
#' The second and third lines contains the strings \code{nind} and \code{nmrk} followed by the number of individuals in 
#' the dataset and the total number of markers, respectively. Lines number 4 and 5 contain the string 
#' \code{mrknames} and \code{indnames} followed by a sequence of the names of the markers and the name of the individuals, 
#' respectively. Lines 6 and 7 contain the strings \code{dosageP} and \code{dosageQ} followed by a sequence of numbers 
#' containing the dosage of all markers in parent \code{P} and \code{Q}. Line 8, contains the string seq followed by 
#' a sequence of integer numbers indicating the sequence each marker belongs. It can be any 'a priori' 
#' information regarding the physical distance between markers. For example, these numbers could refer 
#' to chromosomes, scaffolds or even contigs, in which the markers are positioned. If this information 
#' is not available for a particular marker, NA should be used. If this information is not available for 
#' any of the markers, the string \code{seq} should be followed by a single \code{NA}. Line number 9 contains the string 
#' \code{seqpos} followed by the physical position of the markers into the sequence. The physical position can be 
#' given in any unity of physical genomic distance (base pairs, for instance). However, the user should be 
#' able to make decisions based on these values, such as the occurrence of crossing overs, etc. Line number 10 
#' should contain the string \code{nphen} followed by the number of phenotypic traits. Line number 11 is skipped 
#' (Usually used as a spacer). The next elements are strings containing the name of the phenotypic trait with no space characters
#' followed by the phenotypic values. The number of lines should be the same number of phenotypic traits. 
#' \code{NA} represents missing values. The line number 12 + \code{nphen} is skipped. Finally, the last element is a table 
#' containing the probability distribution for each combination of marker and offspring. The first two columns 
#' represent the marker and the offspring, respectively. The remaining elements represent the probability 
#' associated with each one of the possible dosages. \code{NA} represents missing data.
#' 
#'
#' @param file.in a character string with the name of (or full path to) the input file which contains the data to
#'     be read
#'     
#' @param prob.thres probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than \code{prob.thres} 
#'     are considered as missing data for the dosage calling purposes (default = 0.95)
#'     
#' @param filter.non.conforming if \code{TRUE} (default) converts data points with unexpected 
#'        genotypes (i.e. no double reduction) to 'NA'. See function \code{\link[mappoly]{segreg_poly}} 
#'        for information on expected classes and their respective frequencies.  
#'
#' @param elim.redundant logical. If \code{TRUE} (default), removes redundant markers
#' during map construction, keeping them annotated to export to the final map.
#'
#' @param verbose if \code{TRUE} (default), the current progress is shown; if
#'     \code{FALSE}, no output is produced
#' 
#' @return an object of class \code{mappoly.data} which contains a
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
#'     \item{sequence.pos}{physical position of the markers into the
#'       sequence}
#'     \item{seq.ref}{NULL (unused in this type of data)}
#'     \item{seq.alt}{NULL (unused in this type of data)}
#'     \item{all.mrk.depth}{NULL (unused in this type of data)}
#'     \item{prob.thres}{probability threshold to associate a marker call to a 
#'     dosage. Markers with maximum genotype probability smaller than 'prob.thres' 
#'     were considered as missing data in the 'geno.dose' matrix}
#'     \item{geno.dose}{a matrix containing the dosage for each markers (rows) 
#'       for each individual (columns). Missing data are represented by 
#'       \code{ploidy_level + 1}}
#'     \item{geno}{a data.frame 
#'       containing the probability distribution for each combination of
#'       marker and offspring. The first two columns represent the marker
#'       and the offspring, respectively. The remaining elements represent
#'       the probability associated to each one of the possible
#'       dosages. Missing data are converted from NA to the expected
#'       segregation ratio using function \code{\link[mappoly]{segreg_poly}}}
#'     \item{n.phen}{number of phenotypic traits}
#'     \item{phen}{a matrix containing the phenotypic data. The rows
#'                 correspond to the traits and the columns correspond
#'                 to the individuals}
#'     \item{chisq.pval}{a vector containing p-values related to the chi-squared 
#'     test of Mendelian segregation performed for all markers}
#'     \item{kept}{if elim.redundant=TRUE, holds all non-redundant markers}
#'     \item{elim.correspondence}{if elim.redundant=TRUE, holds all non-redundant markers and
#' its equivalence to the redundant ones}
#'     
#' @examples
#' \donttest{
#' #### Tetraploid Example
#' ft="https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/hexa_sample"
#' tempfl <- tempfile()
#' download.file(ft, destfile = tempfl)
#' SolCAP.dose.prob <- read_geno_prob(file.in  = tempfl)
#' print(SolCAP.dose.prob, detailed = TRUE)
#' plot(SolCAP.dose.prob)
#'}
#'
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
#' @export read_geno_prob

read_geno_prob <- function(file.in, prob.thres = 0.95, filter.non.conforming = TRUE, elim.redundant = TRUE, verbose = TRUE) {
    ## get ploidy level ----------------------
    temp <- scan(file.in, what = character(), sep = " ", nlines = 1, quiet = TRUE)
    m <- na.omit(as.numeric(temp[2]))
    ## get number of individuals -------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 1, nlines = 1, quiet = TRUE)
    n.ind <- na.omit(as.numeric(temp[2]))
    ## get number of markers -----------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 2, nlines = 1, quiet = TRUE)
    n.mrk <- na.omit(as.numeric(temp[2]))
    ## get marker names ----------------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 3, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk)
        stop("\n\t\t--------------------------------------------------
                The number of markers and the length of the marker
                names vector do not match.\n
                Please, check data.
                --------------------------------------------------\n")
    mrk.names <- na.omit(temp[-1])
    ## get individual names ------------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 4, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.ind)
        stop("\n\t\t--------------------------------------------------
                The number of individuals and the length of the
                individual names vector do not match.\n
                Please, check data.
                --------------------------------------------------\n")
    ind.names <- na.omit(temp[-1])
    ## get dosage in parent P ----------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 5, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    dosage.p <- na.omit(as.integer(temp[-1]))
    if (length(dosage.p) != n.mrk)
        stop("\n\t\t--------------------------------------------------
                The number of markers and the length of the dosage
                vector for parent P do not match.\n
                Please, check data.
                --------------------------------------------------\n")
    ## get dosage in parent Q ----------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 6, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    dosage.q <- na.omit(as.integer(temp[-1]))
    if (length(dosage.q) != n.mrk)
        stop("\n\t\t--------------------------------------------------
                The number of markers and the length of the dosage
                vector for parent Q do not match.\n
                Please, check data.
                --------------------------------------------------\n")
    ## monomorphic markers
    dp<-abs(abs(dosage.p-(m/2))-(m/2))
    dq<-abs(abs(dosage.q-(m/2))-(m/2))
    id<-dp+dq!=0
    ## get sequence info ---------------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 7, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk && length(temp) - 1 > 1)
        stop("\n\t\t--------------------------------------------------
                The number of sequence indices and the number of
                markers do not match\n.
                Please, check data.
                --------------------------------------------------\n")
    sequence <- as.integer(temp[-1])
    ## get sequence position info ------------
    temp <- scan(file.in, what = character(), sep = " ", skip = 8, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk && length(temp) - 1 > 1)
        stop("\n\t\t--------------------------------------------------
                The number of sequence positions and the number of
                markers do not match\n.
                Please, check data.
                --------------------------------------------------\n")
    sequencepos <- as.numeric(temp[-1])
    names(sequencepos) <- names(sequence) <- names(dosage.q) <- names(dosage.p) <-  mrk.names
    ## checking for phenotypic info ----------
    temp <- scan(file.in, what = character(), sep = " ", skip = 9, quiet = TRUE)
    nphen <- na.omit(as.numeric(temp[2]))
    phen <- NULL
    if (verbose && nphen > 0) {
      message("Skipping phenotype: information currently ignored")
      #phen <- read.table(file.in, skip = 11, row.names = 1, col.names = c("mrk", ind.names), colClasses = c("character", rep("numeric", n.ind)), nrows = nphen,
      #    comment.char = "")
    }
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
    geno <- read.table(file.in, skip = 12 + nphen, colClasses = c("character", "character", rep("numeric", m + 1)), nrows = n.mrk * n.ind, comment.char = "")
    colnames(geno) <- c("mrk", "ind", as.character(0:m))
    mrk <- NULL
    geno<-subset(geno, mrk%in%mrk.names[id])
    ## transforming na's in expected genotypes using Mendelian segregation
    i.na <- which(apply(geno, 1, function(x) any(is.na(x))))
    if (length(i.na) > 0) {
        m.na <- match(geno[i.na, 1], mrk.names)
        dp.na <- dosage.p[m.na]
        dq.na <- dosage.q[m.na]
        for (i in 1:length(m.na)) geno[i.na[i], -c(1, 2)] <- segreg_poly(m, dp.na[i], dq.na[i])
    }
    ## ordering data frame by individuals
    #if(all(unique(geno$ind)!=ind.names)){
      ind.names<-ind.names[order(ind.names)]
      geno<-geno[order(geno$ind),]      
    #}
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
      geno.dose[is.na(geno.dose)] <- m + 1
    }
    ## returning the 'mappoly.data' object
    if (verbose) cat("\n    Done with reading.\n")
    res<-structure(list(m = m,
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
      Ds <- array(NA, dim = c(m+1, m+1, m+1))
      for(i in 0:m)
        for(j in 0:m)
          Ds[i+1,j+1,] <- segreg_poly(m = m, dP = i, dQ = j)
      Dpop<-cbind(res$dosage.p, res$dosage.q)
      M<-t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
      dimnames(M)<-list(res$mrk.names, c(0:m))
      M<-cbind(M, res$geno.dose)
      res$chisq.pval<-apply(M, 1, mrk_chisq_test, m = m)
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
        res$elim.correspondence$seq.ref = NA
        res$elim.correspondence$seq.alt = NA
        res$elim.correspondence$all.mrk.depth = NA
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
