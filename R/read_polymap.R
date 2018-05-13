#' Data Input
#'
#' Reads an external data file. The format of the file is described in \code{details}
#' section. This fucntion creates an object of class \code{mappoly.data}
#'
#' The first line of the input file contains the word \code{ploidy}
#' followed by the ploidy level of the parents.  The second and third
#' lines contains the words \code{nind} and \code{nmrk} followed by
#' the number of individuals in the dataset and the total number of
#' markers, respectively. The lines number 4 and 5 should contain the
#' word \code{mrknames} and \code{indnames} followed by a sequence of
#' the names of the markers and the name of the individuals,
#' respectively. Lines 6 and 7 contain the words \code{dosageP} and
#' \code{dosageQ} followed by a sequence of numbers containing the
#' dosage of all markers in parent 1 and 2. Line 8, contain the string
#' \code{seq} fallowed by a sequence of integer numbers indicating the
#' sequence each marker belongs. It can be any 'a priori' information
#' regarding the physical distance between markers. For example, these
#' numbers could refer to chromosomes, scaffolds or even contigs, in
#' which the markers are positioned. If this information is not
#' available for a particular marker, NA should be used. If this
#' information is not available for any of the markers, the string
#' \code{seq} should be followed by a single NA. Line number 9 contain
#' the word \code{seqpos} fallowed by the physical position of the
#' markers into the sequence. The physical position can be given in
#' any unity of physical genomic distance (base pairs, for
#' instance). However, the user should be able to take decisions based
#' on these values, such as occurrence of crossing overs. Line number
#' 10 should contains the word \code{nphen} fallowed by the number of
#' phenotypic traits. Line number 11 is skipped (Usually used as
#' spacer).  The next elements are strings containing the name of the
#' phenotypic trait followed by the phenotypic values. The number of
#' lines should be the same number of phenotypic traits. Missing
#' values are represented by \code{NA}. The line number 12 +
#' \code{nphen} is skipped.  Finally, the last element is a table
#' containing the dosage for each markers (rows) for each individual
#' (columns). Missing data are represented by \code{NA}.'
#'
#' @param file.in  the name of the input file which contains the data to
#'     be read.
#'
#' @param x an object of one of the classes \code{mappoly.data}
#'
#' @param detailed if available, print the number of markers per sequence
#'
#' @param ... currently ignored
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
#'       sequence.}
#'     \item{sequence.pos}{Physical position of the markers into the
#'       sequence}
#'     \item{geno}{a data.frame containing probability distribution
#'                 for each combination of marker and offspring.}
#'     \item{n.phen}{number of phenotypic traits}
#'     \item{phen}{a matrix containing the phenotypic data. The rows
#'                 corespond to the trais and the columns correspond
#'                 to the individuals}
#' @examples
#' \dontrun{
#'     hexa.file<-system.file('extdata', 'hexafake', package = 'mappoly')
#'     hexafake<-read_geno(file.in  = hexa.file)
#'     print(hexafake, detailed = TRUE)
#'
#'     ## Same thing
#'     data(hexafake)
#'}
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2017) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_
#'
#' @export read_geno

read_geno <- function(file.in ) {
    ## get ploidy level ----------------------
    temp <- scan(file.in , what = character(), sep = " ", nlines = 1, quiet = TRUE)
    m <- na.omit(as.numeric(temp[2]))
    ## get number of individuals -------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 1, nlines = 1, quiet = TRUE)
    n.ind <- na.omit(as.numeric(temp[2]))
    ## get number of markers -----------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 2, nlines = 1, quiet = TRUE)
    n.mrk <- na.omit(as.numeric(temp[2]))
    ## get marker names ----------------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 3, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk)
        stop("\n\t\t----------------------------------
                Number of markers and length of marker
                names vector do not match.
                Please, check data.
                --------------------------------------\n")
    mrk.names <- na.omit(temp[-1])
    ## get individual names ------------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 4, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.ind)
        stop("\n\t\t----------------------------------
                Number of individuals and length of
                individual names vector do not match.
                Please, check data.
                --------------------------------------\n")
    ind.names <- na.omit(temp[-1])
    ## get dosage in parent P ----------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 5, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk)
        stop("\n\t\t-----------------------------------
                Number of markers and length of dosage
                vector for parent P do not match.
                Please, check data.
                ---------------------------------------\n")
    dosage.p <- na.omit(as.numeric(temp[-1]))
    ## get dosage in parent Q ----------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 6, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk)
        stop("\n\t\t-----------------------------------
                Number of markers and length of dosage
                vector for parent Q do not match.
                Please, check data.
                ---------------------------------------\n")
    dosage.q <- na.omit(as.numeric(temp[-1]))
    ## get sequence info ---------------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 7, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk && length(temp) - 1 > 1)
        stop("\n\t\t-------------------------------------
                Number of sequence indices and number of
                markers do not match.
                Please, check data.
                ------------------------------------------\n")
    sequence <- as.numeric(temp[-1])
    ## get sequence position info ------------
    temp <- scan(file.in , what = character(), sep = " ", skip = 8, nlines = 1, quiet = TRUE)
    temp <- temp[!temp == ""]
    if (length(temp) - 1 != n.mrk && length(temp) - 1 > 1)
        stop("\n\t\t--------------------------------------------------
                The number of sequence positions and the number of
                markers do not match\n.
                Please, check data.
                --------------------------------------------------\n")
    sequencepos <- as.numeric(temp[-1])
    ## checking for phenotypic info ----------
    temp <- scan(file.in , what = character(), sep = " ", skip = 9, quiet = TRUE)
    nphen <- na.omit(as.numeric(temp[2]))
    phen <- NULL
    if (nphen > 0) {
        phen <- read.table(file.in , skip = 11, row.names = 1, col.names = c("mrk", ind.names), colClasses = c("character", rep("numeric", n.ind)), nrows = nphen,
            comment.char = "")
    }
    cat("Reading the following data:")
    cat("\n    Ploidy level:", m)
    cat("\n    No. markers: ", n.mrk)
    cat("\n    No. individuals: ", n.ind)
    if (all(unique(nphen) != 0))
        cat("\n    This dataset contains phenotypic information.")

    if (length(sequence) > 1)
        cat("\n    This dataset contains sequence information.")
    cat("\n   ...")
    ## get genotypic info --------------------
    geno.dose <- read.table(file.in , skip = 12 + nphen)
    dim(geno.dose)
    if(nrow(geno.dose)!=length(mrk.names))
      stop("\n\t\t-------------------------------------
                Number of marker names is different from
                the number of markeres in the dataset.
                Please, check data.
                ------------------------------------------\n")
    if(ncol(geno.dose)!=length(ind.names))
      stop("\n\t\t-------------------------------------
                Number of individual names is different from
                the number of individuals in the dataset.
                Please, check data.
                ------------------------------------------\n")
    dimnames(geno.dose)<-list(mrk.names, ind.names)
    geno.dose[is.na(geno.dose)] <- m + 1
    ## returning the 'mappoly.data' object
    cat("\n    Done with reading.\n")
    structure(list(m = m,
                   n.ind = n.ind,
                   n.mrk = n.mrk,
                   ind.names = ind.names,
                   mrk.names = mrk.names,
                   dosage.p = dosage.p,
                   dosage.q = dosage.q,
                   sequence = sequence,
                   sequence.pos = sequencepos,
                   geno.dose = geno.dose,
                   nphen = nphen,
                   phen = phen),
              class = "mappoly.data")
}

#' @export
print.mappoly.data <- function(x, detailed = FALSE, ...) {
    cat("This is an object of class 'mappoly.data'\n")
    cat("    Ploidy level:   ", x$m, "\n")
    cat("    No. individuals:   ", x$n.ind, "\n")
    cat("    No. markers:       ", x$n.mrk, "\n")
    w <- table(x$sequence)
    if (length(x$sequence) <= 1)
        cat("\n    No. markers per sequence: not available") else if (detailed) {
        cat("\n    ----------\n    No. markers per sequence:\n")
        print(data.frame(seq = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
        cat("    ----------\n")
        cat(paste0("    Markers not allocated: ", sum(is.na(x$sequence))))
    } else cat("\n    This dataset contains sequence information.")
    cat("\n    ----------\n    No. of markers per dosage in both parents:\n")
    freq <- table(paste(x$dosage.p, x$dosage.q, sep = "-"))
    d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
    print(data.frame(P = paste0("    ", d.temp[, 1]), Q = d.temp[, 2], freq = as.numeric(freq)), row.names = FALSE)
    if (x$nphen != 0)
        cat("\n    This dataset contains phenotypic information.\n")
}

#' @export
#' @importFrom graphics barplot
plot.mappoly.data <- function(x, detailed = FALSE, ...)
{
  freq <- table(paste(x$dosage.p, x$dosage.q, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  mrk.dist<-matrix(as.numeric(freq), nrow = 1, dimnames = list("", apply(d.temp, 1 , paste, collapse = "-")))
  barplot(mrk.dist, las = 2)
}
