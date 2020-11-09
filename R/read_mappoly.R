#' Data Input
#'
#' Reads an external data file. The format of the file is described in the \code{Details}
#' section. This function creates an object of class \code{mappoly.data}
#' 
#' The first line of the input file contains the string \code{ploidy} followed by the ploidy level of the parents.
#' The second and third lines contain the strings \code{nind} and \code{nmrk} followed by the number of individuals in 
#' the dataset and the total number of markers, respectively. Lines number 4 and 5 contain the strings
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
#' containing the dosage for each marker (rows) for each individual (columns). \code{NA} represents missing values.
#'
#' @param file.in a character string with the name of (or full path to) the input file
#'  which contains the data to be read
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
#' @param x an object of class \code{mappoly.data}
#'
#' @param detailed if available, print the number of markers per sequence (default = FALSE)
#'
#' @param thresh.line position of a threshold line for p values of the segregation test (default = 10e-06)
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
#' fl1 = "https://raw.githubusercontent.com/mmollina/MAPpoly_vignettes/master/data/SolCAP_dosage"
#' tempfl <- tempfile()
#' download.file(fl1, destfile = tempfl)
#' SolCAP.dose <- read_geno(file.in  = tempfl)
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
#' @export read_geno

read_geno <- function(file.in, filter.non.conforming = TRUE, elim.redundant = TRUE, verbose = TRUE) {
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
  temp <- scan(file.in , what = character(), sep = " ", skip = 7, nlines = 1, quiet = TRUE)
  temp <- temp[!temp == ""]
  if (length(temp) - 1 != n.mrk && length(temp) - 1 > 1)
    stop("\n\t\t-------------------------------------
         Number of sequence indices and number of
         markers do not match.
         Please, check data.
         ------------------------------------------\n")
  sequence <- as.integer(temp[-1])
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
  names(sequencepos) <- names(sequence) <- names(dosage.q) <- names(dosage.p) <-  mrk.names
  ## checking for phenotypic info ----------
  temp <- scan(file.in , what = character(), sep = " ", skip = 9, quiet = TRUE)
  nphen <- na.omit(as.numeric(temp[2]))
  phen <- NULL
  if (nphen > 0) {
    phen <- read.table(file.in , skip = 11, row.names = 1, col.names = c("mrk", ind.names), colClasses = c("character", rep("numeric", n.ind)), nrows = nphen,
                       comment.char = "")
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
      cat("\n    ...\n")
  }

  ## get genotypic info --------------------
  geno.dose <- read.table(file.in , skip = 12 + nphen)
  if(nrow(geno.dose)!=length(mrk.names))
    stop("\n\t\t-------------------------------------
         Number of marker names is different from
         the number of markers in the dataset.
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
    res$dosage.p = res$dosage.p[-c(mrks.rem)]
    res$dosage.q = res$dosage.q[-c(mrks.rem)]
    res$sequence = res$sequence[-c(mrks.rem)]
    res$sequence.pos = res$sequence.pos[-c(mrks.rem)]
    res$chisq.pval = res$chisq.pval[-c(mrks.rem)]
  }
  return(res)
}

#' @rdname read_geno
#' @export
print.mappoly.data <- function(x, detailed = FALSE, ...) {
  cat("This is an object of class 'mappoly.data'\n")
  cat("    Ploidy level:                           ", x$m, "\n")
  cat("    No. individuals:                        ", x$n.ind, "\n")
  cat("    No. markers:                            ", x$n.mrk, "\n")
  miss<-round(100*sum(x$geno.dose==x$m+1)/length(as.matrix(x$geno.dose)),2)
  if(!is.null(x$kept)){
    redundant = round(100*(nrow(x$elim.correspondence)/(length(x$kept)+nrow(x$elim.correspondence))),2)
  }
  ##if no prior probabilities
  if(nrow(x$geno)==x$n.mrk){
    cat("    Missing data:                            ", miss, "%\n", sep = "")  
  } else {
    cat("    Missing data under ", x$prob.thres, " prob. threshold: ", miss, "%\n", sep = "")    
  }
  if(!is.null(x$kept)){
    cat("    Redundant markers:                       ", redundant, "%\n", sep = "")  
  }
  w <- table(x$sequence)
  if (length(w) <= 1)
    cat("\n    No. markers per sequence: not available") else if (detailed) {
      cat("\n    ----------\n    No. markers per sequence:\n")
      print(data.frame(seq = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
      cat("    ----------\n")
      cat(paste0("    Markers with no sequence information: ", sum(is.na(x$sequence))))
    } else cat("\n    This dataset contains sequence information.")
  cat("\n    ----------\n    No. of markers per dosage combination in both parents:\n")
  freq <- table(paste(x$dosage.p, x$dosage.q, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  print(data.frame(P1 = paste0("    ", d.temp[, 1]), P2 = d.temp[, 2], freq = as.numeric(freq)), row.names = FALSE)
  if (x$nphen != 0)
    cat("\n    This dataset contains phenotypic information.\n")
}

#' @rdname read_geno
#' @export
#' @importFrom graphics barplot layout mtext image legend 
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices blues9
plot.mappoly.data <- function(x, thresh.line=10e-6, ...)
{
  freq <- table(paste(x$dosage.p, x$dosage.q, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type<-apply(d.temp, 1, function(x,m) paste0(sort(abs(abs(as.numeric(x)-(m/2))-(m/2))), collapse=""), m = x$m)
  type.names<-names(table(type))
  mrk.dist<-as.numeric(freq)
  names(mrk.dist)<-apply(d.temp, 1 , paste, collapse = "-")
  w <- c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696",
         "#737373", "#525252", "#252525", "#000000")
  pal<-colorRampPalette(w)(length(type.names))
  oldpar <- par(mar = c(5,4,1,2))
  on.exit(par(oldpar))
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,2,3))
  barplot(mrk.dist, las = 2, col = pal[match(type, type.names)], 
          xlab = "Number of markers", 
          ylab = "Dosage combination", horiz = TRUE)
  if(is.null(x$chisq.pval))
  {
    plot(0, 0, axes = FALSE, xlab = "", ylab="", type = "n")
    text(x=0, y=0, labels = "No segregation test", cex = 2)
  } else{
    par(mar=c(1,1,1,2))
    par(xaxs="i")
    plot(log10(x$chisq.pval), axes = FALSE, xlab = "", ylab="", pch = 16, 
         col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2))
    axis(4, line = 1)
    mtext(text = bquote(log[10](P)), side = 4, line = 4, cex = .7)
    lines(x=c(0, x$n.mrk), y = rep(log10(thresh.line),2), col = 2, lty = 2)
  }
  par(mar=c(5,1,0,2))
  
  if(x$m == 2) {
    pal <- c("black", "#FC8D59", "#FFFFBF", "#91CF60")
  }else if(x$m == 4){
    pal <- c("black", "#D7191C", "#FDAE61", "#FFFFBF", "#A6D96A", "#1A9641")
  }else if(x$m == 6){
    pal <- c("black", "#D73027", "#FC8D59", "#FEE08B", "#FFFFBF", "#D9EF8B", "#91CF60", "#1A9850")
  }else if(x$m == 8){
    pal <- c("black", "#D73027", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#D9EF8B", "#A6D96A", "#66BD63", "#1A9850") 
  } else pal <- c("black", gg_color_hue(x$m))
  
  names(pal)<-c(-1:x$m)
  M <- as.matrix(x$geno.dose)
  M[M==x$m+1]<--1
  image(x = 1:nrow(M), z = M, axes = FALSE, xlab = "",
        col = pal[as.character(sort(unique(as.vector(M))))], useRaster = TRUE)
  mtext(text = "Markers", side = 1, line = .4)
  mtext(text = "Individuals", side = 2, line = .2)
  par(mar = c(0,0,0,0))
  plot(0:10,0:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend(0,10, 
         horiz=FALSE, 
         legend=c("missing", 0:x$m),
         pch=22,
         pt.cex = 3,
         pt.bg=pal, pt.lwd = 0,
         bty = "n", xpd=TRUE)
  if(!is.null(x$elim.correspondence)){
    par(mar = c(5,0,2,2))
    red = round(100*nrow(x$elim.correspondence)/(length(x$kept)+nrow(x$elim.correspondence)),1)
    mat = matrix(c(100-red, red), ncol = 1)
    w = barplot(mat, main="",
                xlab="", col=c(blues9[3],blues9[6]),
                axes = F, width = .5, border = NA, xlim = c(0,1)) 
    
    text(w, c((100-red)/2,   100 - red/2),  c(paste0(100 - red, " %"), paste0(red, " %")))
    mtext(text = "Unique vs. Redundant", line = -1, side = 4, cex = .8)
  }
  par(mfrow=c(1,1))
}



