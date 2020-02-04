#' Create a sequence of markers
#'
#' Makes a sequence of markers based on an object of another class.
#'
#' @param input.obj an object of one of the following classes:
#'     \code{mappoly.data}, \code{mappoly.group}, \code{mappoly.unique.seq},
#'     \code{pcmap} or \code{pcmap3d}
#'
#' @param arg can be one of the following objects: i) a string 'all',
#'     resulting in a sequence with all markers in the raw data; ii) a
#'     string or a vector of strings \code{'seqx'}, where \code{x}
#'     is the sequence (\code{x=0} indicates unassigned markers); iii) a
#'     \code{vector} of integers specifying which markers comprise the
#'     sequence; iv) an integer representing linkage group if 
#'     \code{input.object} has class \code{mappoly.group}; or v) NULL if 
#'     \code{input.object} has class \code{pcmap}, \code{pcmap3d} or 
#'     \code{mappoly.unique.seq}
#'
#' @param data.name name of the object of class \code{mappoly.data}
#'
#' @param x an object of the class \code{mappoly.sequence}
#'
#' @param ... currently ignored
#'
#' @return An object of class \code{mappoly.sequence}, which is a
#'     list containing the following components:
#'     \item{seq.num}{a \code{vector} containing the (ordered) indices
#'         of markers in the sequence, according to the input file}
#'     \item{seq.phases}{a \code{list} with the linkage phases between
#'         markers in the sequence, in corresponding positions. \code{-1}
#'         means that there are no defined linkage phases}
#'     \item{seq.rf}{a \code{vector} with the recombination
#'         frequencies between markers in the sequence. \code{-1} means
#'         that there are no estimated recombination frequencies}
#'     \item{loglike}{log-likelihood of the corresponding linkage
#'         map}
#'     \item{data.name}{name of the object of class
#'         \code{mappoly.data} with the raw data}
#'     \item{twopt}{name of the object of class \code{mappoly.twopt}
#'         with the 2-point analyses. \code{-1} means that the twopt
#'         estimates were not computed}
#'
#' @examples
#'     data(hexafake)
#'     all.mrk<-make_seq_mappoly(hexafake, 'all')
#'     seq1.mrk<-make_seq_mappoly(hexafake, 'seq1')
#'     some.mrk.pos<-c(1,4,28,32,45)
#'     (some.mrk.1<-make_seq_mappoly(hexafake, some.mrk.pos))
#'     #same thing
#'     (some.mrk.names<-hexafake$mrk.names[c(1,4,28,32,45)])
#'     some.mrk.2<-make_seq_mappoly(hexafake, some.mrk.names)
#'     identical(some.mrk.1, some.mrk.2)
#'
#'     ## Removing redundant markers and makeing a new sequence
#'     red.mrk<-elim_redundant(all.mrk)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @export

make_seq_mappoly <- function(input.obj, arg = NULL, data.name = NULL) {
  ## checking for correct object
  input_classes <- c("mappoly.data", "mappoly.unique.seq", "pcmap", "pcmap3d", "mappoly.group", "mappoly.chitest.seq")
  if (!inherits(input.obj, input_classes)) {
    stop(deparse(substitute(input.obj)), " is not an object of class 'mappoly.data',
               'mappoly.chitest.seq', 'mappoly.unique.seq', 'pcmap', 'pcmap3d', or 'mappoly.group'")
  }
  ## checking for argument to make a sequence
  if (is.null(arg) && class(input.obj) != "mappoly.chitest.seq" && class(input.obj) != "mappoly.unique.seq" && class(input.obj) != "pcmap"&& class(input.obj) != "pcmap3d") {
    stop("argument 'arg' expected.")
  }
  if (class(input.obj) == "mappoly.data")
  {
    chisq.pval<-input.obj$chisq.pval
    chisq.pval.thres<-NULL
    ## gathering sequence data
    sequence <- sequence.pos <- NULL
    if (any(!is.na(input.obj$sequence)))
      sequence <- input.obj$sequence
    if (any(!is.na(input.obj$sequence.pos)))
      sequence.pos <- input.obj$sequence.pos
    
    ## make sequence with all markers
    if (length(arg) == 1 && arg == "all")
    {
      seq.num <- as.integer(1:input.obj$n.mrk)
    }
    else if (all(is.character(arg)) && length(grep("seq", arg)) == length(arg))
    {
      if (length(input.obj$sequence) == 1 && input.obj$sequence == 0)
        stop("There is no sequence information in ", deparse(substitute(input.obj)))
      seq.num <- as.integer(which(!is.na(match(input.obj$sequence, gsub("[^0-9]", "", arg)))))
      sequence <- input.obj$sequence[seq.num]
      if (length(input.obj$sequence.pos) > 2)
        sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else if (all(is.character(arg)) && (length(arg) == length(arg %in% input.obj$mrk.names)))
    {
      seq.num <- as.integer(match(arg, input.obj$mrk.names))
      sequence <- input.obj$sequence[seq.num]
      if (length(input.obj$sequence.pos) > 2)
        sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else if (is.vector(arg) && all(is.numeric(arg)))
    {
      seq.num <- as.integer(arg)
      sequence <- input.obj$sequence[seq.num]
      if (length(input.obj$sequence.pos) > 2)
        sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else stop("Invalid argument to select markers")
    if (is.null(data.name))
      data.name <- as.character(sys.call())[2]
  }
  if (class(input.obj) == "mappoly.unique.seq")
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the unique sequence instead.")
    return(input.obj$unique.seq)
  }
  if (class(input.obj) == "mappoly.chitest.seq")
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using chi-square filtered markers instead.")
    tmp<-make_seq_mappoly(get(input.obj$data.name, pos = 1), arg = input.obj$keep, data.name = input.obj$data.name)
    tmp$chisq.pval.thres<-input.obj$chisq.pval.thres
    tmp$chisq.pval<-get(input.obj$data.name, pos = 1)$chisq.pval[input.obj$keep]
    return(tmp)
  }
  if (class(input.obj) == "mappoly.group")
  {
    chisq.pval<-input.obj$chisq.pval
    chisq.pval.thres<-input.obj$chisq.pval.thres
    seq.num <- as.numeric(names(which(input.obj$groups.snp == arg)))
    data.name <- input.obj$data.name
    input.obj <- get(data.name, pos = 1)
    if(!all(is.na(input.obj$sequence)) && !all(is.na(input.obj$sequence.pos)))
    {
      sequence <- input.obj$sequence[seq.num]
      sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else
      sequence <- sequence.pos <- NULL
  }
  if (class(input.obj) == "pcmap" | class(input.obj) == "pcmap3d" )
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the MDS order instead.")
    return(make_seq_mappoly(get(input.obj$data.name, pos = 1),
                            arg = as.character(input.obj$locimap$locus),
                            data.name = input.obj$data.name))
  }
  structure(list(m = input.obj$m, seq.num = seq.num, seq.mrk.names = input.obj$mrk.names[seq.num], 
                 seq.dose.p = input.obj$dosage.p[seq.num], seq.dose.q = input.obj$dosage.q[seq.num],
                 seq.phases = -1, seq.rf = -1, loglike = -1, sequence = sequence, sequence.pos = sequence.pos, 
                 data.name = data.name, twopt = -1, chisq.pval = chisq.pval, 
                 chisq.pval.thres =  chisq.pval.thres),
            class = "mappoly.sequence")
}

#' @rdname make_seq_mappoly
#' @keywords internal
#' @export
print.mappoly.sequence <- function(x, ...) {
  cat("This is an object of class 'mappoly.sequence'\n")
  if (x$loglike == -1) {
    cat("    ------------------------\n    Parameters not estimated\n    ------------------------\n")
  }
  n.mrk <- length(x$seq.num)
  cat("    Ploidy level:      ", x$m, "\n")
  cat("    No. individuals:   ", get(x$data.name)$n.ind, "\n")
  cat("    No. markers:       ", length(x$seq.num), "\n")
  w <- table(x$sequence)
  if (all(is.null(x$sequence)) || all(is.na(x$sequence)))
    cat("\n    No. markers per sequence: not available")
  else {
    cat("\n    ----------\n    No. markers per sequence:\n")
    print(data.frame(sequence = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
  }
  cat("\n    ----------\n    No. of markers per dosage in both parents:\n")
  freq <- table(paste(get(x$data.name)$dosage.p[x$seq.num], get(x$data.name)$dosage.q[x$seq.num], sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  print(data.frame(dP = paste0("    ", d.temp[, 1]), dQ = d.temp[, 2], freq = as.numeric(freq)), row.names = FALSE)
  if (x$loglike != -1) {
    cat("\n    ---------------------------------------------\n")
    cat("\n       log-likelihood:\t", x$loglike)
    cat("\n       rec. fraction:\t", round(imf_h(x$seq.rf), 1))
    cat("\n\n")
    M <- matrix("|", n.mrk, x$m * 2)
    for (i in 1:n.mrk) {
      if (all(x$seq.phases$Q[[i]] != 0))
        M[i, c(x$seq.phases$P[[i]], x$seq.phases$Q[[i]] + x$m)] <- "o" else M[i, x$seq.phases$P[[i]]] <- "o"
    }
    M <- cbind(get(x$data.name)$mrk.names[x$seq.num], M)
    format(apply(M, 1, function(y) cat(c("\t", y[1], "\t", y[2:(x$m + 1)], rep(" ", 4), y[(x$m + 2):(x$m * 2 + 1)], "\n"), collapse = "")))
  }
}


#' @rdname make_seq_mappoly
#' @export
#' @keywords internal
#' @importFrom graphics barplot layout mtext image legend 
#' @importFrom grDevices colorRampPalette
plot.mappoly.sequence <- function(x, ...)
{
  m<-x$m
  freq <- table(paste(x$seq.dose.p, x$seq.dose.q, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type<-apply(d.temp, 1, function(x,m) paste0(sort(abs(abs(as.numeric(x)-(m/2))-(m/2))), collapse=""), m = x$m)
  type.names<-names(table(type))
  mrk.dist<-as.numeric(freq)
  names(mrk.dist)<-apply(d.temp, 1 , paste, collapse = "-")
  pal<-colorRampPalette(RColorBrewer::brewer.pal(9,"Greys"))(length(type.names))
  op <- par(mar = c(5,4,1,2))
  layout(matrix(c(1,1,2,3,2,4), 2, 3), widths = c(1.2,3,.5), heights = c(1,2))
  barplot(mrk.dist, las = 2, col = pal[match(type, type.names)], 
          xlab = "dosage combination", 
          ylab = "number of markers", horiz = TRUE)
  pval<-x$chisq.pval[x$seq.mrk.names]
  if(is.null(x$chisq.pval))
  {
    plot(0, 0, axes = FALSE, xlab = "", ylab="", type = "n")
    text(x=0, y=0, labels = "No p values in the input object", cex = 2)
  } else{
    par(mar = c(1,1,1,3))
    plot(log10(pval), axes = FALSE, xlab = "", ylab="", pch = 16, col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2))
    axis(4)
    mtext(text = "log10(p.value)", side = 4, line = -1, cex = .7)
  }
  par(mar = c(5,1,0,4))
  pal<-c("black", RColorBrewer::brewer.pal((x$m+1),"RdYlGn"))
  names(pal)<-c(-1:x$m)
  M<-as.matrix(get(x$data.name, pos = 1)$geno.dose[x$seq.mrk.names,])
  M[M==x$m+1]<--1
  image(M, axes = FALSE,
        col = pal[as.character(sort(unique(as.vector(M))))], useRaster = TRUE)
  mtext(text = "Markers", side = 1)
  mtext(text = "Individuals", side = 2)
  par(mar = c(0,0,0,0))
  plot(0:10,0:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend(0,10, 
         horiz=FALSE, 
         legend=c("missing", 0:x$m),
         pch=22,
         pt.cex = 3,
         pt.bg=pal, pt.lwd = 0,
         bty = "n", xpd=TRUE)
  par(op)
  par(mfrow=c(1,1))
}
