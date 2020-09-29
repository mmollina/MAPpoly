#' Create a sequence of markers
#'
#' Makes a sequence of markers based on an object of another class.
#'
#' @param input.obj an object of one of the following classes:
#'     \code{mappoly.data}, \code{mappoly.map}, \code{mappoly.group}, \code{mappoly.unique.seq},
#'     \code{mappoly.pcmap} or \code{mappoly.pcmap3d}
#'
#' @param arg can be one of the following objects: i) a string 'all',
#'     resulting in a sequence with all markers in the raw data; ii) a
#'     string or a vector of strings \code{'seqx'}, where \code{x}
#'     is the sequence (\code{x=0} indicates unassigned markers); iii) a
#'     \code{vector} of integers specifying which markers comprise the
#'     sequence; iv) an integer representing linkage group if 
#'     \code{input.object} has class \code{mappoly.group}; or v) NULL if 
#'     \code{input.object} has class \code{mappoly.pcmap}, \code{mappoly.pcmap3d} or 
#'     \code{mappoly.unique.seq}
#'
#' @param data.name name of the object of class \code{mappoly.data}
#' 
#' @param genomic.info optional argument applied for \code{mappoly.group} objects only. This argument can be \code{NULL},
#'     or can hold the numeric combination of sequences from genomic information to be used when making the sequences.
#'     When \code{genomic.info = NULL} (default), the function returns a sequence containing all markers defined 
#'     by the grouping function. When \code{genomic.info = 1}, the function returns a sequence with markers
#'     that matched the intersection between grouping function and genomic information, considering the sequence
#'     from genomic information that holds the maximum number of markers matching the group;
#'     when \code{genomic.info = c(1,2)}, the function returns a sequence with markers
#'     that matched the intersection between grouping function and genomic information, considering two sequences
#'     from genomic information that presented the maximum number of markers matching the group; and so on.
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
#'     all.mrk<-make_seq_mappoly(hexafake, 'all')
#'     seq1.mrk<-make_seq_mappoly(hexafake, 'seq1')
#'     plot(seq1.mrk)
#'     some.mrk.pos<-c(1,4,28,32,45)
#'     (some.mrk.1<-make_seq_mappoly(hexafake, some.mrk.pos))
#'     plot(some.mrk.1)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}, with modifications by Gabriel Gesteira, \email{gabrielgesteira@usp.br}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2019) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _G3: Genes, Genomes, Genetics_. 
#'     \url{https://doi.org/10.1534/g3.119.400378} 
#'
#' @export

make_seq_mappoly <- function(input.obj, arg = NULL, data.name = NULL, genomic.info = NULL) {
  ## checking for correct object
  input_classes <- c("mappoly.data", "mappoly.map", "mappoly.unique.seq", "mappoly.pcmap", "mappoly.pcmap3d", 
                     "mappoly.group", "mappoly.chitest.seq")
  if (!inherits(input.obj, input_classes)) {
    stop(deparse(substitute(input.obj)), " is not an object of class 'mappoly.data', 'mappoly.map', 
               'mappoly.chitest.seq', 'mappoly.unique.seq', 'mappoly.pcmap', 'mappoly.pcmap3d', or 'mappoly.group'")
  }
  ## if input object is a map, call 'make_seq_mappoly' recursively
  if(inherits (input.obj, "mappoly.map"))
    return(make_seq_mappoly(get(input.obj$info$data.name, pos = 1), 
                          arg = input.obj$info$mrk.names, 
                          data.name = input.obj$info$data.name))
  ## checking for argument to make a sequence
  if (is.null(arg) && !inherits(input.obj, "mappoly.chitest.seq") && !inherits(input.obj, "mappoly.unique.seq") && 
      !inherits(input.obj, "mappoly.pcmap") && !inherits(input.obj, "mappoly.pcmap3d")) {
    stop("argument 'arg' expected.")
  }
  ## Variables defined to block removing redundant markers
  realkeep = FALSE
  tokeep = FALSE
  # ## Old code to handle redundant markers
  # if ((!is.null(input.obj$kept))){
  #   tokeep = input.obj$kept
  #   realkeep = TRUE
  #   seq.num = match(tokeep,input.obj$mrk.names)
  # }
  if (inherits(input.obj, "mappoly.data"))
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
      if (realkeep) {seq.num = match(tokeep,input.obj$mrk.names)} else {seq.num = as.integer(1:input.obj$n.mrk)}
    }
    else if (all(is.character(arg)) && length(grep("seq", arg)) == length(arg))
    {
      if (length(input.obj$sequence) == 1 && input.obj$sequence == 0)
        stop("There is no sequence information in ", deparse(substitute(input.obj)))
      seq.num1 <- as.integer(which(!is.na(match(input.obj$sequence, gsub("[^0-9]", "", arg)))))
      if(realkeep) {seq.num = intersect(seq.num1, seq.num)} else {seq.num = seq.num1}
      sequence <- input.obj$sequence[seq.num]
      if (length(input.obj$sequence.pos) > 2)
        sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else if (all(is.character(arg)) && (length(arg) == length(arg %in% input.obj$mrk.names)))
    {
      seq.num1 <- as.integer(match(arg, input.obj$mrk.names))
      if(realkeep) {seq.num = intersect(seq.num1, seq.num)} else {seq.num = seq.num1}
      sequence <- input.obj$sequence[seq.num]
      if (length(input.obj$sequence.pos) > 2)
        sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else if (is.vector(arg) && all(is.numeric(arg)))
    {
      seq.num1 <- as.integer(arg)
      if(realkeep) {seq.num = intersect(seq.num1, seq.num)} else {seq.num = seq.num1}
      sequence <- input.obj$sequence[seq.num]
      if (length(input.obj$sequence.pos) > 2)
        sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else stop("Invalid argument to select markers")
    if (is.null(data.name))
      data.name <- as.character(sys.call())[2]
  }
  if (inherits(input.obj, "mappoly.unique.seq"))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the unique sequence instead.")
    return(input.obj$unique.seq)
  }
  if (inherits(input.obj, "mappoly.chitest.seq"))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using chi-square filtered markers instead.")
    tmp<-make_seq_mappoly(get(input.obj$data.name, pos = 1), arg = input.obj$keep, data.name = input.obj$data.name)
    tmp$chisq.pval.thres<-input.obj$chisq.pval.thres
    tmp$chisq.pval<-get(input.obj$data.name, pos = 1)$chisq.pval[input.obj$keep]
    return(tmp)
  }
  if (inherits(input.obj, "mappoly.group"))
  {
    chisq.pval<-input.obj$chisq.pval
    chisq.pval.thres<-input.obj$chisq.pval.thres
    if (!is.null(genomic.info) && is.numeric(genomic.info)){
      seq.num.group = as.numeric(names(which(input.obj$groups.snp == arg)))
      seqs = names(sort(input.obj$seq.vs.grouped.snp[arg,-c(ncol(input.obj$seq.vs.grouped.snp))], decreasing = T))[genomic.info]
    } else {
      seq.num1 <- as.numeric(names(which(input.obj$groups.snp == arg)))
      if(realkeep) seq.num = intersect(seq.num1, seq.num)
      else seq.num = seq.num1
    }
    data.name <- input.obj$data.name
    input.obj <- get(data.name, pos = 1)
    if (!is.null(genomic.info)){
      seq.num.seq = match(input.obj$mrk.names[(input.obj$sequence %in% seqs)], input.obj$mrk.names)
      seq.num1 = intersect(seq.num.group, seq.num.seq)
      if(realkeep) seq.num = intersect(seq.num1, seq.num)
      else seq.num = seq.num1
    }
    if(!all(is.na(input.obj$sequence)) && !all(is.na(input.obj$sequence.pos)))
    {
      sequence <- input.obj$sequence[seq.num]
      sequence.pos <- input.obj$sequence.pos[seq.num]
    }
    else
      sequence <- sequence.pos <- NULL
  }
  if (inherits(input.obj, "mappoly.pcmap") | inherits(input.obj, "mappoly.pcmap3d" ))
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
#' @importFrom graphics barplot layout mtext image legend 
#' @importFrom grDevices colorRampPalette
plot.mappoly.sequence <- function(x, ...)
{
  oldpar <- par(mar = c(5,4,1,2))
  on.exit(par(oldpar))
  m<-x$m
  freq <- table(paste(x$seq.dose.p, x$seq.dose.q, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type<-apply(d.temp, 1, function(x,m) paste0(sort(abs(abs(as.numeric(x)-(m/2))-(m/2))), collapse=""), m = x$m)
  type.names<-names(table(type))
  mrk.dist<-as.numeric(freq)
  names(mrk.dist)<-apply(d.temp, 1 , paste, collapse = "-")
  w <- c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696",
         "#737373", "#525252", "#252525", "#000000")
  pal<-colorRampPalette(w)(length(type.names))
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,4.5,.5))
  barplot(mrk.dist, las = 2, col = pal[match(type, type.names)], 
          xlab = "Number of markers", 
          ylab = "Dosage combination", horiz = TRUE)
  pval<-x$chisq.pval[x$seq.mrk.names]
  if(is.null(x$chisq.pval))
  {
    plot(0, 0, axes = FALSE, xlab = "", ylab="", type = "n")
    text(x=0, y=0, labels = "No segregation test", cex = 2)
  } else{
    par(mar=c(1,1,1,2))
    par(xaxs="i")
    plot(log10(pval), axes = FALSE, xlab = "", ylab="", pch = 16, 
         col = rgb(red=0.2, green=0.2, blue=1.0, alpha=0.2))
    axis(4, line = 1)
    mtext(text = bquote(log[10](P)), side = 4, line = 4, cex = .7)
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
  M<-as.matrix(get(x$data.name, pos = 1)$geno.dose[x$seq.mrk.names,])
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
  par(mfrow=c(1,1))
}
