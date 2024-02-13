#' Create a Sequence of Markers
#'
#' Constructs a sequence of markers based on an object belonging to various specified classes. This
#' function is versatile, supporting multiple input types and configurations for generating marker sequences.
#'
#' @param input.obj An object belonging to one of the specified classes: \code{mappoly.data},
#' \code{mappoly.map}, \code{mappoly.sequence}, \code{mappoly.group}, \code{mappoly.unique.seq},
#' \code{mappoly.pcmap}, \code{mappoly.pcmap3d}, \code{mappoly.geno.ord}, or \code{mappoly.edit.order}.
#'
#' @param arg Specifies the markers to include in the sequence, accepting several formats: a string 'all' for all
#' markers; a string or vector of strings 'seqx' where x is the sequence number (0 for unassigned markers); a
#' vector of integers indicating specific markers; or a vector of integers representing linkage group numbers if
#' \code{input.obj} is of class \code{mappoly.group}. For certain classes (\code{mappoly.pcmap}, \code{mappoly.pcmap3d},
#' \code{mappoly.unique.seq}, or \code{mappoly.geno.ord}), \code{arg} can be \code{NULL}.
#'
#' @param data.name Name of the \code{mappoly.data} class object.
#'
#' @param info.parent Selection criteria based on parental information: \code{'all'} for all dosage combinations,
#' \code{'P1'} for markers informative in parent 1, or \code{'P2'} for markers informative in parent 2. Default
#' is \code{'all'}.
#'
#' @param genomic.info Optional and applicable only to \code{mappoly.group} objects. Specifies the use of genomic
#' information in sequence creation. With \code{NULL} (default), all markers defined by the grouping function are
#' included. Numeric values indicate the use of specific sequences from genomic information, aiming to match the
#' maximum number of markers with the group. Supports single values or vectors for multiple sequence consideration.
#'
#' @param x An object of class \code{mappoly.sequence}.
#'
#' @param ... Currently ignored.
#' 
#' @return Returns an object of class `mappoly.sequence`, comprising:
#'   \item{"seq.num"}{Ordered vector of marker indices according to the input.}
#'   \item{"seq.phases"}{List of linkage phases between markers; -1 for undefined phases.}
#'   \item{"seq.rf"}{Vector of recombination frequencies; -1 for not estimated frequencies.}
#'   \item{"loglike"}{Log-likelihood of the linkage map.}
#'   \item{"data.name"}{Name of the `mappoly.data` object with raw data.}
#'   \item{"twopt"}{Name of the `mappoly.twopt` object with 2-point analyses; -1 if not computed.}
#'
#' @examples
#' all.mrk <- make_seq_mappoly(hexafake, 'all')
#' seq1.mrk <- make_seq_mappoly(hexafake, 'seq1')
#' plot(seq1.mrk)
#' some.mrk.pos <- c(1,4,28,32,45)
#' some.mrk.1 <- make_seq_mappoly(hexafake, some.mrk.pos)
#' plot(some.mrk.1)
#'
#' @author Marcelo Mollinari \email{mmollin@ncsu.edu}, with modifications by Gabriel Gesteira 
#' \email{gdesiqu@ncsu.edu}
#'
#' @references
#' Mollinari, M., and Garcia, A. A. F. (2019). Linkage analysis and haplotype phasing in experimental
#' autopolyploid populations with high ploidy level using hidden Markov models. _G3: Genes|Genomes|Genetics_,
#' \doi{10.1534/g3.119.400378}.
#'
#' @export

make_seq_mappoly <- function(input.obj, 
                             arg = NULL, 
                             data.name = NULL, 
                             info.parent = c('all', 'p1', 'p2'), 
                             genomic.info = NULL) {
  ## checking for correct object
  input_classes <- c("mappoly.data", "mappoly.map", "mappoly.sequence", 
                     "mappoly.unique.seq", "mappoly.pcmap", "mappoly.pcmap3d", 
                     "mappoly.group", "mappoly.chitest.seq", "mappoly.geno.ord",
                     "mappoly.edit.order")
  if (!inherits(input.obj, input_classes)) {
    stop("invalid input object.", call. = FALSE)
  }
  ## if input object is a map, call 'make_seq_mappoly' recursively
  info.parent <- match.arg(info.parent)
  if(inherits (input.obj, "mappoly.map"))
    return(make_seq_mappoly(get(input.obj$info$data.name, pos = 1), 
                            arg = input.obj$info$mrk.names,
                            info.parent = info.parent,
                            data.name = input.obj$info$data.name))
  if(inherits (input.obj, "mappoly.sequence")){
    if(is.null(arg))
      arg = input.obj$seq.mrk.names
    if(!is.character(arg))
      stop("provide marker names when using 'mappoly.sequence' as input object.", 
           call. = FALSE)
    return(make_seq_mappoly(get(input.obj$data.name, pos = 1), 
                            arg = arg,
                            info.parent = info.parent,
                            data.name = input.obj$data.name))
  }
  ## checking for argument to make a sequence
  if (is.null(arg) && !inherits(input.obj, "mappoly.chitest.seq") && 
      !inherits(input.obj, "mappoly.unique.seq") && 
      !inherits(input.obj, "mappoly.pcmap") && 
      !inherits(input.obj, "mappoly.pcmap3d") && 
      !inherits(input.obj, "mappoly.geno.ord") &&
      !inherits(input.obj, "mappoly.edit.order")) {
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
    chisq.pval <- input.obj$chisq.pval
    chisq.pval.thres <- NULL
    ## gathering sequence data
    chrom <- genome.pos <- NULL
    if (any(!is.na(input.obj$chrom)))
      chrom <- input.obj$chrom
    if (any(!is.na(input.obj$genome.pos)))
      genome.pos <- input.obj$genome.pos
    
    ## make sequence with all markers
    if (length(arg)  ==  1 && arg  ==  "all")
    {
      if (realkeep) {
        seq.num = match(tokeep,input.obj$mrk.names)
      } else {
          seq.num = as.integer(1:input.obj$n.mrk)
          }
    }
    else if (all(is.character(arg)) && length(grep("seq", arg))  ==  length(arg))
    {
      if (length(input.obj$chrom)  ==  1 && input.obj$chrom  ==  0)
        stop("There is no chromosome information in ", deparse(substitute(input.obj)))
      seq.num1 <- as.integer(which(!is.na(match(input.obj$chrom, gsub("[^0-9]", "", arg)))))
      if(realkeep) {seq.num = intersect(seq.num1, seq.num)} else {seq.num = seq.num1}
      chrom <- input.obj$chrom[seq.num]
      if (length(input.obj$genome.pos) > 2)
        genome.pos <- input.obj$genome.pos[seq.num]
    }
    else if (all(is.character(arg)) && (length(arg)  ==  length(arg %in% input.obj$mrk.names)))
    {
      seq.num1 <- as.integer(match(arg, input.obj$mrk.names))
      if(realkeep) {
        seq.num = intersect(seq.num1, seq.num)
      } else {
        seq.num = seq.num1
        }
      chrom <- input.obj$chrom[seq.num]
      if (length(input.obj$genome.pos) > 2)
        genome.pos <- input.obj$genome.pos[seq.num]
    }
    else if (is.vector(arg) && all(is.numeric(arg)))
    {
      seq.num1 <- as.integer(arg)
      if(realkeep) {seq.num = intersect(seq.num1, seq.num)} else {seq.num = seq.num1}
      chrom <- input.obj$chrom[seq.num]
      if (length(input.obj$genome.pos) > 2)
        genome.pos <- input.obj$genome.pos[seq.num]
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
    tmp <- make_seq_mappoly(get(input.obj$data.name, pos = 1), 
                            arg = input.obj$keep,
                            info.parent = info.parent, 
                            data.name = input.obj$data.name)
    tmp$chisq.pval.thres <- input.obj$chisq.pval.thres
    tmp$chisq.pval <- get(input.obj$data.name, pos = 1)$chisq.pval[input.obj$keep]
    return(tmp)
  }
  if (inherits(input.obj, "mappoly.group"))
  {
    chisq.pval <- input.obj$chisq.pval
    chisq.pval.thres <- input.obj$chisq.pval.thres
    lgs.idx <- which(input.obj$groups.snp  %in%  arg)
    seq.num.group = as.numeric(names(input.obj$groups.snp)[lgs.idx])
    
    if (!is.null(genomic.info) && is.numeric(genomic.info)){
      seqs = colnames(input.obj$seq.vs.grouped.snp)[genomic.info]
    } else {
      if(realkeep) seq.num = intersect(seq.num.group, seq.num)
      else seq.num = seq.num.group
    }
    data.name <- input.obj$data.name
    input.obj <- get(data.name, pos = 1)
    if (!is.null(genomic.info)){
      seq.num.seq = match(input.obj$mrk.names[(input.obj$chrom %in% seqs)], input.obj$mrk.names)
      seq.num1 = intersect(seq.num.group, seq.num.seq)
      if(realkeep) seq.num = intersect(seq.num1, seq.num)
      else seq.num = seq.num1
    }
    if(!all(is.na(input.obj$chrom)) && !all(is.na(input.obj$genome.pos)))
    {
      chrom <- input.obj$chrom[seq.num]
      genome.pos <- input.obj$genome.pos[seq.num]
    }
    else
      chrom <- genome.pos <- NULL
  }
  if (inherits(input.obj, "mappoly.pcmap") | inherits(input.obj, "mappoly.pcmap3d" ))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the MDS order instead.")
    return(make_seq_mappoly(input.obj = get(input.obj$data.name, pos = 1),
                            arg = as.character(input.obj$locimap$locus),
                            info.parent = info.parent,
                            data.name = input.obj$data.name))
  }
  if (inherits(input.obj, "mappoly.geno.ord"))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the genome order instead.")
    return(make_seq_mappoly(get(input.obj$data.name, pos = 1),
                            arg = as.character(rownames(input.obj$ord)),
                            info.parent = info.parent,
                            data.name = input.obj$data.name))
  }
  if (inherits(input.obj, "mappoly.edit.order"))
  {
    if(!is.null(arg))
      warning("Ignoring argument 'arg' and using the edited sequence order instead.")
    return(make_seq_mappoly(get(input.obj$data.name, pos = 1),
                            arg = input.obj$edited_order,
                            data.name = input.obj$data.name))
  }
  dp1 <- input.obj$dosage.p1[seq.num]
  dp2 <- input.obj$dosage.p2[seq.num]
  if(info.parent == "p1")
    id <- dp2 == 0 | dp2 == input.obj$ploidy
  else if(info.parent == "p2")
    id <- dp1 == 0 | dp1 == input.obj$ploidy
  else 
    id <- seq_along(seq.num)
  structure(list(ploidy = input.obj$ploidy, 
                 seq.num = seq.num[id], 
                 seq.mrk.names = input.obj$mrk.names[seq.num][id], 
                 seq.dose.p1 = input.obj$dosage.p1[seq.num][id], 
                 seq.dose.p2 = input.obj$dosage.p2[seq.num][id],
                 seq.phases = -1, 
                 seq.rf = -1, 
                 loglike = -1, 
                 chrom = chrom[id], 
                 genome.pos = genome.pos[id], 
                 data.name = data.name, 
                 twopt = -1, 
                 chisq.pval = chisq.pval, 
                 chisq.pval.thres = chisq.pval.thres),
            class = "mappoly.sequence")
}

#' @rdname make_seq_mappoly
#' @export
print.mappoly.sequence <- function(x, ...) {
  cat("This is an object of class 'mappoly.sequence'\n")
  if (x$loglike  ==  -1) {
    cat("    ------------------------\n    Parameters not estimated\n    ------------------------\n")
  }
  n.mrk <- length(x$seq.num)
  cat("    Ploidy level:      ", x$ploidy, "\n")
  cat("    No. individuals:   ", get(x$data.name)$n.ind, "\n")
  cat("    No. markers:       ", length(x$seq.num), "\n")
  w <- table(x$chrom)
  if (all(is.null(x$chrom)) || all(is.na(x$chrom)))
    cat("\n    No. markers per sequence: not available")
  else {
    cat("\n    ----------\n    No. markers per sequence:\n")
    print(data.frame(chrom = paste0("       ", names(w)), No.mrk = as.numeric(w)), row.names = FALSE)
  }
  cat("\n    ----------\n    No. of markers per dosage in both parents:\n")
  freq <- table(paste(get(x$data.name)$dosage.p1[x$seq.num], get(x$data.name)$dosage.p2[x$seq.num], sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  print(data.frame(dP1 = paste0("    ", d.temp[, 1]), dP2 = d.temp[, 2], freq = as.numeric(freq)), row.names = FALSE)
  if (x$loglike != -1) {
    cat("\n    ---------------------------------------------\n")
    cat("\n       log-likelihood:\t", x$loglike)
    cat("\n       rec. fraction:\t", round(imf_h(x$seq.rf), 1))
    cat("\n\n")
    M <- matrix("|", n.mrk, x$ploidy * 2)
    for (i in 1:n.mrk) {
      if (all(x$seq.phases$Q[[i]] != 0))
        M[i, c(x$seq.phases$P[[i]], x$seq.phases$Q[[i]] + x$ploidy)] <- "o" else M[i, x$seq.phases$P[[i]]] <- "o"
    }
    M <- cbind(get(x$data.name)$mrk.names[x$seq.num], M)
    format(apply(M, 1, function(y) cat(c("\t", y[1], "\t", y[2:(x$ploidy + 1)], rep(" ", 4), y[(x$ploidy + 2):(x$ploidy * 2 + 1)], "\n"), collapse = "")))
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
  ploidy <- x$ploidy
  freq <- table(paste(x$seq.dose.p1, x$seq.dose.p2, sep = "-"))
  d.temp <- matrix(unlist(strsplit(names(freq), "-")), ncol = 2, byrow = TRUE)
  type <- apply(d.temp, 1, function(x,ploidy) paste0(sort(abs(abs(as.numeric(x)-(ploidy/2))-(ploidy/2))), collapse = ""), ploidy = x$ploidy)
  type.names <- names(table(type))
  mrk.dist <- as.numeric(freq)
  names(mrk.dist) <- apply(d.temp, 1 , paste, collapse = "-")
  #w <- c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696",
  #       "#737373", "#525252", "#252525", "#000000")
  #pal <- colorRampPalette(w)(length(type.names))
  layout(matrix(c(1,1,1,2,3,3,6,4,5), 3, 3), widths = c(1.2,3,.5), heights = c(1.5,4.5,.5))
  barplot(mrk.dist, las = 2, #col = pal[match(type, type.names)], 
          xlab = "Number of markers", 
          ylab = "Dosage combination", horiz = TRUE)
  pval <- x$chisq.pval[x$seq.mrk.names]
  if(is.null(x$chisq.pval))
  {
    plot(0, 0, axes = FALSE, xlab = "", ylab = "", type = "n")
    text(x = 0, y = 0, labels = "No segregation test", cex = 2)
  } else{
    par(mar = c(1,1,1,2))
    par(xaxs = "i")
    plot(log10(pval), axes = FALSE, xlab = "", ylab = "", pch = 16, 
         col = rgb(red = 0.25, green = 0.64, blue = 0.86, alpha = 0.3))
    axis(4, line = 1)
    mtext(text = bquote(log[10](P)), side = 4, line = 4, cex = .7)
  }
  par(mar = c(5,1,0,2))
  pal <- c("black", colorRampPalette(c("#D73027", "#F46D43", "#FDAE61", "#FEE090",
                                       "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1",
                                       "#4575B4"))(x$ploidy + 1))
  names(pal) <- c(-1:x$ploidy)
  M <- as.matrix(get(x$data.name, pos = 1)$geno.dose[x$seq.mrk.names,])
  M[M == x$ploidy+1] <- -1
  image(x = 1:nrow(M), z = M, axes = FALSE, xlab = "",
        col = pal[as.character(sort(unique(as.vector(M))))], useRaster = TRUE)
  mtext(text = "Markers", side = 1, line = .4)
  mtext(text = "Individuals", side = 2, line = .2)
  par(mar = c(0,0,0,0))
  plot(0:10,0:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend(0,10, 
         horiz = FALSE, 
         legend = c("missing", 0:x$ploidy),
         pch = 22,
         pt.cex = 3,
         pt.bg = pal, pt.lwd = 0,
         bty = "n", xpd = TRUE)
  par(mfrow = c(1,1))
}
