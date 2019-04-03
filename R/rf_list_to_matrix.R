#' Recombination fraction list to matrix
#'
#' Transforms the recombination fraction list contained in an object
#' of class \code{poly.est.two.pts.pairwise} into a recombination
#' fraction matrix
#'
#' \code{thresh_LOD_ph} should be set in order to only selects
#'     recombination fractions which have LOD scores associated to the
#'     linkage phase configuration bigger than \code{thresh_LOD_ph}
#'     for the second most likely linkage phase configuration.
#'
#' @param input.twopt an object of class \code{poly.est.two.pts.pairwise}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase
#'     configuration.
#'
#' @param thresh.LOD.rf LOD score threshold for recombination phase
#'
#' @param thresh.rf recombination fraction threshold
#'
#' @param n.clusters Number of parallel processes to spawn
#'
#' @param verbose if \code{TRUE}, current progress is shown; if
#'     \code{FALSE}, no output is produced.
#'
#' @param x an object of class mappoly.rf.matrix
#'
#' @param type type of matrix that shuold be printed. Can be one of the
#'        following: \code{"rf"}, for recombination fraction or \code{"lod"}
#'        for LOD Score
#'
#' @param ord the order in which the markes should be ploted
#'
#' @param rem which markers should be removed from the heatmap
#'
#' @param main.text title of the heatmap
#'
#' @param index should the numbers corresponding to the markers be printed in the diagonal of the heatmap?
#'
#' @param ... currently ignored
#'
#' @return A list containing two matrices. The first one contains the
#'     filtered recombination fraction and the second one contains the
#'     information matrix
#'
#' @examples
#'   \dontrun{
#'     data(hexafake)
#'     all.mrk<-make_seq_mappoly(hexafake, 'all')
#'     red.mrk<-elim_redundant(all.mrk)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     counts.web<-cache_counts_twopt(unique.mrks, get.from.web = TRUE)
#'     all.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                count.cache = counts.web,
#'                                n.clusters = 16,
#'                                verbose=TRUE)
#'
#'     ## Full recombination fraction matrix
#'     mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
#'     plot(mat.full)
#'     plot(mat.full, type = "lod")
#'
#'     ## Fintered matrix
#'     mat.filt<-rf_list_to_matrix(input.twopt=all.pairs,
#'                                 thresh.LOD.ph = 5,
#'                                 thresh.LOD.rf = 5,
#'                                 thresh.rf = 0.5,
#'                                 n.clusters = 1,
#'                                 verbose = TRUE)
#'     plot(mat.filt)
#'     plot(mat.filt, type = "lod")
#'  }
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M., and Garcia, A.  A. F. (2018) Linkage
#'     analysis and haplotype phasing in experimental autopolyploid
#'     populations with high ploidy level using hidden Markov
#'     models, _submited_. \url{https://doi.org/10.1101/415232}
#'
#' @import fields
#' @export rf_list_to_matrix

rf_list_to_matrix <-  function(input.twopt,
                               thresh.LOD.ph = 0,
                               thresh.LOD.rf = 0,
                               thresh.rf = 0.5,
                               n.clusters = 1,
                               verbose = TRUE) {
  ## checking for correct object
  input_classes <-c("poly.est.two.pts.pairwise",
                    "poly.haplo.est.two.pts.pairwise")
  if (!inherits(input.twopt, input_classes)) {
    stop(deparse(substitute(input.twopt)),
         " is not an object of class 'poly.est.two.pts.pairwise'")
  }
    pair_input <- input.twopt$pairwise
    marnames <- rownames(get(input.twopt$data.name, pos = 1)$geno.dose)[sort(input.twopt$seq.num)]
    lod.mat <- rec.mat <- matrix(NA, input.twopt$n.mrk, input.twopt$n.mrk)
    #### UPDATE: instead of recovering the order from names, provide using the object 'input.twopt'
    #seq.num.orig<-unique(sapply(strsplit(x = names(input.twopt$pairwise), split = "-"), function(x) as.numeric(x[1])))
    #seq.num.orig<-c(seq.num.orig, strsplit(tail(names(input.twopt$pairwise), n = 1), "-")[[1]][2])
    #dimnames(lod.mat) = dimnames(rec.mat) = list(seq.num.orig, seq.num.orig)
  if (n.clusters > 1) {
    start <- proc.time()
    if (verbose)
      cat("INFO: Using ", n.clusters, " CPUs.\n")
    cl <- parallel::makeCluster(n.clusters)
    parallel::clusterExport(cl,
                            varlist = c("thresh.LOD.ph", "thresh.LOD.rf", "thresh.rf"),
                            envir = environment())
    rf.lod.mat <- parallel::parSapply(cl, pair_input, "select_rf",
                                      thresh.LOD.ph, thresh.LOD.rf, thresh.rf)
    parallel::stopCluster(cl)
    end <- proc.time()
    if (verbose) {
      cat("INFO: Done with",
          length(pair_input),
          " pairs of markers \n")
      cat("INFO: Operation took:",
          round((end - start)[3],
                digits = 3),
          "seconds\n")
    }
  } else {
    if (verbose) {
      cat("INFO: Going singlemode. Using one CPU.\n")
    }
    rf.lod.mat <- sapply(pair_input, function(x, thresh.LOD.ph, thresh.LOD.rf, thresh.rf)
    {
      if(any(is.na(x)))
        return(c(NA,NA))
      if((nrow(x) == 1 || abs(x[2 , 1]) >= thresh.LOD.ph) &&
         abs(x[1, 3]) >= thresh.LOD.rf &&
         abs(x[1, 2]) <= thresh.rf)
        return(x[1,2:3])
      else return (c(NA,NA))
    }, thresh.LOD.ph, thresh.LOD.rf, thresh.rf)
  }
  rec.mat[lower.tri(rec.mat)] <- as.numeric(rf.lod.mat[1,])
  rec.mat[upper.tri(rec.mat)] <- t(rec.mat)[upper.tri(rec.mat)]
  lod.mat[lower.tri(lod.mat)] <- as.numeric(rf.lod.mat[2,])
  lod.mat[upper.tri(lod.mat)] <- t(lod.mat)[upper.tri(lod.mat)]
  dimnames(rec.mat)<-dimnames(lod.mat)<-list(marnames, marnames)
  structure(list(thresh.LOD.ph = thresh.LOD.ph,
                 thresh.LOD.rf = thresh.LOD.rf,
                 thresh.rf = thresh.rf,
                 #id = id,
                 rec.mat = rec.mat,
                 lod.mat = abs(lod.mat),
                 data.name  = input.twopt$data.name,
                 cl = class(input.twopt)),
            class = "mappoly.rf.matrix")
}
#' @rdname rf_list_to_matrix
#' @keywords internal 
#' @export
print.mappoly.rf.matrix <- function(x, ...) {
  ## checking for correct object
  if (!any(class(x) == "mappoly.rf.matrix"))
    stop(deparse(substitute(x)), " is not an object of class 'group'")

  cat("  This is an object of class 'mappoly.rf.matrix'\n\n")
  ## criteria
  cat("  Criteria used to filter markers:\n\n")
  cat(
    "      Configuration phase LOD:           ",
    x$thresh.LOD.ph,
    "\n      Recombination fraction LOD:        ",
    x$thresh.LOD.rf,
    "\n      Maximum recombination fraction:    ",
    x$thresh.rf,
    "\n"
  )
  n.mrk <- ncol(x$rec.mat)
  per.fill <-
    round(100 * sum(!is.na(x$rec.mat)) / (length(x$rec.mat)-ncol(x$rec.mat)), 1)
  ## printing summary
  cat("\n  No. markers:             ", n.mrk, "\n")
  cat("  Percentage filled:       ", per.fill, "%\n")
}

#' @rdname rf_list_to_matrix
#' @keywords internal 
#' @export
plot.mappoly.rf.matrix <- function(x, type = c("rf", "lod"), ord = NULL, rem = NULL, main.text = NULL, index = TRUE, ...)
{
  type<-match.arg(type)
  if(type == "rf"){
    w<-x$rec.mat
    if(is.null(main.text))
      main.text<-"Recombination fraction matrix"
    col.range <-
      na.omit(rev(fields::tim.colors())[1:(ceiling(128 * max(x$rec.mat, na.rm = TRUE)) + 1)])
    brks<-NULL
  } else if(type == "lod")
  {
    w<-log10(x$lod.mat)
    if(is.null(main.text))
      main.text<-"log(LOD) Score matrix"
    col.range <- na.omit(fields::tim.colors()[1:(ceiling(128 * max(x$lod.mat, na.rm = TRUE)) + 1)])
    col.range <- col.range[ceiling(seq(1, length(col.range), length.out = 10))]
    brks<-seq(min(w, na.rm = TRUE), max(w, na.rm = TRUE), length.out = 11)
    brks<-round(exp(brks/log10(exp(1))),1)
  } else stop("Invalid matrix type.")
  if(!is.null(ord))
  {
    w<-w[ord,ord]
  }
  if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem) == 0))
  {
    o<-which(colnames(x$rec.mat)%in%rem)
    w<-w[-o,-o]
  }

  fields::image.plot(
    w,
    col = col.range,
    lab.breaks = brks,
    main = main.text,
    useRaster = FALSE,
    axes = FALSE
  )
  if(ncol(w) < 100)
    ft<-.7
  else
    ft<-100/ncol(w)
  if(index)
    text(x=seq(0,1, length.out=ncol(w)), y=seq(0,1, length.out=ncol(w)),
         labels = colnames(w), cex=ft)
}


#' Select rf adn lod based on thresholds
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
select_rf <- function(x, thresh.LOD.ph, thresh.LOD.rf, thresh.rf)
{
  if(any(is.na(x)))
    return(c(NA,NA))
  if((nrow(x) == 1 || abs(x[2 , 1]) >= thresh.LOD.ph) &&
     abs(x[1, 3]) >= thresh.LOD.rf &&
     abs(x[1, 2]) <= thresh.rf)
    return(x[1,2:3])
  else return (c(NA,NA))
}
