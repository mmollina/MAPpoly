#' Recombination fraction list to matrix
#'
#' Transforms the recombination fraction list contained in an object
#' of class \code{poly.est.two.pts.pairwise} into a recombination
#' fraction matrix
#'
#' \code{thresh_LOD_ph} should be set in order to only select
#'     recombination fractions that have LOD scores associated to the
#'     linkage phase configuration higher than \code{thresh_LOD_ph}
#'     when compared to the second most likely linkage phase configuration.
#'
#' @param input.twopt an object of class \code{poly.est.two.pts.pairwise}
#'
#' @param thresh.LOD.ph LOD score threshold for linkage phase configurations (default = 0)
#'
#' @param thresh.LOD.rf LOD score threshold for recombination fractions (default = 0)
#'
#' @param thresh.rf the threshold used for recombination fraction filtering (default = 0.5)
#'
#' @param ncpus number of parallel processes (i.e. cores) to spawn (default = 1)
#' 
#' @param shared.alleles if \code{TRUE}, computes two matrices (for both parents) indicating 
#'                       the number of homologues that share alleles (default = FALSE) 
#'
#' @param verbose if \code{TRUE} (default), current progress is shown; if
#'     \code{FALSE}, no output is produced
#'
#' @param x an object of class \code{mappoly.rf.matrix}
#'
#' @param type type of matrix that should be printed. Can be one of the
#'        following: \code{"rf"}, for recombination fraction or \code{"lod"}
#'        for LOD Score
#'
#' @param ord the order in which the markers should be plotted (default = NULL)
#'
#' @param rem which markers should be removed from the heatmap (default = NULL)
#'
#' @param main.text a character string as the title of the heatmap (default = NULL)
#'
#' @param index \code{logical} should the name of the markers be printed in the 
#' diagonal of the heatmap? (default = FALSE)
#' 
#' @param fact positive integer. factor expressed as number of cells to be aggregated 
#' (default = 1, no aggregation)
#' 
#' @param ... currently ignored
#'
#' @return A list containing two matrices. The first one contains the
#'     filtered recombination fraction and the second one contains the
#'     information matrix
#'
#' @examples
#'     all.mrk<-make_seq_mappoly(hexafake, 1:20)
#'     red.mrk<-elim_redundant(all.mrk)
#'     unique.mrks<-make_seq_mappoly(red.mrk)
#'     all.pairs<-est_pairwise_rf(input.seq = unique.mrks,
#'                                ncpus = 1,
#'                                verbose=TRUE)
#'
#'     ## Full recombination fraction matrix
#'     mat.full<-rf_list_to_matrix(input.twopt=all.pairs)
#'     plot(mat.full)
#'     plot(mat.full, type = "lod")
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
#' @import fields
#' @export rf_list_to_matrix

rf_list_to_matrix <- function(input.twopt,
                              thresh.LOD.ph = 0,
                              thresh.LOD.rf = 0,
                              thresh.rf = 0.5,
                              ncpus = 1L,
                              shared.alleles = FALSE,
                              verbose = TRUE) {
  ## checking for correct object
  if (!inherits(input.twopt, "poly.est.two.pts.pairwise")) {
    stop(deparse(substitute(input.twopt)),
         " is not an object of class 'poly.est.two.pts.pairwise'")
  }
  pair_input <- input.twopt$pairwise
  marnames <- rownames(get(input.twopt$data.name, pos = 1)$geno.dose)[sort(input.twopt$seq.num)]
  marindex <- sort(input.twopt$seq.num)
  lod.mat <- rec.mat <- matrix(NA, input.twopt$n.mrk, input.twopt$n.mrk)
  if(shared.alleles)
    ShP <- ShQ <- lod.mat
  #### UPDATE: instead of recovering the order from names, provide using the object 'input.twopt'
  #seq.num.orig<-unique(sapply(strsplit(x = names(input.twopt$pairwise), split = "-"), function(x) as.numeric(x[1])))
  #seq.num.orig<-c(seq.num.orig, strsplit(tail(names(input.twopt$pairwise), n = 1), "-")[[1]][2])
  #dimnames(lod.mat) = dimnames(rec.mat) = list(seq.num.orig, seq.num.orig)
  if (ncpus > 1) {
    start <- proc.time()
    if (verbose)
      cat("INFO: Using ", ncpus, " CPUs.\n")
    cl <- parallel::makeCluster(ncpus)
    parallel::clusterExport(cl,
                            varlist = c("thresh.LOD.ph", "thresh.LOD.rf", "thresh.rf", "shared.alleles"),
                            envir = environment())
    rf.lod.mat <- parallel::parSapply(cl, pair_input, "select_rf",
                                      thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles)
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
    rf.lod.mat <- sapply(pair_input, function(x, thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles)
    {
      if(any(is.na(x)))
      {
        if(shared.alleles){return(c(NA,NA,NA,NA))} else return(c(NA,NA))
      }
      if((nrow(x) == 1 || abs(x[2 , 1]) >= thresh.LOD.ph) &&
         abs(x[1, 3]) >= thresh.LOD.rf &&
         abs(x[1, 2]) <= thresh.rf)
      {
        if(shared.alleles){
          y <- strsplit(rownames(x), "-")
          return(c(x[1,2:3], as.numeric(y[[1]][1]), as.numeric(y[[1]][2])))
        } else {
          return(x[1,2:3])          
        }
      }
      else{
        {
          if(shared.alleles){return(c(NA,NA,NA,NA))} else return(c(NA,NA))
        }
      }
    }, thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles)
  }
  rec.mat[lower.tri(rec.mat)] <- as.numeric(rf.lod.mat[1,])
  rec.mat[upper.tri(rec.mat)] <- t(rec.mat)[upper.tri(rec.mat)]
  lod.mat[lower.tri(lod.mat)] <- as.numeric(rf.lod.mat[2,])
  lod.mat[upper.tri(lod.mat)] <- t(lod.mat)[upper.tri(lod.mat)]
  dimnames(rec.mat)<-dimnames(lod.mat)<-list(marnames, marnames)
  if(shared.alleles){
    ShP[lower.tri(ShP)] <- as.numeric(rf.lod.mat[3,])
    ShP[upper.tri(ShP)] <- t(ShP)[upper.tri(ShP)]
    ShQ[lower.tri(ShQ)] <- as.numeric(rf.lod.mat[4,])
    ShQ[upper.tri(ShQ)] <- t(ShQ)[upper.tri(ShQ)]
    dimnames(ShP)<-dimnames(ShQ)<-list(marindex, marindex)
  } else{
    ShP<-ShQ<-NULL
  }
  structure(list(thresh.LOD.ph = thresh.LOD.ph,
                 thresh.LOD.rf = thresh.LOD.rf,
                 thresh.rf = thresh.rf,
                 rec.mat = rec.mat,
                 lod.mat = abs(lod.mat),
                 ShP = ShP,
                 ShQ = ShQ,
                 data.name  = input.twopt$data.name,
                 chisq.pval.thres = input.twopt$chisq.pval.thres,
                 chisq.pval = input.twopt$chisq.pval,
                 cl = class(input.twopt)),
            class = "mappoly.rf.matrix")
}

#' @rdname rf_list_to_matrix
#' @export
print.mappoly.rf.matrix <- function(x, ...) {
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
#' @export
plot.mappoly.rf.matrix <- function(x, type = c("rf", "lod"), ord = NULL, rem = NULL, 
                                   main.text = NULL, index = FALSE, fact = 1, ...){
  type<-match.arg(type)
  if(type == "rf"){
    w<-x$rec.mat
    if(!is.null(ord))
    {
      w<-w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem) == 0))
    {
      o<-which(colnames(x$rec.mat)%in%rem)
      w<-w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    if(is.null(main.text))
      main.text<-"Recombination fraction matrix"
    col.range <-
      na.omit(rev(fields::tim.colors())[1:(ceiling(128 * max(x$rec.mat, na.rm = TRUE)) + 1)])
    brks<-NULL
  } else if(type == "lod")
  {
    w<-x$lod.mat
    if(!is.null(ord))
    {
      w<-w[ord,ord]
    }
    if(!(is.null(rem) || sum(colnames(x$rec.mat)%in%rem) == 0))
    {
      o<-which(colnames(x$rec.mat)%in%rem)
      w<-w[-o,-o]
    }
    if(fact > 1)
      w <- aggregate_matrix(w, fact)
    w[w < 1e-4]<-1e-4
    w<-log10(w)
    if(is.null(main.text))
      main.text<-"log(LOD) Score matrix"
    col.range <- na.omit(fields::tim.colors()[1:(ceiling(128 * max(x$lod.mat, na.rm = TRUE)) + 1)])
    col.range <- col.range[ceiling(seq(1, length(col.range), length.out = 10))]
    brks<-seq(min(w, na.rm = TRUE), max(w, na.rm = TRUE), length.out = 11)
    brks<-round(exp(brks/log10(exp(1))),1)
  } else stop("Invalid matrix type.")
  
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


#' Select rf and lod based on thresholds
#'
#' @param void internal function to be documented
#' @keywords internal
select_rf <- function(x, thresh.LOD.ph, thresh.LOD.rf, thresh.rf, shared.alleles = FALSE)
{
  if(any(is.na(x)))
  {
    if(shared.alleles){return(c(NA,NA,NA,NA))} else return(c(NA,NA))
  }
  if((nrow(x) == 1 || abs(x[2 , 1]) >= thresh.LOD.ph) &&
     abs(x[1, 3]) >= thresh.LOD.rf &&
     abs(x[1, 2]) <= thresh.rf)
  {
    if(shared.alleles){
      y <- strsplit(rownames(x), "-")
      return(c(x[1,2:3], as.numeric(y[[1]][1]), as.numeric(y[[1]][2])))
    } else {
      return(x[1,2:3])          
    }
  }
  else{
    {
      if(shared.alleles){return(c(NA,NA,NA,NA))} else return(c(NA,NA))
    }
  }
}
