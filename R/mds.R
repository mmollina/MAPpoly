#' Estimates loci position using Multidimensional Scaling
#'
#' Estimates loci position using Multidimensional Scaling proposed by
#' \cite{Preedy and Hackett (2016)}. The code is an adaptation from
#' the package \code{TetraploidSNPMap}, available under GNU GENERAL PUBLIC LICENSE,
#' Version 3, at
#' \url{https://github.com/BiomathematicsAndStatisticsScotland/TetraploidSNPMap}
#'
#' @param input.mat an object of class \code{mappoly.input.matrix}
#'
#' @param p integer. The smoothing parameter for the principal curve.
#'   If \code{NULL} (default) this will be done using the leave-one-out cross validation
#'
#' @param n vector of integers or strings containing loci to be omitted from the analysis
#'
#' @param ndim number of dimensions to be considered in the multidimensional scaling procedure (default = 2)
#'
#' @param weight.exponent the exponent that should be used in the LOD score values to weight the
#'        MDS procedure (default = 2)
#'
#' @param verbose if \code{TRUE} (default), display information about the analysis
#'
#' @param displaytext logical. If \code{TRUE}, display the name of the markers in the
#'   diagnostic plot
#'   
#' @param x an object of class \code{mappoly.mds}
#'
#' @param ... currently ignored
#'
#' @return A list containing:
#' \item{M}{the input distance map}
#' \item{sm}{the unconstrained wMDS results}
#' \item{pc}{the principal curve results}
#' \item{distmap}{a matrix of pairwise distances between
#' loci where the columns are in the estimated order}
#' \item{locimap}{a data frame of the loci containing the name
#' and position of each locus in order of increasing distance}
#' \item{length}{integer giving the total length of the segment}
#' \item{removed}{a vector of the names of loci removed from the analysis}
#' \item{scale}{the scaling factor from the MDS}
#' \item{locikey}{a data frame showing the number associated with each
#' locus name for interpreting the MDS configuration plot}
#' \item{confplotno}{a data frame showing locus name associated
#' with each number on the MDS configuration plots}
#'
#' @examples
#'  \dontrun{
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
#'
#'     lgs <- group_mappoly(input.mat = mat.full,
#'                          input.seq = unique.mrks,
#'                          expected.groups = 3,
#'                          inter = TRUE,
#'                          comp.mat = TRUE, #this data has physical information
#'                          verbose=TRUE)
#'     lgs
#'     plot(lgs)
#'     lg1 <- make_seq_mappoly(lgs, 1)
#'     lg2 <- make_seq_mappoly(lgs, 2)
#'     lg3 <- make_seq_mappoly(lgs, 3)
#'
#'     ##Plot matrices
#'     p1<-make_pairs_mappoly(input.seq = lg1, input.twopt = all.pairs)
#'     p2<-make_pairs_mappoly(input.seq = lg2, input.twopt = all.pairs)
#'     p3<-make_pairs_mappoly(input.seq = lg3, input.twopt = all.pairs)
#'
#'     m1<-rf_list_to_matrix(input.twopt = p1)
#'     m2<-rf_list_to_matrix(input.twopt = p2)
#'     m3<-rf_list_to_matrix(input.twopt = p3)
#'
#'     op<-par(mfrow = c(1,3), pty = "s")
#'     plot(m1, main.text = "LG1")
#'     plot(m2, main.text = "LG2")
#'     plot(m3, main.text = "LG3")
#'     par(op)
#'
#'     ## Removing disruptive SNPs
#'
#'     lg1.filt<-rf_snp_filter(p1, 5, 5, 0.15, thresh.perc = 0.05)
#'     lg2.filt<-rf_snp_filter(p2, 5, 5, 0.15, thresh.perc = 0.05)
#'     lg3.filt<-rf_snp_filter(p3, 5, 5, 0.15, thresh.perc = 0.05)
#'
#'     p1.filt<-make_pairs_mappoly(input.seq = lg1.filt, input.twopt = all.pairs)
#'     p2.filt<-make_pairs_mappoly(input.seq = lg2.filt, input.twopt = all.pairs)
#'     p3.filt<-make_pairs_mappoly(input.seq = lg3.filt, input.twopt = all.pairs)
#'
#'     m1.filt<-rf_list_to_matrix(input.twopt = p1.filt)
#'     m2.filt<-rf_list_to_matrix(input.twopt = p2.filt)
#'     m3.filt<-rf_list_to_matrix(input.twopt = p3.filt)
#'
#'     op<-par(mfrow = c(2,3), pty = "s")
#'     plot(m1, main.text = "LG1")
#'     plot(m2, main.text = "LG2")
#'     plot(m3, main.text = "LG3")
#'     plot(m1.filt, main.text = "LG1.filt")
#'     plot(m2.filt, main.text = "LG2.filt")
#'     plot(m3.filt, main.text = "LG3.filt")
#'     par(op)
#'     mds.ord1 <- mds_mappoly(input.mat = m1.filt, p = NULL, n = NULL, ndim = 2)
#'     plot(mds.ord1)
#'     mds.ord1 <- mds_mappoly(input.mat = m1.filt, p = NULL,
#'                             n = c(346, 333), ndim = 2)
#'     plot(mds.ord1)
#'     mds.ord1 <- mds_mappoly(input.mat = m1.filt,
#'                            p = NULL, n = c(346, 333, 437, 445, 275,31),
#'                            ndim = 2)
#'     plot(mds.ord1)
#'
#'     mds.ord2 <- mds_mappoly(input.mat = m2.filt, p = NULL, n = NULL, ndim = 2)
#'     plot(mds.ord2)
#'     mds.ord3 <- mds_mappoly(input.mat = m3.filt, p = NULL, n = NULL, ndim = 2)
#'     plot(mds.ord3)
#'
#'     mds.seq.ord1 <- make_seq_mappoly(mds.ord1)
#'     mds.seq.ord2 <- make_seq_mappoly(mds.ord2)
#'     mds.seq.ord3 <- make_seq_mappoly(mds.ord3)
#'
#'
#'    }
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} mostly adapted from TetraploidSNPMap codes
#'
#' @references
#'  Preedy, K. F., & Hackett, C. A. (2016). A rapid marker ordering approach for
#'  high-density genetic linkage maps in experimental autotetraploid populations
#'  using multidimensional scaling. _Theoretical and Applied Genetics_, 129(11),
#'  2117-2132. \url{http://doi.org/10.1007/s00122-016-2761-8}
#'
#' @importFrom smacof smacofSym
#' @importFrom princurve principal.curve
#' @importFrom stats runif 
#' @importFrom utils read.csv write.csv
#' @importFrom MDSMap calc.nnfit
#' @export mds_mappoly
mds_mappoly<-function(input.mat,
                      p = NULL,
                      n = NULL,
                      ndim = 2,
                      weight.exponent = 2,
                      verbose = TRUE)
{
  o<-is.na(input.mat$rec.mat)
  input.mat$rec.mat[o]<-0.0000001
  input.mat$lod.mat[o]<-0.0000001
  if(weight.exponent != 1)
    input.mat$lod.mat<-input.mat$lod.mat^weight.exponent
  diag(input.mat$lod.mat)<-diag(input.mat$rec.mat)<-NA
  locinames <- rownames(input.mat$rec.mat)
  lodrf <- list(rf = input.mat$rec.mat, lod = input.mat$lod.mat, nloci = ncol(input.mat$rec.mat), locinames = locinames)
  confplotno<-1:lodrf$nloci
  if(!is.null(n)){
    if(!is.numeric(n))n<-which(lodrf$locinames%in%n)    
    r<-lodrf$rf[-n,-n]
    lod<-lodrf$lod[-n,-n]
    confplotno<-confplotno[-n]
  } else {
    r<-lodrf$rf
    lod<-lodrf$lod
  }
  #M <- MDSMap::dmap(r,"haldane")
  M <- imf_h(r)
  nloci=length(confplotno)
  smacofsym<-smacof::smacofSym(M,ndim=ndim,weightmat=lod,itmax=100000)
  pc1<-princurve::principal_curve(smacofsym$conf,maxit=150,spar=p,smoother="smooth_spline")
  scale<-sum(smacofsym$delta)/sum(smacofsym$dhat) 
  # Configuration dissim are based on the normalized observed diss - dhat. 
  # True observed dissimilarities are delta
  maporder<-pc1$ord
  estpos<-pc1$lambda[maporder]*scale*100
  # gives the estimated length from the beginning of the line
  rownames<-lodrf$locinames[maporder]
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder, Vectorize(function(i,j)lod[i,j]))
  rownames(distmap)<-rownames;colnames(distmap)<-rownames
  rownames(lodmap)<-rownames;colnames(lodmap)<-rownames
  if(!is.null(n))  {
    locikey<-data.frame(locus=lodrf$locinames[-n],confplotno=confplotno)
  } else {
    locikey<-data.frame(locus=lodrf$locinames,confplotno=confplotno)
  }
  nnfit<-MDSMap::calc.nnfit(distmap,lodmap,estpos)
  locimap<-data.frame(confplotno=confplotno[maporder],locus=locikey$locus[maporder],position=estpos,nnfit=nnfit$pointfits,row.names=1:nloci)
  if(!is.null(n)) {
    removedloci<-data.frame(n,lodrf$locinames[n],row.names=NULL)
  } else {
    removedloci<-n
  }
  map<-list(smacofsym=smacofsym,pc=pc1,distmap=distmap,lodmap=lodmap,locimap=locimap,length=max(estpos),removed=n,locikey=locikey,meannnfit=nnfit$meanfit)
  if(verbose)
  {
    cat(paste('Stress:', round(map$smacofsym$stress,5)))
    cat(paste('\nMean Nearest Neighbour Fit:', round(map$meannnfit,5)))
  }
  map$data.name<-input.mat$data.name
  if(ndim == 2) {
    return(structure(map, class="mappoly.pcmap"))
  } else {
    return(structure(map, class="mappoly.pcmap3d"))
  }
}

#' @rdname mds_mappoly
#' @keywords internal
#' @export
mappoly.print.pcmap<-function(x, ...)
{
  cat("\nThis is an object of class 'mappoly.mds'")
  cat("\nNumber of markers: ", nrow(x$locimap))
  cat("\nNumber of dimensions used: ", x$sm$ndim)
  cat("\nStress: ", x$smacofsym$stress)
  cat("\nMean Nearest Neighbour Fit:", x$meannnfit)
}

#' @rdname mds_mappoly
#' @keywords internal
#' @export
mappoly.print.pcmap3d<-function(x, ...)
{
  cat("\nThis is an object of class 'mappoly.mds'")
  cat("\nNumber of markers: ", nrow(x$locimap))
  cat("\nNumber of dimensions used: ", x$sm$ndim)
  cat("\nStress: ", x$smacofsym$stress)
  cat("\nMean Nearest Neighbour Fit:", x$meannnfit)
}

#' @rdname mds_mappoly
#' @keywords internal
#' @importFrom utils getFromNamespace
#' @export
mappoly.plot.pcmap<-function(x, D1lim = NULL, D2lim = NULL, displaytext = FALSE,...)
{
  fun <- getFromNamespace("plot.pcmap", "MDSMap")
  fun(x, D1lim = D1lim, D2lim = D2lim, displaytext = displaytext, ...)
}

#' @rdname mds_mappoly
#' @keywords internal
#' @export
mappoly.plot.pcmap3d<-function(x, D1lim = NULL, D2lim = NULL, displaytext = FALSE,...)
{
  fun <- getFromNamespace("plot.pcmap3d", "MDSMap")
  fun(x, D1lim = D1lim, D2lim = D2lim, displaytext = displaytext, ...)
}
