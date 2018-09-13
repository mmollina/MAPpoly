#' Estimates loci position using Multidimensional Scaling
#'
#' Estimates loci position using Multidimensional Scaling proposed by
#' Preedy and Hackett (2016). The code is an adaptation from
#' TetraploidSNPMap code, available under GNU GENERAL PUBLIC LICENSE,
#' Version 3, at
#' \url{https://github.com/BiomathematicsAndStatisticsScotland/TetraploidSNPMap}
#'
#' @param input.mat an object of class \code{mappoly.input.matrix}.
#'
#' @param p integer. The smoothing parameter for the principal curve.
#'   If \code{NULL} this will be done using leave one out cross validation.
#'
#' @param n vector of integers or strings containing loci to be omitted from the analysis
#'
#' @param ndim number of dimensions to be considered in the multidimensional scaling procedure
#'
#' @param x an object of class \code{mappoly.mds}
#'
#' @param weight.exponent the exponent that should be used in the LOD score values to weight the
#'        MDS procedure.
#'
#' @param verbose if \code{TRUE},  display information about the analysis
#'
#' @param displaytext logical. If \code{TRUE}, display the name of the markers in the
#'   diagnostic plot.
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
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu} mostly adapding from TetraploidSNPMap codes
#'
#' @references
#'  Preedy, K. F., & Hackett, C. A. (2016). A rapid marker ordering approach for
#'  high-density genetic linkage maps in experimental autotetraploid populations
#'  using multidimensional scaling. Theoretical and Applied Genetics, 129(11),
#'  2117-2132. \url{http://doi.org/10.1007/s00122-016-2761-8}
#'
#' @importFrom smacof smacofSym
#' @importFrom princurve principal.curve
#' @importFrom stats runif
#' @importFrom utils read.csv write.csv
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
  locinames <- get(input.mat$data.name, pos=1)$mrk.names[as.numeric(colnames(input.mat$rec.mat))]
  X <- list(rf = input.mat$rec.mat, lod = input.mat$lod.mat, nloci = ncol(input.mat$rec.mat), locinames = locinames)
  if(!is.null(n)){
    t<-table(n)[(table(n)>1)]
    if(length(t)>0){
      write("Warning: the list of omitted loci is not unique, please remove duplicates of the loci listed below","")
      write(names(t),"")
      return(0)
    }
  }
  map<-calc_maps_load_and_trim_prin_curve(lodrf = X ,spar=p,n,ndim=ndim)
  if(verbose)
  {
    write(paste('Stress:',map$sm$stress),"")
    write(paste('Nearest Neighbour Fit:',sum(map$locimap$nnfit)),"")
    write(paste('Mean Nearest Neighbour Fit:',mean(map$locimap$nnfit)),"")
  }
  map$data.name<-input.mat$data.name
  structure(map, class="mappoly.mds")
}
#' @rdname mds_mappoly
#' @export
print.mappoly.mds<-function(x, ...)
{
  cat("\nThis is an object of class 'mappoly.mds'")
  cat("\nNumber of markers: ", nrow(x$locimap))
  cat("\nNumber of dimensions used: ", x$sm$ndim)
}
#' @rdname mds_mappoly
#' @export
plot.mappoly.mds<-function(x, displaytext = FALSE, ...)
{
  if(x$sm$ndim==2)
    plot_diag_pc(x, displaytext = displaytext)
  else
    plot_diag_pc_3d(x, displaytext = displaytext)
}

#' Get nearest informative marker
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_nearest_informative<-function(loci,lodmap){
  #split matrix by loci
  neighbours<-NULL

  if(loci>1) {
    locileft<-lodmap[loci,(loci-1):1]
    if(length(which(locileft!=0))>0)    neighbours<-loci-min(which(locileft!=0))
  }
  if(loci<dim(lodmap)[2]){
    lociright<-lodmap[loci,(loci+1):dim(lodmap)[2]]
    if(length(which(lociright!=0))>0)  neighbours<-c(neighbours,loci+min(which(lociright!=0)))
  }
  neighbours
}

#' Get nearest informative marker
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
calc_nnfit_loci<-function(loci,distmap,lodmap,estmap){
  nns<-get_nearest_informative(loci,lodmap)
  obs<-distmap[loci,nns]
  est<-estmap[loci]-estmap[nns]
  nn.fit<-sum(abs(obs-est))
  nn.fit
}

#' calc nnfit
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
calc_nnfit<-function(distmap,lodmap,estmap){
  pointfits<-unlist(lapply(1:dim(distmap)[2],calc_nnfit_loci,distmap=distmap,lodmap=lodmap,estmap=estmap))
  fit<-sum(pointfits)
  meanfit<-mean(pointfits)
  list(fit=fit,pointfits=pointfits,meanfit=meanfit)
}

#'  Estimates loci position using Principal curves
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @importFrom smacof smacofSym
#' @importFrom princurve principal.curve
#' @export
calc_maps_load_and_trim_prin_curve<-function(lodrf, spar=NULL, n=NULL, ndim=2, weightfn='lod')
  {
  confplotno<-1:lodrf$nloci
  if(!is.null(n)){
    if(!is.numeric(n))n<-which(lodrf$locinames%in%n)
    r<-lodrf$rf[-n,-n]
    lod<-lodrf$lod[-n,-n]
    confplotno<-confplotno[-n]
    locinames<-lodrf$locinames[-n]
  }
  else
  {
    r<-lodrf$rf
    lod<-lodrf$lod
    locinames<-lodrf$locinames
  }
  M<-imf_h(r) # Haldane
  nloci=length(confplotno)
  sm <- smacof::smacofSym(delta = M, ndim = ndim,
                        weightmat = lod,
                        itmax = 100000)
  pc1 <- princurve::principal.curve(x = sm$conf,
                                    maxit = 150,
                                    spar = spar)
  if(ndim==2) {
    conf<-data.frame(plotno=confplotno,
                     name=locinames,
                     MDS1=sm$conf[,1],
                     MDS2=sm$conf[,2],
                     Smooth1=pc1$s[,"D1"],
                     Smooth2=pc1$s[,"D2"])
  } else {
    if(ndim > 3) warning("Using 3 dimensions on MDS")
    conf<-data.frame(plotno=confplotno,
                     name=locinames,
                     MDS1=sm$conf[,1],
                     MDS2=sm$conf[,2],
                     MDS3=sm$conf[,3],
                     Smooth1=pc1$s[,"D1"],
                     Smooth2=pc1$s[,"D2"],
                     Smooth3=pc1$s[,"D3"])
  }
  # gives the projection of points onto the line
  # configuration dissim are basedon the normalized observed distances - dhat.
  # True observed dissimilarities are delta
  scale<-sum(sm$delta)/sum(sm$dhat)
  maporder<-pc1$tag
  estpos<-pc1$lambda[maporder]*scale#*100                #gives the estimated length from the beginning of the line
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder, Vectorize(function(i,j)lod[i,j]))

  if(!is.null(n)){
    locikey<-data.frame(locus=lodrf$locinames[-n],confplotno=confplotno)
  } else
    locikey<-data.frame(locus=lodrf$locinames,confplotno=confplotno)
  if(!is.null(n)){
    removedloci<-data.frame(n,lodrf$locinames[n],row.names=NULL)
  } else removedloci<-n
  l<-dim(sm$conf)[1]
  nnfit<-calc_nnfit(distmap,lodmap,estpos)
  locimap<-data.frame(confplotno=confplotno[maporder],
                      name=locikey$locus[maporder],
                      position=estpos,
                      nnfit=nnfit$pointfits)
  phasemap<-data.frame(locikey$locus[maporder],
                       estpos)
  #write.table(nloci,file=paste(st,'_phasemap.txt',sep=""),sep=",",row.names=FALSE,quote=FALSE,col.names=FALSE)
  #write.table(phasemap,file=paste(st,'_phasemap.txt',sep=""),sep=",",row.names=FALSE,quote=FALSE,col.names=FALSE,append=TRUE)
  #write.table(conf,file=paste(st,'_conf.txt',sep=""),sep=",",row.names=FALSE)
  #write.table(locimap,file=paste(st,'_estimatedmap.txt',sep=""),sep=",",row.names=FALSE,quote=FALSE)
  list(M=M,
       lod=lod,
       sm=sm,
       pc=pc1,
       locimap=locimap,
       length=max(estpos),
       removed=n,
       locikey=locikey,
       scale=scale,
       confplotno=confplotno,
       nnfit=nnfit)
}

#'  Plot disgnostic graphics (2D)
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
plot_diag_pc<-function(mappc,D1lim=NULL,D2lim=NULL,displaytext=TRUE){

  with(mappc,{
    if (displaytext==TRUE) labels=locikey$locus else labels=locikey$confplotno
    #png(paste(st,'_diagplot.png',sep=""),width=960,height=480)
    op<-par(mfrow=c(1,2))
    plot(sm$conf,type="n",main='MDS with principal curve',xlim=D1lim,ylim=D2lim,xlab='Dim 1',ylab='Dim 2')
    text(sm$conf,labels=labels,cex=0.8)
    lines(pc)
    if (displaytext==TRUE) labels1=locimap$name else labels1=locimap$confplotno
    plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit')
    text(locimap$position,locimap$nnfit,labels1)
    par(op)
    #dev.off()
    })
}

#'  Plot disgnostic graphics (3D)
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
plot_diag_pc_3d<-function(mappc,D1lim=NULL,D2lim=NULL,D3lim=NULL,displaytext=TRUE){
  numscale = 1.0
  with(mappc,{
    if (displaytext==TRUE) labels=locikey$locus else labels=locikey$confplotno
    #png(paste(st,'_diagplot.png',sep=""),height=1200,width=1200)
    op<-par(mfrow=c(2,2))
    plot(sm$conf[,'D1'],sm$conf[,'D2'],type="n",main='MDS with principal curve',xlab='Dimension 1',ylab='Dimension 2',xlim=D1lim,ylim=D2lim)
    text(sm$conf[,'D1'],sm$conf[,'D2'],labels=labels,cex=numscale)
    lines(pc$s[,'D1'][pc$tag],pc$s[,'D2'][pc$tag])
    plot(sm$conf[,'D1'],sm$conf[,'D3'],type="n",main='MDS with principal curve',xlab='Dimension 1',ylab='Dimension 3',xlim=D1lim,ylim=D3lim)
    text(sm$conf[,'D1'],sm$conf[,'D3'],labels=labels,cex=numscale)
    lines(pc$s[,'D1'][pc$tag],pc$s[,'D3'][pc$tag])
    plot(sm$conf[,'D2'],sm$conf[,'D3'],type="n",main='MDS with principal curve',xlab='Dimension 2',ylab='Dimension 3',xlim=D2lim,ylim=D3lim)
    text(sm$conf[,'D2'],sm$conf[,'D3'],labels=labels,cex=numscale)
    lines(pc$s[,'D2'][pc$tag],pc$s[,'D3'][pc$tag])
    if (displaytext==TRUE)
      labels1=locimap$name
    else
      labels1=locimap$confplotno
    plot(locimap$position,locimap$nnfit,type='n',xlab='Position',ylab='nnfit')
    text(locimap$position,locimap$nnfit,labels1)
    par(op)
    #dev.off()
    #write.table(sm$conf, "sm.conf")
    #write.table(locikey, "locikey")
    #write.table(pc$s[pc$tag,], "pc")
    #plot3d(sm$conf,type="n")
    #text3d(sm$conf,text=labels)
    #lines3d(pc$s[pc$tag,])
  }
  )
}

#' Calculates a new nnfit based on a new order after estimating a map.
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
recalc_nnfit_from_map<-function(estmap,mapobject,fname=NULL){
  estmap<-read.csv(estmap,header=FALSE)
  M<-mapobject$map$M
  lod<-mapobject$map$lod
  lnames<-colnames(M)
  names<-estmap[,1]
  maporder<- sapply(1:length(names),function(i)which(lnames==names[i]))
  distmap<-outer(maporder,maporder,Vectorize(function(i,j)M[i,j]))
  lodmap<-outer(maporder,maporder,Vectorize(function(i,j)lod[i,j]))
  nnfit<-calc_nnfit(distmap,lodmap,estmap[,2])
  newmap<-data.frame(name=estmap[,1],position=estmap[,2],nnfit=nnfit$pointfits)
  if(!is.null(fname)) write.csv(newmap,file=fname)
  nnfit
}

#' get_dist_loci
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_dist_loci<-function(loci,estmap,realmap){
  l<-estmap$name[loci]
  dist<-estmap[estmap$name==l,]$position-realmap[which(realmap[,1]==l),2]
  dist
}

#' Calculates the mean of the square distances of points from their real position
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
mean_dist_from_truth<-function(estmap,realmap){
  dist<-unlist(lapply(1:dim(estmap)[1],get_dist_loci,estmap=estmap,realmap=realmap))
  meansquaredist<-mean(dist^2)
  pointdist<-cbind(name=estmap$name,dist=dist)
  list(pointdist=pointdist, meansquaredist=meansquaredist)
}
