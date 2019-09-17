#' Compute homolog probabilities 
#' 
#' Compute homolog probabilities for all individuals in the full-sib
#' offspring
#'
#' @param input.probs an object of class \code{"mappoly.genoprob"} 
#' @param x an object of class \code{"mappoly.homoprob"}
#' @param stack loagical. If \code{TRUE}, the probabilities are stacked in the 
#'              plot. 
#' @param lg indicates which linkage group should be plotted. If \code{NULL} 
#'           (default), the function plots the first linkage group
#' @param ind indicates which individual should be plotted. It can be the 
#'            position of the individuals in the data set or its name. 
#'            If \code{NULL} (default), the function plots the first 
#'            individual.
#' @param use.plotly if \code{TRUE}, it uses plotly interactive graphics
#' @param ... unused argument
#' 
#'@examples
#' \dontrun{
#'   ## hexaploid example
#'   w <- lapply(maps.hexafake, calc_genoprob)
#'   h.prob <- calc_homoprob(w)
#'   print(h.prob)
#'   plot(h.prob)
#'   plot(h.prob, lg = 1, ind = 5, use.plotly = TRUE)
#'   
#'   ## tetraploid solcap example
#'   #load("~/repos/MAPpoly_vignettes/vignette_1/maps.rda")
#'   h.prob.solcap<-calc_homoprob(genoprob.err)
#'   print(h.prob.solcap)
#'   plot(h.prob.solcap, stack = TRUE, ind = "ind_10")
#'}
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M. et al. (2019) Unraveling the hexaploid 
#'     sweetpotato inheritance using ultra-dense multilocus mapping, 
#'     _submited_. \url{https://doi.org/10.1101/689638}
#'     
#' @export
#' @importFrom ggplot2 ggplot geom_density ggtitle facet_grid theme_minimal ylab xlab aes vars
#' @importFrom plotly ggplotly
#' 
calc_homoprob<-function(input.probs){
  if(class(input.probs) == "mappoly.genoprob")
    input.probs <- list(input.probs)
  if(class(input.probs) == "list"){
    if(!all(sapply(input.probs, class) == "mappoly.genoprob"))
      stop(deparse(substitute(input.probs)), " is not an object of class 'mappoly.sequence' neither a list.")
  }
  df.res <- NULL
  for(j in 1:length(input.probs)){
    cat("\nLinkage group ", j, "...")
    stt.names<-dimnames(input.probs[[j]]$probs)[[1]] ## state names
    mrk.names<-dimnames(input.probs[[j]]$probs)[[2]] ## mrk names
    ind.names<-dimnames(input.probs[[j]]$probs)[[3]] ## individual names
    v<-c(2,4,6,8,10,12)
    names(v)<-choose(v,v/2)^2
    m<-v[as.character(length(stt.names))]
    hom.prob<-array(NA, dim = c(m*2, length(mrk.names), length(ind.names)))
    dimnames(hom.prob)<-list(letters[1:(2*m)], mrk.names, ind.names)
    for(i in letters[1:(2*m)])
      hom.prob[i,,] <- apply(input.probs[[j]]$probs[grep(stt.names, pattern = i),,], c(2,3), function(x) round(sum(x, na.rm = TRUE),4))
    df.hom<-reshape::melt(hom.prob)
    map<-data.frame(map.position = input.probs[[j]]$map, marker = names(input.probs[[j]]$map))
    colnames(df.hom)<-c("homolog", "marker", "individual", "probability")
    df.hom<-merge(df.hom, map, sort = FALSE)
    df.hom$LG <- j
    df.res<-rbind(df.res, df.hom)
  }
  cat("\n")
  structure(list(homoprob = df.res, co.pos = NULL), class = "mappoly.homoprob")
}

#' @rdname calc_homoprob
#' @keywords internal
#' @export
print.mappoly.homoprob<-function(x, ...){
 head(x$homoprob, 20)
}

#' @rdname calc_homoprob
#' @keywords internal
#' @export
plot.mappoly.homoprob<-function(x, stack = FALSE, lg = NULL, ind = NULL, use.plotly = TRUE, ...){
  all.ind<-as.character(unique(x$homoprob$individual))
  if(is.null(ind)){
    ind<-as.character(all.ind[1])
    df.pr1<-subset(x$homoprob, individual == ind)  
  } else if(is.numeric(ind)) {
    if(ind > length(all.ind))
      stop("Please chose an individual number between 1 and ", length(all.ind))
    ind<-as.character(all.ind[ind])
    df.pr1<-subset(x$homoprob, individual == ind)  
  } else if (is.character(ind)){
    if(!ind%in%all.ind)
      stop("Invalid individual name")
  } else stop("Invalid individual name")
  
  LG<-individual<-map.position<-probability<-homolog<-NULL
  if(stack){
    ##subset linkage group
    if(!is.null(lg)){
      df.pr1<-subset(x$homoprob, LG == lg)
    }
    df.pr1<-subset(x$homoprob, individual == ind)  
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7, position = "stack") + 
      ggplot2::ggtitle(paste(ind, "   LG", lg)) + 
      ggplot2::facet_grid(rows = ggplot2::vars(LG)) + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  } else {
    ##subset linkage group
    if(is.null(lg)){
      lg<-1
      df.pr1<-subset(x$homoprob, LG == lg)
    } else df.pr1<-subset(x$homoprob, LG == lg)
    df.pr1<-subset(df.pr1, individual == ind)  
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7) + 
      ggplot2::ggtitle(paste(ind, "   LG", lg)) + 
      ggplot2::facet_grid(rows = ggplot2::vars(homolog)) + 
      ggplot2::theme_minimal() + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  }
  if(use.plotly)
    p <- plotly::ggplotly(p)
  return(p)
}





