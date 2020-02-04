#' Homolog probabilities 
#' 
#' Compute homolog probabilities for all individuals in the full-sib
#' population
#'
#' @param input.genoprobs an object of class \code{mappoly.genoprob}
#' 
#' @param x an object of class \code{mappoly.homoprob}
#' 
#' @param stack logical. If \code{TRUE}, probability profiles of all homologs
#'              are stacked in the plot (default = FALSE)
#'              
#' @param lg indicates which linkage group should be plotted. If \code{NULL} 
#'           (default), the function plots the first linkage group. If 
#'           \code{"all"}, the function plots all linkage groups
#'           
#' @param ind indicates which individuals should be plotted. It can be the 
#'            position of the individuals in the data set or it's name. 
#'            If \code{NULL} (default), the function plots the first 
#'            individual
#'            
#' @param use.plotly if \code{TRUE} (default), it uses plotly interactive graphics
#' 
#' @param ... unused arguments
#' 
#'@examples
#' \dontrun{
#'   ## hexaploid example
#'   w1 <- lapply(maps.hexafake, calc_genoprob)
#'   h.prob <- calc_homoprob(w1)
#'   print(h.prob)
#'   plot(h.prob)
#'   plot(h.prob, lg = 1, ind = 5, use.plotly = FALSE)
#'   plot(h.prob, lg = c(1,3), ind = 15, use.plotly = FALSE)
#'   plot(h.prob, lg = "all")
#'   
#'   ## tetraploid solcap example
#'   w2<-lapply(solcap.dose.map, calc_genoprob)
#'   h.prob.solcap<-calc_homoprob(w2)
#'   print(h.prob.solcap)
#'   plot(h.prob.solcap, ind = "ind_10")
#'   plot(h.prob.solcap, stack = TRUE, ind = 5)
#'   
#'   w3<-lapply(solcap.err.map, calc_genoprob_error, error = 0.05)
#'   h.prob.solcap.err<-calc_homoprob(w3)
#'   plot(h.prob.solcap, lg = 1, ind = 100, use.plotly = FALSE)
#'   plot(h.prob.solcap.err, lg = 1, ind = 100, use.plotly = FALSE)
#'}
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#'
#' @references
#'     Mollinari, M. et al. (2019) Unraveling the hexaploid 
#'     sweetpotato inheritance using ultra-dense multilocus mapping. 
#'     _G3: Genes, Genomes, Genetics_. \url{http://dx.doi.org/10.1534/g3.119.400620}
#'     
#' @export
#' @importFrom ggplot2 ggplot geom_density ggtitle facet_grid theme_minimal ylab xlab aes vars
#' @importFrom plotly ggplotly
#' 
calc_homoprob<-function(input.genoprobs){
  if(class(input.genoprobs) == "mappoly.genoprob")
    input.genoprobs <- list(input.genoprobs)
  if(class(input.genoprobs) == "list"){
    if(!all(sapply(input.genoprobs, class) == "mappoly.genoprob"))
      stop(deparse(substitute(input.genoprobs)), " is not an object of class 'mappoly.sequence' neither a list containing 'mappoly.sequence' objects.")
  }
  df.res <- NULL
  for(j in 1:length(input.genoprobs)){
    cat("\nLinkage group ", j, "...")
    stt.names<-dimnames(input.genoprobs[[j]]$probs)[[1]] ## state names
    mrk.names<-dimnames(input.genoprobs[[j]]$probs)[[2]] ## mrk names
    ind.names<-dimnames(input.genoprobs[[j]]$probs)[[3]] ## individual names
    v<-c(2,4,6,8,10,12)
    names(v)<-choose(v,v/2)^2
    m<-v[as.character(length(stt.names))]
    hom.prob<-array(NA, dim = c(m*2, length(mrk.names), length(ind.names)))
    dimnames(hom.prob)<-list(letters[1:(2*m)], mrk.names, ind.names)
    for(i in letters[1:(2*m)])
      hom.prob[i,,] <- apply(input.genoprobs[[j]]$probs[grep(stt.names, pattern = i),,], c(2,3), function(x) round(sum(x, na.rm = TRUE),4))
    df.hom<-reshape::melt(hom.prob)
    map<-data.frame(map.position = input.genoprobs[[j]]$map, marker = names(input.genoprobs[[j]]$map))
    colnames(df.hom)<-c("homolog", "marker", "individual", "probability")
    df.hom<-merge(df.hom, map, sort = FALSE)
    df.hom$LG <- j
    df.res<-rbind(df.res, df.hom)
  }
  cat("\n")
  structure(list(info = list(m = m, nind = length(ind.names)) , homoprob = df.res), class = "mappoly.homoprob")
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
plot.mappoly.homoprob<-function(x, stack = FALSE, lg = NULL, 
                                ind = NULL, use.plotly = TRUE, ...){
  all.ind<-as.character(unique(x$homoprob$individual))
  #### Individual handling ####
  if(length(ind) > 1){
    warning("More than one individual provided: using the first one")
    ind<-ind[1]
  }
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
  #### LG handling ####
  if(is.null(lg))
    lg <- 1
  if(lg=="all")
    lg <- unique(x$homoprob$LG)
  LG<-individual<-map.position<-probability<-homolog<-NULL
  if(length(lg) > 1 & !stack)
  {
    message("Using 'stack = TRUE' to plot multiple linkage groups")
    stack <- TRUE
  }
  if(stack){
    ##subset linkage group
    if(!is.null(lg)){
      df.pr1<-subset(x$homoprob, LG%in%lg)
      df.pr1<-subset(df.pr1, individual == ind)
    } else 
      df.pr1<-subset(x$homoprob, individual == ind)
    p <- ggplot2::ggplot(df.pr1, ggplot2::aes(x = map.position, y = probability, fill = homolog, color  = homolog)) +
      ggplot2::geom_density(stat = "identity", alpha = 0.7, position = "stack") + 
      ggplot2::ggtitle(ind) + 
      ggplot2::facet_grid(rows = ggplot2::vars(LG)) + 
      ggplot2::ylab(label = "Homologs probabilty") +
      ggplot2::xlab(label = "Map position")
  } else {
    ##subset linkage group
    if(is.null(lg)){
      lg<-1
      df.pr1<-subset(x$homoprob, LG %in% lg)
    } else df.pr1<-subset(x$homoprob, LG %in% lg)
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